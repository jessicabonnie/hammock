#!/usr/bin/env python3
"""
Benchmark hammock-cpp (Mode B) vs bedtools jaccard.

Adapted from hammock/benchmarks/benchmark_cpp_vs_bedtools.py for the
hammock_claude refactor. Differences from the original:

  * Mode B only (A/C aren't compared to bedtools here).
  * New binary path: hammock_claude/build/*/hammock-cpp.
  * Parses the new --verbose stderr ("Sketching: X ms" / "Pairwise+write: Y ms")
    instead of the original "TIMING:" line.
  * Real memory measurement via /usr/bin/time -v -> Maximum RSS, replacing the
    original's dead psutil monitor.
  * BED files are pre-sorted before timing begins (bedtools requires sorted
    input; sorting time isn't charged to bedtools).
"""

import argparse
import concurrent.futures
import contextlib
import csv
import glob
import os
import platform
import random
import re
import resource
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from typing import Any, Callable, Dict, List, Optional

import numpy as np  # type: ignore


def _sort_one(path: str) -> None:
    sorted_path = path + ".sorted"
    with open(sorted_path, "w") as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", path], stdout=out, check=True)
    os.rename(sorted_path, path)

# ---------- paths ----------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
BEDTOOLS_SCRIPT = os.path.join(SCRIPT_DIR, "bedtools.sh")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")


MIN_BINARY_VERSION = (0, 7, 0)


def check_binary_version(binary: str) -> str:
    """Refuse a binary older than the flags and timers this harness assumes.

    Call this on the *resolved* path, whatever produced it. The stale-binary
    case that actually bites is `pip install -e .` skipped after a rebuild: the
    build tree updates and the site-packages copy that HAMMOCK_CPP_BIN points at
    does not, so a probe that only covered the build-tree glob would check the
    one path where nothing can go wrong.
    """
    proc = subprocess.run([binary, "--version"], capture_output=True, text=True)
    text = (proc.stdout or "").strip()
    if proc.returncode != 0 or not text.startswith("hammock-cpp "):
        raise RuntimeError(
            f"{binary} does not understand --version, so it predates 0.7.0. "
            "This harness passes --register-equality (or --no-metrics on "
            "older binaries) and parses microsecond timings, neither of "
            "which that binary has. Rebuild AND reinstall: "
            "pip install -e . --no-build-isolation")
    got = text.split()[1]
    parts = tuple(int(x) for x in got.split(".")[:3])
    if parts < MIN_BINARY_VERSION:
        raise RuntimeError(
            f"{binary} is hammock-cpp {got}; this harness needs "
            f">= {'.'.join(map(str, MIN_BINARY_VERSION))}.")
    return got


_HELP_CACHE: Dict[str, str] = {}


def _probe_binary_help(binary: str) -> str:
    """Cached `<binary> --help` text.

    Used by `_metrics_off_flag` to check flag support directly rather than
    trust `--version`, which can lag behind the actual CLI contract (see
    that function's docstring for why this matters specifically for the
    metrics-column restructure). Cached per binary PATH, not content --
    `fusion_ab.py`'s `--pre-binary`/`--post-binary` are different paths by
    construction, so each gets its own independent cache entry no matter how
    many times `run_hammock()` is called against them.

    Raises RuntimeError on a bad/non-executable path instead of letting a
    raw OSError/PermissionError escape -- this probe did not exist before
    the metrics-column restructure (docs/seed-metrics-column-restructure.md),
    so a bad `--pre-binary`/`--post-binary` value surfaces here for the
    first time now.
    """
    if binary not in _HELP_CACHE:
        try:
            proc = subprocess.run([binary, "--help"], capture_output=True, text=True)
        except OSError as e:
            raise RuntimeError(f"could not run {binary} --help: {e}") from e
        _HELP_CACHE[binary] = (proc.stdout or "") + (proc.stderr or "")
    return _HELP_CACHE[binary]


def _metrics_off_flag(binary: str) -> str:
    """Flag string for the "cheap/register-equality" side of the
    --metrics/(--no-metrics | --register-equality) pair, chosen per BINARY,
    not per caller -- `fusion_ab.py` drives two binaries of different
    vintages through the same `run_hammock()` call in one process, so the
    choice cannot be a single hardcoded string or module-level constant.

    Checked by CAPABILITY (does --help mention --register-equality), not by
    --version number: the metrics-column restructure
    (docs/seed-metrics-column-restructure.md) removed --no-metrics in favour
    of --register-equality/--re without pyproject.toml's version bump
    landing yet (that's Part 7's job) -- a binary built off this exact
    worktree already lacks --no-metrics while still self-reporting
    `--version 0.7.1`, the identical string a genuinely pre-restructure
    binary reports. A version-number gate cannot tell those two apart;
    checking --help for the literal flag can, and stays correct regardless
    of when the version bump actually lands.
    """
    return ("--register-equality" if "--register-equality" in _probe_binary_help(binary)
            else "--no-metrics")


def find_hammock_cpp() -> str:
    """Locate the hammock-cpp binary built under build/<plat>/hammock-cpp."""
    env = os.environ.get("HAMMOCK_CPP_BIN")
    if env and os.path.exists(env):
        return env
    candidates = glob.glob(os.path.join(REPO_ROOT, "build", "*", "hammock-cpp"))
    if candidates:
        return max(candidates, key=os.path.getmtime)
    raise FileNotFoundError(
        "hammock-cpp not found. Build the project first (pip install -e . --no-build-isolation) "
        "or set HAMMOCK_CPP_BIN."
    )


# ---------- defaults ----------
NUM_INTERVALS_PER_FILE = 10000
NUM_FILES_LIST = [2, 4, 8, 16, 32, 64, 128, 256, 512]
NUM_RUNS = 3

TIME_CMD = "/usr/bin/time"

# hammock-cpp >= 0.7.0 reports microseconds and decomposes the pairwise phase.
# PAIR_RE keeps matching "Pairwise+write", whose meaning is unchanged, so
# comparison_time stays comparable with every archived sweep. The colon in
# PAIRONLY_RE/WRITE_RE is load-bearing: "^Write" alone also matches the
# "Wrote <path>" line that follows.
SKETCH_RE = re.compile(r"^Sketching:\s+(\d+)\s+us")
PAIR_RE = re.compile(r"^Pairwise\+write:\s+(\d+)\s+us")
PAIRONLY_RE = re.compile(r"^Pairwise:\s+(\d+)\s+us")
WRITE_RE = re.compile(r"^Write:\s+(\d+)\s+us")
MAXRSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")


def get_system_info() -> Dict[str, Any]:
    """Where a measurement was taken, in enough detail to judge comparability.

    `cpu_model` and `hostname` are load-bearing, not decoration: CMakeLists
    bakes in `-march=native`, so a timing taken on one CPU model is not
    strictly comparable to one taken on another, and this cluster mixes node
    types. An archived CSV that records neither cannot be checked against --
    which is exactly the gap this closes.

    `slurm_job_id` distinguishes a run inside an allocation (cores actually
    reserved) from one on a shared login/dev node, where a co-tenant job
    silently inflates wall times. "none" means no allocation.
    """
    info = {
        "hostname": platform.node(),
        # Two different numbers, and the affinity one is what a timing depends
        # on. os.cpu_count() reports the node; inside --cpus-per-task=16 it still
        # says 48, which would contradict the job's own log line. Note this also
        # means a bare cpu_count cannot distinguish "dev node, no allocation"
        # from "48-core node, allocated" -- don't use it as evidence of either.
        "cpu_count": len(os.sched_getaffinity(0)),
        "cpu_count_node": os.cpu_count(),
        "cpu_model": "unknown",
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "none"),
    }
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    info["cpu_model"] = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    # The version string cannot identify the code. pyproject's version is bumped
    # on release, not per commit, so a pre- and post-change binary both report
    # 0.7.0 -- which is exactly the pair a re-run exists to tell apart. The git
    # SHA can. Recorded with a -dirty marker because an uncommitted working tree
    # is a real state for a benchmark to be run from.
    info["git_sha"] = "unknown"
    try:
        repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        sha = subprocess.run(["git", "-C", repo, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        if sha.returncode == 0:
            dirty = subprocess.run(["git", "-C", repo, "status", "--porcelain"],
                                   capture_output=True, text=True, timeout=10)
            suffix = "-dirty" if (dirty.returncode == 0 and dirty.stdout.strip()) else ""
            info["git_sha"] = sha.stdout.strip() + suffix
    except (OSError, subprocess.SubprocessError):
        pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    info["memory_total_gb"] = int(line.split()[1]) / (1024 * 1024)
                    break
    except OSError:
        pass
    return info


def derive_seed(base: Optional[int], *parts: int) -> Optional[int]:
    """Per-file seed from a run-level base, or None to keep the legacy behaviour.

    Mixing the coordinates in (rather than seeding the global RNG once) makes a
    corpus reproducible independently of generation order, so parallelising the
    generation loop later cannot silently change the data.
    """
    if base is None:
        return None
    h = base & 0xFFFFFFFF
    for x in parts:
        h = (h * 1_000_003 + x) & 0xFFFFFFFF
    return h


def generate_bed_file(num_intervals: int, output_file: str,
                      seed: Optional[int] = None) -> None:
    """Random BED. `seed` makes it reproducible; None keeps the global RNG.

    Unseeded is the historical behaviour and stays the default, because every
    archived CSV was produced that way. It is also the reason a cross-run
    comparison of these benchmarks is only good to a few percent on the bedtools
    leg: hammock's cost barely depends on the interval content (it hashes a fixed
    number of positions) while bedtools' does, so a fresh draw moves the two legs
    by different amounts and the *ratio* drifts without either tool changing.
    """
    rng = random.Random(seed) if seed is not None else random
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    with open(output_file, "w") as f:
        for _ in range(num_intervals):
            chrom = rng.choice(chroms)
            start = rng.randint(0, 10_000_000)
            end = rng.randint(start + 100, start + 10_000)
            f.write(f"{chrom}\t{start}\t{end}\n")


def run_with_time(cmd: List[str], env: Dict[str, str] = None) -> Dict[str, Any]:
    """Run cmd under /usr/bin/time -v, capturing wall, child CPU, and max RSS.

    /usr/bin/time -v writes its report to stderr; we tee that into a temp file
    so we can parse maxrss while still seeing the child's own stderr.
    """
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".time", delete=False) as tf:
        time_log = tf.name
    try:
        wall_start = time.time()
        ru_start = resource.getrusage(resource.RUSAGE_CHILDREN)
        wrapped = [TIME_CMD, "-v", "-o", time_log] + cmd
        result = subprocess.run(wrapped, capture_output=True, text=True, check=True,
                                env=env)
        ru_end = resource.getrusage(resource.RUSAGE_CHILDREN)
        wall_end = time.time()

        cpu_time = (ru_end.ru_utime - ru_start.ru_utime) + (ru_end.ru_stime - ru_start.ru_stime)

        maxrss_kb = None
        with open(time_log) as f:
            for line in f:
                m = MAXRSS_RE.search(line)
                if m:
                    maxrss_kb = int(m.group(1))
                    break
        return {
            "wall_time": wall_end - wall_start,
            "cpu_time": cpu_time,
            "max_rss_mb": (maxrss_kb / 1024) if maxrss_kb is not None else None,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "returncode": result.returncode,
        }
    finally:
        try:
            os.remove(time_log)
        except OSError:
            pass


def run_bedtools(file1_list: str, file2_list: str, num_threads: int) -> Dict[str, Any]:
    """Time one full pairwise bedtools sweep.

    The environment below is what keeps the small-N cells meaningful. bedtools.sh
    used to run `ml bedtools/2.30.0` (~0.28 s), `ml parallel` and two
    `bedtools --version` execs inside the timed region, against a fixed floor of
    ~0.83 s measured end to end. At N=512 that is 0.04% of the cell; at N=2 it is
    99% of it, and it inverted the comparison -- Panel A reported hammock 3.40x
    faster at N=2 when bedtools, given the same work without our Lmod, is roughly
    5x faster. Resolving the binary once per process and handing bedtools.sh the
    path removes the module load from every timed call.

    This does not touch the large-N numbers in any meaningful way; it is the
    small-N end, and the location of the crossover, that it makes honest.
    """
    env = dict(os.environ)
    binary = _resolve_bedtools()
    if binary:
        env["HAMMOCK_BEDTOOLS_BIN"] = binary
        if _BEDTOOLS_LDPATH:
            env["LD_LIBRARY_PATH"] = _BEDTOOLS_LDPATH
    env["HAMMOCK_BEDTOOLS_QUIET"] = "1"
    return run_with_time(
        ["bash", BEDTOOLS_SCRIPT, file1_list, file2_list, str(num_threads)], env=env)


_BEDTOOLS_BIN: Any = "unset"
_BEDTOOLS_LDPATH: Any = None


def _resolve_bedtools():
    """Path to the same bedtools bedtools.sh will use, resolved once per process.

    Cached because resolving it costs a login shell plus an Lmod load (~0.8 s),
    which is 60x the thing being measured -- doing that per rep is what made the
    first version of the calibration meaningless.
    """
    global _BEDTOOLS_BIN, _BEDTOOLS_LDPATH
    if _BEDTOOLS_BIN != "unset":
        return _BEDTOOLS_BIN
    _BEDTOOLS_BIN = None
    module = os.environ.get("HAMMOCK_BEDTOOLS_MODULE", "bedtools/2.30.0")
    probe = (f"ml {module} 2>/dev/null; command -v bedtools; echo '---'; "
             f"ldd $(command -v bedtools) 2>/dev/null") if module else \
            "command -v bedtools; echo '---'; ldd $(command -v bedtools) 2>/dev/null"
    try:
        r = subprocess.run(["bash", "-lc", probe], capture_output=True, text=True, timeout=120)
        out = (r.stdout or "")
        head, _, tail = out.partition("---")
        path = head.strip().splitlines()[-1] if head.strip() else ""
        if path and os.path.exists(path):
            _BEDTOOLS_BIN = path
            # MINIMAL library path, derived from what ldd actually resolved.
            #
            # This is a correctness fix AND the single largest timing fix in this
            # harness. The bedtools module exports an LD_LIBRARY_PATH of 17
            # directories, of which bedtools uses 4; the other 13 are searched
            # fruitlessly by the dynamic linker on EVERY exec, and they live on
            # GPFS. A pairwise workflow is N^2 execs, so that tax is multiplied
            # by 262,144 at N=512: measured 9.13 s vs 3.80 s on 1024 pairs, and
            # 1978 s vs 716 s at N=512. It inflated bedtools by 2.4-2.8x.
            #
            # Passing only the resolved directories is faster, produces
            # byte-identical jaccards (gated below in the benchmark), and unlike
            # relying on the ambient environment it cannot fail with
            # "libbz2.so.1.0: cannot open shared object file" on a node whose
            # login environment happens not to supply a compatible libbz2.
            dirs = []
            for line in tail.splitlines():
                parts = line.split("=>")
                if len(parts) == 2:
                    lib = parts[1].strip().split(" ")[0]
                    if lib.startswith("/"):
                        d = os.path.dirname(lib)
                        if d not in dirs:
                            dirs.append(d)
            _BEDTOOLS_LDPATH = ":".join(dirs) if dirs else None
    except (OSError, subprocess.SubprocessError, IndexError):
        pass
    return _BEDTOOLS_BIN


def harness_floor_ms(binary: str, precision: int, num_threads: int, reps: int = 5):
    """Zero-work cost of each arm: what a run costs when there is nothing to compute.

    WHY THIS EXISTS. Panel A once reported hammock 3.40x faster at N=2 when
    bedtools was in fact faster, because bedtools.sh ran `ml bedtools` (0.28 s),
    `ml parallel` and two `--version` execs INSIDE the timed region while the
    hammock arm was a bare exec of the binary. A fixed ~0.83 s charged to one arm
    only is 0.04% of the N=512 cell and 99% of the N=2 cell, so it inverted the
    small-N comparison while leaving the headline untouched -- the worst possible
    shape for a bug, because the number everyone checks stays right.

    The tell was in the data the whole time and was not read: the bedtools curve
    DROPPED from N=2 to N=4. Cost cannot fall while work quadruples, so a flat
    floor was the only explanation. Measuring the floor directly turns that from
    something a reader has to notice into a number in the CSV.

    Returns (hammock_floor_ms, bedtools_floor_ms), each the median of `reps` runs
    on a 1x1 list of a single-interval BED -- as close to zero work as the
    interfaces allow. Compare each against the cell it is a component of; when
    the floor is a large fraction of a cell, that cell is measuring startup and
    no speedup should be quoted from it.
    """
    tmpdir = tempfile.mkdtemp(prefix="floor_")
    try:
        bed = os.path.join(tmpdir, "one.bed")
        with open(bed, "w") as f:
            f.write("chr1\t100\t150\n")
        lst = os.path.join(tmpdir, "one.txt")
        with open(lst, "w") as f:
            f.write(bed + "\n")
        out = os.path.join(tmpdir, "o")

        def med(fn):
            xs = []
            for _ in range(reps + 1):          # +1 warm-up, discarded
                t0 = time.time()
                try:
                    fn()
                except Exception:
                    return None
                xs.append((time.time() - t0) * 1000.0)
            if not xs[1:]:
                return None      # never let np.median([]) hand back a silent nan
            return float(np.median(xs[1:]))

        hm = med(lambda: subprocess.run(
            [binary, lst, lst, "--mode", "B", "-p", str(precision), "-o", out,
             "--threads", str(num_threads), _metrics_off_flag(binary)],
            capture_output=True, check=True))

        env = dict(os.environ)
        b = _resolve_bedtools()
        if b:
            env["HAMMOCK_BEDTOOLS_BIN"] = b
            if _BEDTOOLS_LDPATH:
                env["LD_LIBRARY_PATH"] = _BEDTOOLS_LDPATH
        env["HAMMOCK_BEDTOOLS_QUIET"] = "1"
        bt = med(lambda: subprocess.run(
            ["bash", BEDTOOLS_SCRIPT, lst, lst, str(num_threads)],
            capture_output=True, check=True, env=env))
        return hm, bt
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def bedtools_serial_ms(file1_list: str, file2_list: str, reps: int = 5):
    """Median wall time of ONE `bedtools jaccard` on the first pair, in ms.

    This is the calibration that makes the bedtools leg auditable. `bedtools
    jaccard` has no batch mode, so a pairwise workflow launches one process per
    pair -- N^2 of them -- and on these nodes process creation is capped near
    123 exec/s and does not scale with cores. The consequence is that
    "bedtools at t=16" can silently mean "bedtools at t~1.5": measured, 1024
    pairs cost the same at --jobs 16 as at --jobs 1.

    Recording a serial per-pair cost lets the CSV carry bedtools' ACHIEVED
    parallel efficiency

        eff = n_pairs * serial_per_pair / (wall * threads)

    instead of leaving a reader to assume it was near 1. Without it, a speedup
    computed against a wrapper that failed to parallelize is indistinguishable
    from one against a baseline that used its cores -- and the first inflates
    our result by up to ~6x.

    Costs `reps` extra invocations per cell (~50-100 ms), which is noise next to
    the N^2 loop it calibrates. Returns None if it cannot be measured, so a
    failure here degrades the audit rather than the benchmark.
    """
    binary = _resolve_bedtools()
    if binary is None:
        return None
    try:
        with open(file1_list) as f:
            a = f.readline().strip()
        with open(file2_list) as f:
            b = f.readline().strip()
        if not a or not b:
            return None
        # Invoke the resolved binary DIRECTLY. An earlier version ran
        # `bash -lc 'ml bedtools/...; bedtools jaccard ...'` per rep and
        # measured 810 ms/pair -- that is a login shell plus an Lmod module
        # load, not bedtools, and it inflated the implied efficiency to a
        # physically impossible 4.02x-of-4. The module is resolved once, in
        # _resolve_bedtools, and never again.
        cmd = [binary, "jaccard", "-a", a, "-b", b]
        # One warm-up: the first exec pays page-in for the binary and its
        # shared libs, a cost the N^2 loop pays once rather than per pair.
        subprocess.run(cmd, capture_output=True, timeout=120)
        times = []
        for _ in range(reps):
            t0 = time.time()
            r = subprocess.run(cmd, capture_output=True, timeout=120)
            if r.returncode != 0:
                return None
            times.append((time.time() - t0) * 1000.0)
        times.sort()
        return times[len(times) // 2]
    except Exception:
        return None


def run_hammock(
    binary: str,
    file1_list: str,
    file2_list: str,
    precision: int,
    num_threads: int,
    keep_output: bool = False,
    sub_b: float = 1.0,
    sub_b_method: str = "mixed-stride",
    metrics: bool = False,
    ie_only: bool = False,
) -> Dict[str, Any]:
    """Run hammock-cpp Mode B and return timing/RSS.

    With keep_output=True, the output CSV is left in place and its path is
    returned as result["output_csv"]; the caller is responsible for deleting it.

    sub_b < 1.0 enables point subsampling using sub_b_method (default
    mixed-stride, matching the post-9778ef8 binary default). At sub_b == 1.0
    we omit the flags to keep the cmd line byte-identical to pre-subB runs.

    metrics/ie_only together select the output shape (three-way, matching
    docs/seed-metrics-column-restructure.md's contract; ie_only takes
    precedence if both happen to be True):
      - metrics=True:  passes --metrics (the full 8-column block).
      - ie_only=True:  passes NO flag (the bare default) -- 1 column,
        jaccard_similarity_ie alone. Added 2026-08-11 (Part 9): this is now
        the cheap-ish way to obtain jaccard_similarity_ie, superseding the
        old "--metrics is the only way to get it" state of affairs. It still
        pays the fused union pass (same compute cost as --metrics), so it is
        NOT interchangeable with the register-equality arm below -- it is
        cheaper than --metrics only in write cost (fewer columns), not in
        compute. Default False so every pre-existing call site is unaffected.
      - both False: passes whatever `_metrics_off_flag(binary)` resolves to
        for THIS binary -- --register-equality on a binary built after the
        metrics-column restructure, --no-metrics on an older one (see that
        function's docstring; this is not a version-number check, since the
        two vintages can self-report the identical --version string). This
        is the genuinely cheap arm: it skips the union pass entirely.

    The metrics/ie_only block costs a union + cardinality per pair (same cost
    family), so a run with either on is NOT timing-comparable to the
    published register-equality numbers in RESULTS.md.

    Requires hammock-cpp >= 0.7.0; an older binary rejects the disable-flag
    with exit 2 rather than silently mistiming, which is the intended failure.

    Part 9 (docs/seed-metrics-column-restructure.md), 2026-08-11: --metrics-arm/
    --metrics-all (see arms_for/ie_tool_name_for_subb) used to pass
    metrics=True to obtain jaccard_similarity_ie, because --metrics was the
    only way to get that column. Retargeted to ie_only=True instead -- see
    arms_for's docstring and the small paired measurement recorded in
    experiments/bedtools_benchmark/RESULTS.md / docs/bedtools-baseline-retraction.md.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        out_prefix = tmp.name
    try:
        cmd = [
            binary,
            file1_list,
            file2_list,
            "--mode", "B",
            "-p", str(precision),
            "-o", out_prefix,
            "--threads", str(num_threads),
            "--verbose",
        ]
        if sub_b < 1.0:
            cmd += ["--subB", f"{sub_b:g}", "--subB-method", sub_b_method]
        # Explicit for --metrics/register-equality; ie_only stays flagless on
        # purpose. Note this is a v0.8.0+ (post metrics-restructure) comment:
        # omitting every flag used to mean "the always-full block" on 0.7.x
        # binaries (hence the old warning here that a bare invocation was
        # secretly expensive), but since v0.8.0 the bare default is the
        # 1-column jaccard_similarity_ie shape -- cheaper to WRITE than
        # --metrics, though not cheaper to COMPUTE (still pays the union
        # pass; only the register-equality flag skips that). ie_only takes
        # precedence over metrics if a caller somehow set both. The
        # disable-flag spelling for the register-equality arm depends on
        # THIS binary's own vintage, not on the caller's -- see
        # _metrics_off_flag's docstring.
        if ie_only:
            pass  # bare default: no flag
        elif metrics:
            cmd += ["--metrics"]
        else:
            cmd += [_metrics_off_flag(binary)]
        r = run_with_time(cmd)

        sketch_s: Optional[float] = None
        pair_s: Optional[float] = None
        pair_only_s: Optional[float] = None
        write_s: Optional[float] = None
        for line in r["stderr"].splitlines():
            m = SKETCH_RE.match(line)
            if m:
                sketch_s = int(m.group(1)) / 1e6
                continue
            m = PAIR_RE.match(line)
            if m:
                pair_s = int(m.group(1)) / 1e6
                continue
            m = PAIRONLY_RE.match(line)
            if m:
                pair_only_s = int(m.group(1)) / 1e6
                continue
            m = WRITE_RE.match(line)
            if m:
                write_s = int(m.group(1)) / 1e6
        # Fail loudly. A silent None here becomes a blank cell in the CSV via
        # aggregate(), which reads downstream as "this run had no timing" rather
        # than "the harness could not parse the binary's output" -- most likely
        # cause being a pre-0.7.0 binary still emitting milliseconds.
        if sketch_s is None or pair_s is None:
            raise RuntimeError(
                "could not parse hammock-cpp timing lines (need >= 0.7.0, which "
                f"reports microseconds). stderr was:\n{r['stderr']}")
        r["sketch_creation_time"] = sketch_s
        r["comparison_time"] = pair_s
        r["pair_time"] = pair_only_s
        r["write_time"] = write_s

        if keep_output:
            # out_prefix comes from tempfile.NamedTemporaryFile (line ~554),
            # so it is unique per call -- at most one CSV can ever match this
            # exact prefix, regardless of which tag (_ie/_re/_full, or none
            # on an archival pre-restructure binary) hammock actually wrote.
            # That's what makes an unqualified glob safe HERE, unlike the
            # shared-directory glob[0] bugs fixed elsewhere in this repo
            # (estimator_ie_crossref.py/estimator_ie_topology.py/
            # estimator_ie_tissue.py, which glob a directory many separate
            # runs write into). sorted() only for determinism if the glob
            # library ever returns more than the expected single match.
            candidates = sorted(f for f in glob.glob(out_prefix + "*") if f.endswith(".csv"))
            r["output_csv"] = candidates[0] if candidates else None
            r["_out_prefix"] = out_prefix  # caller cleans up
        return r
    finally:
        if not keep_output:
            for f in glob.glob(out_prefix + "*"):
                try:
                    os.remove(f)
                except OSError:
                    pass


def aggregate(runs: List[Dict[str, Any]], keys: List[str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k in keys:
        vals = [r[k] for r in runs if r.get(k) is not None]
        if not vals:
            out[f"mean_{k}"] = None
            out[f"std_{k}"] = None
            out[f"min_{k}"] = None
            out[f"max_{k}"] = None
            continue
        out[f"mean_{k}"] = float(np.mean(vals))
        out[f"std_{k}"] = float(np.std(vals)) if len(vals) > 1 else 0.0
        out[f"min_{k}"] = float(np.min(vals))
        out[f"max_{k}"] = float(np.max(vals))
    return out


def tool_name_for_subb(sub_b: float) -> str:
    """Tool-column identifier for a hammock run at a given subB.

    subB == 1.0 stays as "hammock_cpp_B" so legacy CSVs are byte-identical;
    anything else gets a "_subB{val}" suffix.
    """
    if sub_b >= 1.0:
        return "hammock_cpp_B"
    return f"hammock_cpp_B_subB{sub_b:g}"


# Deliberately NOT "hammock_cpp_B_metrics". Three R consumers filter
# grepl("^hammock_cpp_B", tool) and then key on sub_b rather than tool
# (plot_pairwise_scaling.R, docs/scripts/synthetic_nscaling.R, make_graphs.R), so
# a second arm sharing subB=1.0 under that prefix would silently double rows
# through their joins. A label that fails the prefix makes all three immune with
# no R changes at all.
IE_ARM_TOOL = "hammock_ie_B"

# Sentinel arm label for the bedtools leg. It is not a CSV tool name (bedtools
# columns are written from `entry["bedtools"]`, a separate key); it exists only
# so bedtools can ride the same _rotate as the hammock arms.
BEDTOOLS_TOOL = "bedtools"


def ie_tool_name_for_subb(sub_b: float) -> str:
    """Tool-column identifier for an IE-shape (bare/no-flag default,
    jaccard_similarity_ie) hammock run at a given subB, for the --metrics-all
    sweep mode (every subB arm run with that shape, not just a single extra
    subB=1.0 arm). Retargeted from --metrics to the bare default 2026-08-11
    (Part 9) -- the name keeps saying "metrics"/"ie" since it's still the
    arm that selects the IE column, just via a cheaper flag now.

    subB == 1.0 reuses IE_ARM_TOOL exactly ("hammock_ie_B") -- byte-compatible
    with every archived CSV and every consumer that does `tool ==
    "hammock_ie_B"` (plot_pairwise_scaling.R:146, plot_pairwise_scaling_supplement.R:169,232).
    subB < 1.0 gets a distinct "hammock_ie_B_subB<val>" label: it fails those
    exact-match filters (so it can't silently fan out a join that assumes one
    hammock_ie_B row per (num_files, threads, precision), the guard at
    plot_pairwise_scaling_supplement.R:185) and fails the
    `grepl("^hammock_cpp_B", tool)` prefix the no-metrics consumers use, for
    the same reason IE_ARM_TOOL does.
    """
    if sub_b >= 1.0:
        return IE_ARM_TOOL
    return f"hammock_ie_B_subB{sub_b:g}"


def arms_for(sub_b_list: List[float], metrics_arm: bool, metrics_all: bool = False):
    """(tool_label, sub_b, use_metrics) per hammock arm in a run.

    Replaces sub_b_list as the anchor for "which hammock runs happened": the
    metrics arm shares subB=1.0 with the baseline arm, so keying downstream
    consumers on sub_b alone can no longer distinguish them.

    metrics_all takes precedence over metrics_arm: every subB value in
    sub_b_list gets its own --metrics arm (ie_tool_name_for_subb), instead of
    the default no-metrics arm plus one extra fixed subB=1.0 metrics arm. Use
    this when the figure being built needs jaccard_similarity_ie wall times
    at subB < 1.0, which metrics_arm alone cannot produce.
    """
    if metrics_all:
        return [(ie_tool_name_for_subb(s), s, True) for s in sub_b_list]
    arms = [(tool_name_for_subb(s), s, False) for s in sub_b_list]
    if metrics_arm:
        arms.append((IE_ARM_TOOL, 1.0, True))
    return arms


def arms_of(entry: Dict[str, Any]):
    """(tool_label, sub_b) pairs for one result entry, for report/plot/CSV.

    Falls back to reconstructing from sub_b_list so entries produced before the
    metrics arm existed still render.
    """
    got = entry.get("hammock_arms")
    if got:
        return [(label, sub_b) for label, sub_b in got]
    return [(tool_name_for_subb(s), s) for s in entry.get("sub_b_list", [1.0])]


def arm_legend(label: str, sub_b: float) -> str:
    if label == IE_ARM_TOOL or label.startswith("hammock_ie_B_subB"):
        return f"hammock-cpp subB={sub_b:g} +IE"
    return f"hammock-cpp subB={sub_b:g}"


def _rotate(seq, k: int):
    """Rotate arm order by run index so arm identity is not confounded with
    position: the last arm in a run has just followed every other arm's cache
    and thermal state."""
    if not seq:
        return seq
    k %= len(seq)
    return seq[k:] + seq[:k]


def run_benchmark(
    binary: str,
    num_files_list: List[int],
    num_runs: int,
    num_intervals: int,
    num_threads: int,
    precision: int,
    sub_b_list: List[float],
    metrics_arm: bool = False,
    metrics_all: bool = False,
    corpus_seed: Optional[int] = None,
    with_bedtools: bool = True,
    checkpoint_cb: Optional[Callable[[List[Dict[str, Any]]], None]] = None,
) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    print("\nBenchmark configuration:")
    print(f"  hammock-cpp binary: {binary}")
    print(f"  intervals/file:     {num_intervals}")
    print(f"  corpus seed:        {corpus_seed if corpus_seed is not None else 'unseeded (not reproducible)'}")
    print(f"  file counts:        {num_files_list}")
    print(f"  runs/config:        {num_runs}")
    print(f"  threads:            {num_threads}")
    print(f"  HLL precision:      {precision}")
    print(f"  subB values:        {sub_b_list}")
    print(f"  metrics arm:        {metrics_arm}")
    print(f"  metrics all arms:   {metrics_all}")
    print(f"  system:             {get_system_info()}")

    metric_keys = ["wall_time", "cpu_time", "max_rss_mb", "sort_time"]
    # bedtools-only: serial per-pair cost and the parallel efficiency it
    # implies. Blank on hammock rows, which are one process and whose
    # threading is a different object entirely (OpenMP inside one address
    # space, not N^2 process launches under a wrapper).
    bedtools_keys = metric_keys + ["bedtools_serial_ms", "bedtools_parallel_eff"]
    # pair_time/write_time decompose comparison_time and are hammock-only, so
    # bedtools rows legitimately leave those cells blank.
    hammock_keys = metric_keys + ["sketch_creation_time", "comparison_time",
                                  "pair_time", "write_time"]

    for num_files in num_files_list:
        print(f"\n{'=' * 60}\n{num_files} files × {num_intervals} intervals\n{'=' * 60}")
        bedtools_runs: List[Dict[str, Any]] = []
        arms = arms_for(sub_b_list, metrics_arm, metrics_all=metrics_all)
        runs_by_tool: Dict[str, List[Dict[str, Any]]] = {label: [] for label, _, _ in arms}

        for run_i in range(num_runs):
            print(f"\nRun {run_i + 1}/{num_runs}: generating {2 * num_files} BED files...")
            with tempfile.TemporaryDirectory() as tmp_dir:
                file1_list, file2_list = [], []
                for i in range(num_files):
                    a = os.path.join(tmp_dir, f"set1_{i}.bed")
                    b = os.path.join(tmp_dir, f"set2_{i}.bed")
                    # run_i is mixed in deliberately: replicates must keep
                    # drawing DIFFERENT corpora or std_wall_time would collapse
                    # to machine noise alone and silently change meaning. The
                    # sequence of draws is what becomes reproducible.
                    generate_bed_file(num_intervals, a,
                                      derive_seed(corpus_seed, run_i, i, 0))
                    generate_bed_file(num_intervals, b,
                                      derive_seed(corpus_seed, run_i, i, 1))
                    file1_list.append(a)
                    file2_list.append(b)

                # Pre-sort outside the per-tool timing (bedtools requires sorted input;
                # hammock is indifferent). We measure sort wall time separately so the
                # "workflow including sort" comparison can be reconstructed downstream.
                # Sort is parallelized across num_threads workers — see
                # docs/bedtools-parallelism-caveat.md for the fairness framing.
                sort_start = time.time()
                all_paths = file1_list + file2_list
                sort_workers = max(1, min(num_threads, len(all_paths)))
                with concurrent.futures.ThreadPoolExecutor(max_workers=sort_workers) as ex:
                    for _ in ex.map(_sort_one, all_paths):
                        pass
                sort_time = time.time() - sort_start

                file1_list_path = os.path.join(tmp_dir, "file1_list.txt")
                file2_list_path = os.path.join(tmp_dir, "file2_list.txt")
                with open(file1_list_path, "w") as f:
                    f.write("\n".join(file1_list) + "\n")
                with open(file2_list_path, "w") as f:
                    f.write("\n".join(file2_list) + "\n")

                print(f"  sort:    {sort_time:.2f}s wall (parallel, {sort_workers} workers)")

                # bedtools rotates with the hammock arms rather than running
                # first every time. It used to sit permanently in position 0,
                # immediately after the pre-sort had just walked every input
                # file and left it in page cache -- so bedtools alone was
                # measured against a warm cache while the hammock arms followed
                # it. That is the same confound _rotate already exists to
                # remove; there was no reason to exempt the baseline from it.
                # (Symptom in the archived data: bedtools at N=2 reads
                # 0.772 +/- 0.475 s, CV 61.5%, and is *slower* than at N=4.)
                leg = ([(BEDTOOLS_TOOL, None, None)] if with_bedtools else []) + arms
                for tool, sub_b, use_metrics in _rotate(leg, run_i):
                    if tool == BEDTOOLS_TOOL:
                        print("  bedtools...", end=" ", flush=True)
                        r = run_bedtools(file1_list_path, file2_list_path, num_threads)
                        r["sort_time"] = sort_time
                        # Calibrate the baseline so its achieved parallelism is
                        # recorded rather than assumed -- see bedtools_serial_ms.
                        ser = bedtools_serial_ms(file1_list_path, file2_list_path)
                        r["bedtools_serial_ms"] = ser
                        n_pairs = num_files * num_files
                        r["bedtools_parallel_eff"] = (
                            (n_pairs * ser / 1000.0) / (r["wall_time"] * num_threads)
                            if ser and r.get("wall_time") else None)
                        r["run_index"] = run_i
                        bedtools_runs.append(r)
                    else:
                        shown = f"subB={sub_b:g}" + (" +IE" if use_metrics else "")
                        print(f"  hammock-cpp Mode B {shown}...", end=" ", flush=True)
                        # Part 9 (2026-08-11): use_metrics=True used to mean
                        # "pass --metrics" -- the only way to get
                        # jaccard_similarity_ie pre-v0.8.0. It now means "use
                        # the bare/ie_only default", which gives the same
                        # column for less write cost (same union-pass compute
                        # cost, see run_hammock's docstring). The tuple/print
                        # label keeps saying "use_metrics"/"+IE" since it's
                        # still selecting the IE-column arm, just via a
                        # cheaper flag than before.
                        r = run_hammock(binary, file1_list_path, file2_list_path,
                                        precision, num_threads, sub_b=sub_b,
                                        ie_only=use_metrics)
                        r["sort_time"] = sort_time
                        r["run_index"] = run_i
                        runs_by_tool[tool].append(r)
                    rss = r["max_rss_mb"]
                    rss_str = f"{rss:.1f} MB" if rss is not None else "n/a"
                    print(f"{r['wall_time']:.2f}s wall, {rss_str} max RSS")

        entry: Dict[str, Any] = {
            "num_files": num_files,
            "num_intervals_per_file": num_intervals,
            "num_threads": num_threads,
            "precision": precision,
            "sub_b_list": sub_b_list,
            # Anchor for every downstream consumer. sub_b alone no longer
            # identifies an arm: the IE arm shares subB=1.0 with the baseline.
            "hammock_arms": [(label, sub_b) for label, sub_b, _ in arms],
            # None, not aggregate([]), when --no-bedtools: an empty aggregate
            # would be a dict of Nones that reads as "measured, and missing"
            # rather than "not run". Consumers branch on the key being None.
            "bedtools": aggregate(bedtools_runs, bedtools_keys) if with_bedtools else None,
            # Raw per-replicate rows, kept alongside the aggregate rather than
            # instead of it. aggregate() above reduces to mean/std/min/max and
            # is what every existing consumer (write_csv, plot, the R scripts)
            # reads -- that stays unchanged. But a wrong SIGN can hide behind a
            # correct-looking std at n=3: Panel A's N=2 cell aggregated to
            # "hammock 2.02x faster" from three reps of 0.14/0.50/1.22s, and a
            # 20-rep re-check of the same cell found the true, tight answer is
            # the opposite (bedtools 1.53x faster, 20/20). That could not be
            # recovered from the aggregate CSV after the fact -- only a rerun
            # could show it, because the raw values were never written down.
            # Leading underscore: private to this module, not part of the
            # public `results` schema anything else should read.
            "_raw_bedtools_runs": list(bedtools_runs),
            "_raw_hammock_runs": {t: list(rs) for t, rs in runs_by_tool.items()},
        }
        for tool, runs_for in runs_by_tool.items():
            entry[tool] = aggregate(runs_for, hammock_keys)
        results.append(entry)

        # Checkpoint after every num_files block, not just at the very end.
        # A long sweep (hours, many N values) previously kept everything in
        # memory and wrote nothing until the whole list finished -- an
        # interruption anywhere (timeout, scancel, node failure) lost the
        # entire run, including N values that had already completed cleanly.
        # See docs/seed-subsampling-synthetic-supplement.md: job 29758041 lost
        # all of N=2..1024 this way when it was cancelled partway into N=2048.
        if checkpoint_cb is not None:
            checkpoint_cb(results)

    return results


# os.umask() has no read-only form -- setting and immediately restoring it is
# the standard idiom to read the process umask without a side effect.
_UMASK = os.umask(0)
os.umask(_UMASK)


@contextlib.contextmanager
def _atomic_open_write(path: str, **open_kwargs):
    """Like open(path, "w", ...) but the target file is only ever replaced by
    a complete write. Writes to a sibling temp file in the same directory
    (same filesystem, so os.replace is atomic) and renames it onto `path`
    only after the `with` block exits cleanly.

    Exists because write_text_report/write_csv/write_runs_csv are now called
    once per num_files block (see run_benchmark's checkpoint_cb) specifically
    so an interrupted long sweep still leaves usable output on disk. A plain
    open(path, "w") truncates the file up front, so a kill (SLURM timeout,
    scancel, node failure -- the exact events checkpointing exists to survive)
    landing mid-write destroys the last good checkpoint instead of just
    failing to add a new one. Confirmed missing by review after job 29758041's
    incident writeup (docs/seed-subsampling-synthetic-supplement.md).
    """
    fd, tmp_path = tempfile.mkstemp(
        dir=os.path.dirname(os.path.abspath(path)) or ".",
        prefix=os.path.basename(path) + ".tmp")
    try:
        # mkstemp always creates the file mode 0600, unlike open(path, "w")
        # which follows the umask (0644 under the standard 022) -- match that
        # so checkpointed results stay group-readable on shared storage like
        # every other file this script writes.
        os.chmod(tmp_path, 0o666 & ~_UMASK)
        with os.fdopen(fd, "w", **open_kwargs) as f:
            yield f
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def write_text_report(results: List[Dict[str, Any]], path: str) -> None:
    with _atomic_open_write(path) as f:
        f.write("hammock-cpp (Mode B) vs bedtools\n")
        f.write("=" * 80 + "\n\n")
        if results:
            f.write(f"intervals/file: {results[0]['num_intervals_per_file']}\n")
            f.write(f"threads:        {results[0]['num_threads']}\n")
            f.write(f"HLL precision:  {results[0]['precision']}\n")
            f.write(f"subB values:    {results[0].get('sub_b_list', [1.0])}\n")
            f.write(f"system:         {get_system_info()}\n\n")
        for r in results:
            bt = r["bedtools"]
            f.write(f"num_files={r['num_files']}\n")
            f.write("-" * 80 + "\n")
            if bt is None:
                f.write("bedtools: not run (--no-bedtools)\n")
                for tool, sub_b in arms_of(r):
                    hm = r[tool]
                    ie = ", +IE columns" if (tool == IE_ARM_TOOL or tool.startswith("hammock_ie_B_subB")) else ""
                    f.write(f"hammock-cpp Mode B [subB={sub_b:g}, mixed-stride{ie}]:\n")
                    f.write(f"  wall:    {hm['mean_wall_time']:.3f} +/- {hm['std_wall_time']:.3f} s\n")
                    f.write(f"  cpu:     {hm['mean_cpu_time']:.3f} s\n")
                    if hm["mean_max_rss_mb"] is not None:
                        f.write(f"  max RSS: {hm['mean_max_rss_mb']:.1f} MB\n")
                    if hm["mean_sketch_creation_time"] is not None:
                        f.write(f"  sketch:  {hm['mean_sketch_creation_time']:.3f} s\n")
                    if hm["mean_comparison_time"] is not None:
                        f.write(f"  pairs:   {hm['mean_comparison_time']:.3f} s\n")
                f.write("\n")
                continue
            f.write("bedtools:\n")
            f.write(f"  wall:    {bt['mean_wall_time']:.3f} ± {bt['std_wall_time']:.3f} s "
                    f"(min {bt['min_wall_time']:.3f}, max {bt['max_wall_time']:.3f})\n")
            f.write(f"  cpu:     {bt['mean_cpu_time']:.3f} ± {bt['std_cpu_time']:.3f} s\n")
            if bt["mean_max_rss_mb"] is not None:
                f.write(f"  max RSS: {bt['mean_max_rss_mb']:.1f} MB "
                        f"(min {bt['min_max_rss_mb']:.1f}, max {bt['max_max_rss_mb']:.1f})\n")
            if bt.get("mean_sort_time") is not None:
                f.write(f"  sort:    {bt['mean_sort_time']:.3f} ± {bt['std_sort_time']:.3f} s "
                        f"(pre-sort, not in wall above; bedtools-workflow wall = wall + sort)\n")
            for tool, sub_b in arms_of(r):
                hm = r[tool]
                ie = ", +IE columns" if (tool == IE_ARM_TOOL or tool.startswith("hammock_ie_B_subB")) else ""
                f.write(f"hammock-cpp Mode B [subB={sub_b:g}, mixed-stride{ie}]:\n")
                f.write(f"  wall:    {hm['mean_wall_time']:.3f} ± {hm['std_wall_time']:.3f} s\n")
                f.write(f"  cpu:     {hm['mean_cpu_time']:.3f} ± {hm['std_cpu_time']:.3f} s\n")
                if hm["mean_max_rss_mb"] is not None:
                    f.write(f"  max RSS: {hm['mean_max_rss_mb']:.1f} MB\n")
                if hm["mean_sketch_creation_time"] is not None:
                    f.write(f"  sketch:  {hm['mean_sketch_creation_time']:.3f} s\n")
                if hm["mean_comparison_time"] is not None:
                    f.write(f"  pairs:   {hm['mean_comparison_time']:.3f} s\n")
                if hm["mean_wall_time"] and bt["mean_wall_time"]:
                    f.write(f"  speedup (wall): {bt['mean_wall_time'] / hm['mean_wall_time']:.2f}x\n")
                if hm["mean_cpu_time"] and bt["mean_cpu_time"]:
                    f.write(f"  speedup (cpu):  {bt['mean_cpu_time'] / hm['mean_cpu_time']:.2f}x\n")
            f.write("\n")


# The same fields the pairwise-cost benchmark carries, for the same reason: an
# archived CSV that records neither the hardware nor the allocation cannot be
# checked for comparability afterwards, and -march=native makes the CPU model
# load-bearing. get_system_info() is also written to the sibling .txt report,
# but a report is not joinable and nothing reads it back.
# binary_version is NOT redundant with git_sha: HAMMOCK_CPP_BIN and --binary both
# bypass the build tree, so the binary that ran can be a stale site-packages copy
# while the repo sits at a newer SHA -- the exact case check_binary_version's
# docstring warns about.
PROVENANCE_COLS = ["hostname", "cpu_model", "cpu_count", "slurm_job_id", "git_sha",
                   "binary_version", "corpus_seed"]


def write_csv(results: List[Dict[str, Any]], path: str,
              provenance: Optional[Dict[str, Any]] = None) -> None:
    if provenance is None:
        provenance = get_system_info()
    prov_row = [provenance.get(c, "unknown") for c in PROVENANCE_COLS]
    cols = [
        "num_files", "num_threads", "precision", "sub_b", "tool",
        "mean_wall_time", "std_wall_time", "min_wall_time", "max_wall_time",
        "mean_cpu_time", "std_cpu_time", "min_cpu_time", "max_cpu_time",
        "mean_max_rss_mb", "std_max_rss_mb", "min_max_rss_mb", "max_max_rss_mb",
        "mean_sort_time", "std_sort_time", "min_sort_time", "max_sort_time",
        "mean_sketch_creation_time", "std_sketch_creation_time",
        "mean_comparison_time", "std_comparison_time",
        "mean_pair_time", "std_pair_time",
        "mean_write_time", "std_write_time",
        # bedtools-only; blank on hammock rows. Without these a reader has to
        # ASSUME the baseline used the cores it was given, and on these nodes it
        # does not: `bedtools jaccard` has no batch mode, so the workflow
        # launches one process per pair, and process creation caps near
        # 123 exec/s regardless of core count. mean_bedtools_parallel_eff is the
        # fraction of `num_threads` the baseline actually achieved -- measured
        # ~1.7x out of 16x. A speedup quoted without it is not interpretable.
        # Read it only at large N. At N=2 or 4 the whole bedtools leg is a few
        # pairs behind a fixed ~0.5 s of module load and process setup, so the
        # ratio reads ~0.01-0.04 and means "startup-dominated", not "no
        # parallelism". It becomes meaningful once n_pairs * serial >> startup,
        # i.e. from about N=32 up.
        "mean_bedtools_serial_ms", "std_bedtools_serial_ms",
        "mean_bedtools_parallel_eff", "std_bedtools_parallel_eff",
        "min_bedtools_parallel_eff", "max_bedtools_parallel_eff",
    ]
    metric_cols = cols[5:]          # resolved from the aggregate dict
    cols = cols + PROVENANCE_COLS   # resolved from provenance, not from `d`
    with _atomic_open_write(path, newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in results:
            tools = [("bedtools", "bedtools", "")] if r.get("bedtools") is not None else []
            for tool, sub_b in arms_of(r):
                tools.append((tool, tool, f"{sub_b:g}"))
            for tool, key, sub_b_str in tools:
                d = r[key]
                row = [r["num_files"], r["num_threads"], r["precision"], sub_b_str, tool]
                for c in metric_cols:
                    v = d.get(c)
                    row.append(f"{v:.6f}" if isinstance(v, float) else ("" if v is None else v))
                row.extend(prov_row)
                w.writerow(row)


RUNS_METRIC_COLS = [
    "wall_time", "cpu_time", "max_rss_mb", "sort_time",
    "sketch_creation_time", "comparison_time", "pair_time", "write_time",
    "bedtools_serial_ms", "bedtools_parallel_eff",
]


def write_runs_csv(results: List[Dict[str, Any]], path: str,
                   provenance: Optional[Dict[str, Any]] = None) -> None:
    """One row per (num_files, tool, replicate) -- the data aggregate() throws
    away. Not a replacement for the aggregate CSV: every existing consumer
    (write_csv, plot, the R scripts) keeps reading that one unchanged. This is
    what makes a confidence interval, a bootstrap, or a "was that cell's mean
    actually reliable" check possible without rerunning the job -- previously
    none of those were possible on an already-collected CSV at all.
    """
    if provenance is None:
        provenance = get_system_info()
    prov_row = [provenance.get(c, "unknown") for c in PROVENANCE_COLS]
    cols = (["num_files", "num_threads", "precision", "sub_b", "tool", "run_index"]
            + RUNS_METRIC_COLS + PROVENANCE_COLS)
    with _atomic_open_write(path, newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in results:
            for run in r.get("_raw_bedtools_runs", []):
                row = [r["num_files"], r["num_threads"], r["precision"], "",
                      "bedtools", run.get("run_index", "")]
                for c in RUNS_METRIC_COLS:
                    v = run.get(c)
                    row.append(f"{v:.6f}" if isinstance(v, float) else ("" if v is None else v))
                row.extend(prov_row)
                w.writerow(row)
            for tool, runs_for in r.get("_raw_hammock_runs", {}).items():
                sub_b = next((s for t, s in arms_of(r) if t == tool), None)
                sub_b_str = f"{sub_b:g}" if sub_b is not None else ""
                for run in runs_for:
                    row = [r["num_files"], r["num_threads"], r["precision"],
                          sub_b_str, tool, run.get("run_index", "")]
                    for c in RUNS_METRIC_COLS:
                        v = run.get(c)
                        row.append(f"{v:.6f}" if isinstance(v, float) else ("" if v is None else v))
                    row.extend(prov_row)
                    w.writerow(row)


def plot(results: List[Dict[str, Any]], png_path: str) -> None:
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        print("matplotlib not available, skipping plot")
        return

    if any(r.get("bedtools") is None for r in results):
        # Every panel of this plot is drawn against the bedtools series.
        # A hammock-only run (--no-bedtools) has no baseline to plot, so
        # skip rather than emit axes with one curve missing and no note.
        print("no bedtools arm in these results, skipping the comparison plot")
        return
    nfiles = [r["num_files"] for r in results]
    bt_wall = [r["bedtools"]["mean_wall_time"] for r in results]
    bt_cpu = [r["bedtools"]["mean_cpu_time"] for r in results]
    bt_rss = [r["bedtools"]["mean_max_rss_mb"] for r in results]
    threads = results[0]["num_threads"] if results else 1
    arms = arms_of(results[0]) if results else []

    hm_colors = ["#ff7f0e", "#1f77b4", "#2ca02c", "#d62728"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].plot(nfiles, bt_wall, "ko-", label=f"bedtools (t={threads})", linewidth=2)
    for i, (tool, sub_b) in enumerate(arms):
        ys = [r[tool]["mean_wall_time"] for r in results]
        axes[0].plot(nfiles, ys, "s--", color=hm_colors[i % len(hm_colors)],
                     label=arm_legend(tool, sub_b), linewidth=2)
    axes[0].set_xlabel("Number of files")
    axes[0].set_ylabel("Wall time (s)")
    axes[0].set_title("Wall time")
    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(nfiles, bt_cpu, "ko-", label="bedtools", linewidth=2)
    for i, (tool, sub_b) in enumerate(arms):
        ys = [r[tool]["mean_cpu_time"] for r in results]
        axes[1].plot(nfiles, ys, "s--", color=hm_colors[i % len(hm_colors)],
                     label=arm_legend(tool, sub_b), linewidth=2)
    axes[1].set_xlabel("Number of files")
    axes[1].set_ylabel("CPU time (s)")
    axes[1].set_title("CPU time")
    axes[1].set_xscale("log", base=2)
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    hm_rss_any = any(r[tool]["mean_max_rss_mb"] is not None
                     for r in results for tool, _ in arms)
    if any(v is not None for v in bt_rss) or hm_rss_any:
        axes[2].plot(nfiles, bt_rss, "ko-", label="bedtools", linewidth=2)
        for i, (tool, sub_b) in enumerate(arms):
            ys = [r[tool]["mean_max_rss_mb"] for r in results]
            axes[2].plot(nfiles, ys, "s--", color=hm_colors[i % len(hm_colors)],
                         label=arm_legend(tool, sub_b), linewidth=2)
        axes[2].set_xlabel("Number of files")
        axes[2].set_ylabel("Max RSS (MB)")
        axes[2].set_title("Peak memory")
        axes[2].set_xscale("log", base=2)
        axes[2].grid(True, alpha=0.3)
        axes[2].legend()
    else:
        axes[2].set_visible(False)

    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark hammock-cpp Mode B vs bedtools")
    parser.add_argument("--threads", "-t", type=int, default=1)
    parser.add_argument("--precision", "-p", type=int, default=16)
    parser.add_argument("--num-intervals", type=int, default=NUM_INTERVALS_PER_FILE)
    parser.add_argument("--num-files", type=str, default=",".join(map(str, NUM_FILES_LIST)),
                        help="Comma-separated file counts")
    parser.add_argument("--runs", type=int, default=NUM_RUNS)
    parser.add_argument("--subB-list", dest="sub_b_list", type=str, default="1.0",
                        help="Comma-separated subB values (default '1.0'). Each value < 1.0 "
                             "adds a hammock variant invoked with --subB <val> --subB-method "
                             "mixed-stride. subB=1.0 emits 'hammock_cpp_B' for backwards compat; "
                             "other values emit 'hammock_cpp_B_subB<val>'.")
    parser.add_argument("--metrics-arm", action="store_true",
                        help="Add a fourth hammock arm at subB=1.0 run with the bare/no-flag "
                             "default (1 column, jaccard_similarity_ie) -- instead of the "
                             "baseline arm's cheap register-equality flag. Still pays the fused "
                             "union pass (same compute cost as --metrics), just not the extra "
                             "write cost of 7 unused columns. "
                             f"Labelled '{IE_ARM_TOOL}' so it fails the "
                             "'^hammock_cpp_B' filter the R consumers use, and cannot double "
                             "rows in their joins. Ignored if --metrics-all is also given. "
                             "Retargeted 2026-08-11 (Part 9, docs/seed-metrics-column-restructure.md): "
                             "used to pass --metrics (the full 8-column block) because that was "
                             "the only way to get jaccard_similarity_ie before v0.8.0; now uses "
                             "ie_only in run_hammock (bare default) instead. See "
                             "experiments/bedtools_benchmark/RESULTS.md for the paired "
                             "measurement backing the switch.")
    parser.add_argument("--metrics-all", action="store_true",
                        help="Run EVERY arm in --subB-list with the bare/no-flag default (1 "
                             "column, jaccard_similarity_ie) instead of the cheap "
                             "register-equality flag (replaces the default arms, rather "
                             "than adding one extra arm the way --metrics-arm does). Use this when "
                             "the figure needs jaccard_similarity_ie wall times at subB < 1.0, "
                             "e.g. a subB-vs-N line plot -- --metrics-arm alone can only give "
                             "that column at subB=1.0. Labels via ie_tool_name_for_subb(): "
                             "'hammock_ie_B' at subB=1.0 (byte-compatible with --metrics-arm's "
                             "label), 'hammock_ie_B_subB<val>' otherwise -- both fail every "
                             "existing R consumer's filter, so this cannot corrupt an existing "
                             "figure's data even if pointed at the same results/ directory.")
    parser.add_argument("--no-bedtools", dest="with_bedtools", action="store_false",
                        default=True,
                        help="Skip the bedtools arm entirely. For extending the N axis past "
                             "the point where an exact baseline is affordable: bedtools is "
                             "Theta(N^2) with a large constant (~2.2 h per replicate at "
                             "N=1024, ~8.6 h at N=2048), while hammock is near-linear in "
                             "this regime because at p=18 it is sketching-dominated. The "
                             "resulting CSV has NO bedtools rows, so plot_pairwise_scaling.R "
                             "cannot consume it -- it is for measuring hammock at catalog "
                             "scale and projecting bedtools from its own fitted curve, "
                             "which must be labelled as a projection wherever it is shown.")
    parser.add_argument("--test", action="store_true", help="Quick test (small files, few runs)")
    parser.add_argument("--corpus-seed", type=int, default=None,
                        help="Seed the synthetic BED corpus so a re-run draws the "
                             "same files. Default (unseeded) matches every archived "
                             "run, but makes cross-run comparison good only to a few "
                             "percent on the bedtools leg, whose cost -- unlike "
                             "hammock's -- depends on the interval content. "
                             "Replicates still draw different corpora; it is the "
                             "sequence that becomes reproducible.")
    parser.add_argument("--binary", type=str, default=None, help="Path to hammock-cpp")
    parser.add_argument("--output-dir", type=str, default=RESULTS_DIR,
                        help="Directory for txt and csv reports")
    parser.add_argument("--figures-dir", type=str, default=FIGURES_DIR,
                        help="Directory for PNG plots")
    args = parser.parse_args()

    binary = args.binary or find_hammock_cpp()
    if not os.path.exists(binary):
        print(f"hammock-cpp not found at {binary}", file=sys.stderr)
        return 1
    # Probe the resolved path, not just the glob: --binary and HAMMOCK_CPP_BIN
    # both bypass find_hammock_cpp()'s search, and the env var is what the
    # sbatch scripts set.
    try:
        binary_version = check_binary_version(binary)
        print(f"hammock-cpp version: {binary_version}")
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 1
    if not os.path.exists(BEDTOOLS_SCRIPT):
        print(f"bedtools.sh not found at {BEDTOOLS_SCRIPT}", file=sys.stderr)
        return 1

    if args.test:
        num_files_list = [2, 4]
        num_runs = 2
        num_intervals = 1000
    else:
        num_files_list = [int(x) for x in args.num_files.split(",")]
        num_runs = args.runs
        num_intervals = args.num_intervals

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.figures_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = f"cpp_vs_bedtools_t{args.threads}_{timestamp}"
    txt_path = os.path.join(args.output_dir, stem + ".txt")
    csv_path = os.path.join(args.output_dir, stem + ".csv")
    png_path = os.path.join(args.figures_dir, stem + ".png")

    sub_b_list = [float(x) for x in args.sub_b_list.split(",") if x.strip()]
    if not sub_b_list:
        print("--subB-list must contain at least one value", file=sys.stderr)
        return 1
    for s in sub_b_list:
        if not (0.0 < s <= 1.0):
            print(f"subB values must lie in (0, 1.0]; got {s}", file=sys.stderr)
            return 1

    # Built before run_benchmark so the checkpoint callback below can write a
    # complete, correctly-provenanced CSV/report after every num_files block,
    # not just once at the end (see run_benchmark's checkpoint_cb docstring
    # note).
    prov = get_system_info()
    prov["binary_version"] = binary_version
    # Empty, not the string "unseeded": the column is otherwise an integer, and a
    # mixed character/numeric column makes a seeded and an unseeded CSV fail to
    # bind_rows in R. Empty reads as NA in both R and pandas.
    prov["corpus_seed"] = args.corpus_seed if args.corpus_seed is not None else ""
    runs_csv_path = re.sub(r"\.csv$", "_runs.csv", csv_path)

    def _checkpoint(partial_results: List[Dict[str, Any]]) -> None:
        write_text_report(partial_results, txt_path)
        write_csv(partial_results, csv_path, prov)
        write_runs_csv(partial_results, runs_csv_path, prov)

    results = run_benchmark(
        binary=binary,
        num_files_list=num_files_list,
        num_runs=num_runs,
        num_intervals=num_intervals,
        num_threads=args.threads,
        precision=args.precision,
        sub_b_list=sub_b_list,
        metrics_arm=args.metrics_arm,
        metrics_all=args.metrics_all,
        corpus_seed=args.corpus_seed,
        with_bedtools=args.with_bedtools,
        checkpoint_cb=_checkpoint,
    )

    # Final write is not redundant with the last checkpoint: it's identical in
    # content (same results list) but guarantees a write happened even if
    # num_files_list were ever empty, and keeps this call site independent of
    # checkpointing being wired up correctly above.
    write_text_report(results, txt_path)
    write_csv(results, csv_path, prov)
    write_runs_csv(results, runs_csv_path, prov)
    plot(results, png_path)
    print(f"\nReports: {txt_path}\n         {csv_path}\n         {runs_csv_path} (raw per-replicate rows)\nFigure:  {png_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
