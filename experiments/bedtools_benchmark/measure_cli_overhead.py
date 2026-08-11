#!/usr/bin/env python3
"""Quantify the Python-CLI-vs-hammock-cpp overhead gap.

Every benchmark figure in this paper times `hammock-cpp` (the standalone
binary, no Python interpreter in the loop), not the `hammock` Python CLI end
users actually run -- see the language drafted for the methods/limitations
text. This script measures the gap directly instead of leaving it as a named
but unquantified caveat: same corpus, same flags (mode B, precision,
threads), same machine/allocation, alternating which tool goes first each
replicate (the same _rotate-style anti-confound benchmark_cpp_vs_bedtools.py
uses elsewhere).

The `hammock` CLI has no --metrics/--no-metrics opt-out (CLAUDE.md: it always
emits the full 9-column block, matching hammock-cpp's default since v0.7.0).
hammock-cpp's arm is controlled by --cpp-metrics-arm: 'no-metrics' (the
default, since 2026-08-11) matches what benchmark_cpp_vs_bedtools.py actually
passes for every headline figure; 'metrics' reproduces the original
apples-to-apples 9-column comparison (job 29758101,
results/cli_overhead_1786473951.csv) for continuity with that data. See
docs/seed-hammock-cpp-file-dispatch.md Part 2's "Important scope note": the
crossover was originally measured only on the --metrics arm, which is NOT the
arm the headline figures use, so a --no-metrics remeasurement was flagged as
the natural next step. Note the CLI's own wall time is unaffected either way
(it always emits 9 columns) -- only the cpp side's arm changes here.

Extended to N up to 2048, where a single run can take hours, this picked up
two hardenings mirrored from benchmark_cpp_vs_bedtools.py's checkpoint/pin
fix (commits 4f288d3/b8a082a), which exist because job 29758041 lost an
entire multi-hour sweep -- including N values that had already finished
cleanly -- to a concurrent binary rebuild plus a single end-of-run CSV write:

- Checkpointed, atomically-written CSV after every num_files block (see
  _atomic_open_write, reused from that module) instead of one write at the
  very end. A killed/timed-out run now costs at most the in-progress N.
- A provenance fingerprint (git HEAD, dirty status of the source paths that
  feed either binary, and both binaries' mtime/size) captured once at start
  and re-checked at every checkpoint. IMPORTANT ASYMMETRY, discovered while
  wiring this up: hammock-cpp is a single self-contained binary and CAN be
  fully pinned to a job-local copy (see sbatch_cli_overhead.sh) the same way
  the sibling script pins it. The `hammock` CLI CANNOT be pinned the same
  way -- `claude-ref-comparison` is a pip editable install, and empirically
  (checked by hand, not assumed) prepending a job-local copy to PYTHONPATH
  does NOT shadow it: `hammock.__file__` still resolves to the live repo's
  `python/hammock/` tree even with a copy earlier on PYTHONPATH, because the
  editable install's meta-path finder intercepts before normal sys.path
  resolution. So the CLI's pure-Python behavior is live-linked to whatever
  is on disk in this repo for the whole duration of the run, by construction
  -- an edit to python/hammock/*.py elsewhere needs no rebuild at all to
  change what's being measured mid-sweep. The fingerprint below is a
  *detection* mechanism for that risk, not a fix for it: it cannot stop
  drift, it can only flag rows measured after drift was observed.

Small and quick by design (order-of-magnitude characterization for a
limitations sentence, not a headline figure): a handful of N points, one
subB, one precision, alternated tool order per replicate.
"""
import argparse
import concurrent.futures
import glob
import os
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
sys.path.insert(0, SCRIPT_DIR)
from benchmark_cpp_vs_bedtools import (  # noqa: E402
    generate_bed_file, derive_seed, get_system_info, check_binary_version,
    _sort_one, _atomic_open_write,
)


def run_with_wall(cmd: List[str]) -> Dict[str, Any]:
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed (exit {r.returncode}): {r.stderr[-2000:]}")
    return {"wall_time": wall}


def _gen_one(job: Tuple[int, str, Optional[int]]) -> None:
    """Top-level (picklable) wrapper for ProcessPoolExecutor.map."""
    num_intervals, path, seed = job
    generate_bed_file(num_intervals, path, seed)


def _cli_core_so_glob(hammock_cli_bin: str) -> Optional[str]:
    """Locate the CLI's compiled `hammock._core` extension from its bin path.

    `<env>/bin/hammock` -> `<env>/lib/python3.*/site-packages/hammock/_core*.so`.
    Best-effort: returns None (not an error) if the layout doesn't match,
    since this only feeds the provenance fingerprint, not correctness.

    Canonicalized via realpath: this conda env has a `python3.1 -> python3.10`
    symlink alongside the real `python3.10` directory, so the glob can match
    both and `glob.glob`'s ordering isn't guaranteed -- without this, two runs
    against an unchanged install could print two different-looking paths for
    the same file (harmless for the fingerprint's mtime/size comparison,
    which follows symlinks either way, but confusing to read and worth
    pinning down since it was found while smoke-testing, not assumed away).
    """
    env_root = os.path.dirname(os.path.dirname(os.path.abspath(hammock_cli_bin)))
    hits = glob.glob(os.path.join(env_root, "lib", "python3.*", "site-packages",
                                   "hammock", "_core*.so"))
    if not hits:
        return None
    return sorted({os.path.realpath(h) for h in hits})[0]


def _stat_or_none(path: Optional[str]) -> Optional[Tuple[float, int]]:
    if not path:
        return None
    try:
        st = os.stat(path)
        return (st.st_mtime, st.st_size)
    except OSError:
        return None


def provenance_fingerprint(hammock_cpp_bin: str, hammock_cli_bin: str) -> Dict[str, Any]:
    """What could silently change either tool's behavior mid-run, cheap to
    check. Not a pin -- see module docstring -- a detector: compared against
    the start-of-run fingerprint once per replicate, a mismatch means some
    rows may not be comparable to others in the same CSV.

    Point-in-time, not a log: `git status --porcelain` only sees the state at
    the moment it runs. An edit to python/hammock/*.py (or cpp/, bindings/)
    that lands and is reverted between two checks -- including one that
    genuinely changed CLI behavior for whatever ran during that window -- is
    invisible to this detector. A provenance_drift=False run is "no
    *persistent* change was observed at replicate granularity," not "verified
    uncontaminated."
    """
    fp: Dict[str, Any] = {}
    try:
        fp["git_head"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True,
            text=True, check=True).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--", "python", "cpp", "bindings"],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout
        fp["dirty_relevant_paths"] = sorted(
            line[3:] for line in dirty.splitlines() if line.strip())
    except (subprocess.CalledProcessError, OSError) as e:
        fp["git_error"] = str(e)
    fp["hammock_cpp_bin_stat"] = _stat_or_none(hammock_cpp_bin)
    fp["cli_core_so_path"] = _cli_core_so_glob(hammock_cli_bin)
    fp["cli_core_so_stat"] = _stat_or_none(fp["cli_core_so_path"])
    return fp


def fingerprint_diff(a: Dict[str, Any], b: Dict[str, Any]) -> List[str]:
    """Human-readable list of what changed between two fingerprints, empty if none.

    A transient `git` failure (e.g. a concurrent commit/gc briefly locking
    .git on shared GPFS during a multi-hour job) leaves git_head/
    dirty_relevant_paths absent from that side's fingerprint rather than
    populated -- comparing via plain .get() would then read as "None ->
    real SHA", a spurious drift alarm for a run where nothing actually
    changed. Skip a key on either side that failed to capture (git_error
    set) instead of treating "unknown" as "different"; a real git_error is
    reported separately, not folded into this diff.
    """
    changes = []
    for key in ("git_head", "dirty_relevant_paths", "hammock_cpp_bin_stat",
                "cli_core_so_stat"):
        if "git_error" in a or "git_error" in b:
            if key in ("git_head", "dirty_relevant_paths"):
                continue  # can't trust either side's git-derived value here
        if a.get(key) != b.get(key):
            changes.append(f"{key}: {a.get(key)!r} -> {b.get(key)!r}")
    if "git_error" in b and "git_error" not in a:
        changes.append(f"git_error (new): {b['git_error']!r} -- git-derived "
                        f"fields skipped above, not compared")
    return changes


def write_checkpoint(rows: List[Dict[str, Any]], out_path: str, fieldnames: List[str]) -> None:
    import csv
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with _atomic_open_write(out_path, newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def print_n_summary(num_files: int, by_cell: Dict[Tuple[int, str], List[float]],
                     header: bool) -> None:
    import statistics
    if header:
        print(f"\n{'N':>6} {'cpp_med':>10} {'cpp_[min,max]':>16} "
              f"{'cli_med':>10} {'cli_[min,max]':>16} {'cli/cpp':>9}")
    cpp = by_cell[(num_files, "hammock_cpp")]
    cli = by_cell[(num_files, "hammock_cli")]
    if not cpp or not cli:
        return
    cpp_m, cli_m = statistics.median(cpp), statistics.median(cli)
    print(f"{num_files:>6} {cpp_m:>10.4f} "
          f"[{min(cpp):>6.3f},{max(cpp):>6.3f}] "
          f"{cli_m:>10.4f} [{min(cli):>6.3f},{max(cli):>6.3f}] "
          f"{cli_m / cpp_m:>9.3f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--num-files", default="2,8,32,128")
    ap.add_argument("--num-intervals", type=int, default=10000)
    ap.add_argument("--precision", type=int, default=18)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--corpus-seed", type=int, default=20260811)
    ap.add_argument("--cpp-metrics-arm", choices=["no-metrics", "metrics"],
                     default="no-metrics",
                     help="hammock-cpp arm to time (default: no-metrics, "
                          "matching what benchmark_cpp_vs_bedtools.py passes "
                          "for every headline figure; 'metrics' reproduces "
                          "the original job-29758101 comparison). The CLI "
                          "always emits the 9-column block either way.")
    ap.add_argument("--hammock-cpp-bin", default=os.environ.get(
        "HAMMOCK_CPP_BIN",
        "/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp"))
    ap.add_argument("--hammock-cli-bin", default=
        "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock")
    ap.add_argument("--output", default=None)
    ap.add_argument("--gen-workers", type=int, default=None,
                     help="Process-pool size for corpus generation (default: "
                          "--threads). Generation is pure-Python/GIL-bound, "
                          "so this uses processes, not threads, to actually "
                          "parallelize -- previously a fully serial "
                          "generation+sort loop (sort was the bigger of the "
                          "two, ~5-7x generation's per-file cost) ate up to "
                          "~50%% of a cell's wall time at N=32 before both "
                          "were parallelized.")
    args = ap.parse_args()

    for b in (args.hammock_cpp_bin, args.hammock_cli_bin):
        if not os.path.exists(b):
            print(f"not found: {b}", file=sys.stderr)
            return 1
    check_binary_version(args.hammock_cpp_bin)

    print(f"hammock-cpp: {args.hammock_cpp_bin}")
    print(f"hammock CLI: {args.hammock_cli_bin}")
    sysinfo = get_system_info()
    print(f"system: {sysinfo}")
    if args.runs % 2 != 0:
        print(f"WARNING: --runs {args.runs} is not a multiple of the 2-tool "
              f"alternation cycle -- the first-mover split will be uneven "
              f"(docs/seed-benchmark-methodology.md item 7). Prefer an even "
              f"--runs.", file=sys.stderr)

    start_fp = provenance_fingerprint(args.hammock_cpp_bin, args.hammock_cli_bin)
    print(f"provenance at start: {start_fp}")
    drifted = False

    num_files_list = [int(x) for x in args.num_files.split(",")]
    print(f"planned N sweep for this arm ({args.cpp_metrics_arm}): "
          f"{num_files_list} x {args.runs} runs each. A CSV snapshot mid-run "
          f"does not self-describe how far along it is -- check this log's "
          f"most recent '[checkpoint after N=...]' / 'still to come' line, "
          f"or squeue/sacct for the job's state, before treating any given "
          f"snapshot as complete.")
    rows: List[Dict[str, Any]] = []
    fieldnames = ["num_files", "run", "tool", "wall_time", "cpp_metrics_arm",
                  "hostname", "cpu_model", "slurm_job_id", "git_sha",
                  "provenance_drift"]

    arm_tag = "_nometrics" if args.cpp_metrics_arm == "no-metrics" else "_metrics"
    out_path = args.output or os.path.join(
        SCRIPT_DIR, "results", f"cli_overhead{arm_tag}_{int(time.time())}.csv")

    import statistics  # noqa: F401 (used in print_n_summary)
    from collections import defaultdict
    by_cell: Dict[Tuple[int, str], List[float]] = defaultdict(list)

    gen_workers = args.gen_workers or max(1, args.threads)
    with concurrent.futures.ProcessPoolExecutor(max_workers=gen_workers) as gen_pool:
        for n_idx, num_files in enumerate(num_files_list):
            print(f"\n{'=' * 60}\n{num_files} files x {args.num_intervals} intervals\n{'=' * 60}")
            for run_i in range(args.runs):
                with tempfile.TemporaryDirectory() as tmp_dir:
                    file1_list, file2_list = [], []
                    gen_jobs = []
                    for i in range(num_files):
                        a = os.path.join(tmp_dir, f"set1_{i}.bed")
                        b = os.path.join(tmp_dir, f"set2_{i}.bed")
                        gen_jobs.append((args.num_intervals, a, derive_seed(args.corpus_seed, run_i, i, 0)))
                        gen_jobs.append((args.num_intervals, b, derive_seed(args.corpus_seed, run_i, i, 1)))
                        file1_list.append(a)
                        file2_list.append(b)
                    list(gen_pool.map(_gen_one, gen_jobs))

                    # sort (both tools require it identically; not timed) --
                    # parallel threads, same helper benchmark_cpp_vs_bedtools.py
                    # uses (subprocess.run releases the GIL while the child
                    # runs, so threads parallelize this fine; it's the
                    # CPU-bound generation above that needed processes).
                    all_paths = file1_list + file2_list
                    sort_workers = max(1, min(args.threads, len(all_paths)))
                    with concurrent.futures.ThreadPoolExecutor(max_workers=sort_workers) as sp:
                        for _ in sp.map(_sort_one, all_paths):
                            pass

                    l1 = os.path.join(tmp_dir, "l1.txt")
                    l2 = os.path.join(tmp_dir, "l2.txt")
                    with open(l1, "w") as f:
                        f.write("\n".join(file1_list) + "\n")
                    with open(l2, "w") as f:
                        f.write("\n".join(file2_list) + "\n")

                    cpp_out = os.path.join(tmp_dir, "cpp_out")
                    cli_out = os.path.join(tmp_dir, "cli_out")
                    cpp_metrics_flag = ("--metrics" if args.cpp_metrics_arm == "metrics"
                                         else "--no-metrics")
                    cpp_cmd = [args.hammock_cpp_bin, l1, l2, "--mode", "B",
                               "-p", str(args.precision), "-o", cpp_out,
                               "--threads", str(args.threads), cpp_metrics_flag]
                    cli_cmd = [args.hammock_cli_bin, l1, l2, "--mode", "B",
                               "--precision", str(args.precision), "-o", cli_out,
                               "--threads", str(args.threads)]

                    # Alternate order by run_i so neither tool systematically goes
                    # first (warm page cache / thermal state confound).
                    pair = [("hammock_cpp", cpp_cmd), ("hammock_cli", cli_cmd)]
                    if run_i % 2:
                        pair = pair[::-1]
                    for label, cmd in pair:
                        r = run_with_wall(cmd)
                        print(f"  {label}: {r['wall_time']:.4f}s wall")
                        rows.append({"num_files": num_files, "run": run_i,
                                     "tool": label, "wall_time": r["wall_time"],
                                     "cpp_metrics_arm": args.cpp_metrics_arm,
                                     "hostname": sysinfo.get("hostname"),
                                     "cpu_model": sysinfo.get("cpu_model"),
                                     "slurm_job_id": sysinfo.get("slurm_job_id"),
                                     "git_sha": sysinfo.get("git_sha"),
                                     "provenance_drift": drifted})
                        by_cell[(num_files, label)].append(r["wall_time"])

                # Fingerprint check once per run_i (both tools already ran),
                # not once per num_files block. The check itself is cheap
                # (git rev-parse + git status + two stats, sub-100ms measured)
                # against per-replicate runtimes of seconds to minutes, so
                # this is affordable at every granularity that matters: a
                # block-level check left the ENTIRE block in which drift
                # actually happened permanently mismarked clean (drift is
                # only visible at the *next* check, and the previous check
                # was up to `--runs` replicates in the past); per-run_i
                # narrows that blind spot to at most one run_i (2 tool
                # invocations) instead of up to `--runs` of them.
                cur_fp = provenance_fingerprint(args.hammock_cpp_bin, args.hammock_cli_bin)
                changes = fingerprint_diff(start_fp, cur_fp)
                if changes and not drifted:
                    drifted = True
                    print(f"WARNING: provenance drift detected after "
                          f"N={num_files} run {run_i} -- either binary's "
                          f"underlying source/artifact changed during or "
                          f"before this replicate: {changes}. Rows measured "
                          f"before this replicate are known pre-drift. The "
                          f"N={num_files} run={run_i} rows just appended are "
                          f"of UNKNOWN provenance (the change could have "
                          f"landed between the two tool invocations) despite "
                          f"being stamped provenance_drift=False -- detection "
                          f"granularity is per-replicate, not per-tool-call. "
                          f"Rows from the next replicate onward are stamped "
                          f"provenance_drift=True. This detector cannot see "
                          f"an edit that was reverted before this check ran "
                          f"(e.g. python/hammock/*.py touched and restored "
                          f"between two checks) -- it is a point-in-time "
                          f"comparison, not a log.", file=sys.stderr)

            # Checkpoint after every num_files block, not just at the very
            # end -- see module docstring. A long sweep previously kept
            # everything in memory and wrote nothing until the whole N list
            # finished; an interruption anywhere lost results for N values
            # that had already completed cleanly.
            write_checkpoint(rows, out_path, fieldnames)
            print(f"  [checkpoint after N={num_files}: wrote {out_path} "
                  f"({len(rows)} rows; still to come this arm: "
                  f"{num_files_list[n_idx + 1:] or 'none -- this is the last N'})]")
            print_n_summary(num_files, by_cell, header=(n_idx == 0))

    print(f"\nfinal write: {out_path}")
    if drifted:
        print(f"WARNING: this run experienced provenance drift -- see stderr "
              f"above and the provenance_drift column. Treat any cross-N "
              f"comparison spanning the drift point with caution.",
              file=sys.stderr)

    print(f"\n{'N':>6} {'cpp_med':>10} {'cpp_[min,max]':>16} "
          f"{'cli_med':>10} {'cli_[min,max]':>16} {'cli/cpp':>9}")
    for n in num_files_list:
        print_n_summary(n, by_cell, header=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
