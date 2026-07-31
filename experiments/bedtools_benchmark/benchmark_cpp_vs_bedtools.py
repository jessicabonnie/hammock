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
from typing import Any, Dict, List, Optional

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

SKETCH_RE = re.compile(r"^Sketching:\s+(\d+)\s+ms")
PAIR_RE = re.compile(r"^Pairwise\+write:\s+(\d+)\s+ms")
MAXRSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s+(\d+)")


def get_system_info() -> Dict[str, Any]:
    info = {
        "cpu_count": os.cpu_count(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
    }
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    info["memory_total_gb"] = int(line.split()[1]) / (1024 * 1024)
                    break
    except OSError:
        pass
    return info


def generate_bed_file(num_intervals: int, output_file: str) -> None:
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    with open(output_file, "w") as f:
        for _ in range(num_intervals):
            chrom = random.choice(chroms)
            start = random.randint(0, 10_000_000)
            end = random.randint(start + 100, start + 10_000)
            f.write(f"{chrom}\t{start}\t{end}\n")


def run_with_time(cmd: List[str]) -> Dict[str, Any]:
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
        result = subprocess.run(wrapped, capture_output=True, text=True, check=True)
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
    return run_with_time(["bash", BEDTOOLS_SCRIPT, file1_list, file2_list, str(num_threads)])


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
) -> Dict[str, Any]:
    """Run hammock-cpp Mode B and return timing/RSS.

    With keep_output=True, the output CSV is left in place and its path is
    returned as result["output_csv"]; the caller is responsible for deleting it.

    sub_b < 1.0 enables point subsampling using sub_b_method (default
    mixed-stride, matching the post-9778ef8 binary default). At sub_b == 1.0
    we omit the flags to keep the cmd line byte-identical to pre-subB runs.

    metrics=True adds --metrics, which emits jaccard_similarity_ie and the
    containment/cosketch block. It costs a union + cardinality per pair, so a
    run with it on is NOT timing-comparable to the published numbers in
    RESULTS.md -- use it only on untimed accuracy passes.
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
        if metrics:
            cmd += ["--metrics"]
        r = run_with_time(cmd)

        sketch_s: Optional[float] = None
        pair_s: Optional[float] = None
        for line in r["stderr"].splitlines():
            m = SKETCH_RE.match(line)
            if m:
                sketch_s = int(m.group(1)) / 1000.0
                continue
            m = PAIR_RE.match(line)
            if m:
                pair_s = int(m.group(1)) / 1000.0
        r["sketch_creation_time"] = sketch_s
        r["comparison_time"] = pair_s

        if keep_output:
            csvs = [f for f in glob.glob(out_prefix + "*") if f.endswith(".csv")]
            r["output_csv"] = csvs[0] if csvs else None
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


def run_benchmark(
    binary: str,
    num_files_list: List[int],
    num_runs: int,
    num_intervals: int,
    num_threads: int,
    precision: int,
    sub_b_list: List[float],
) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    print("\nBenchmark configuration:")
    print(f"  hammock-cpp binary: {binary}")
    print(f"  intervals/file:     {num_intervals}")
    print(f"  file counts:        {num_files_list}")
    print(f"  runs/config:        {num_runs}")
    print(f"  threads:            {num_threads}")
    print(f"  HLL precision:      {precision}")
    print(f"  subB values:        {sub_b_list}")
    print(f"  system:             {get_system_info()}")

    metric_keys = ["wall_time", "cpu_time", "max_rss_mb", "sort_time"]
    hammock_keys = metric_keys + ["sketch_creation_time", "comparison_time"]

    for num_files in num_files_list:
        print(f"\n{'=' * 60}\n{num_files} files × {num_intervals} intervals\n{'=' * 60}")
        bedtools_runs: List[Dict[str, Any]] = []
        hammock_runs_by_subb: Dict[float, List[Dict[str, Any]]] = {s: [] for s in sub_b_list}

        for run_i in range(num_runs):
            print(f"\nRun {run_i + 1}/{num_runs}: generating {2 * num_files} BED files...")
            with tempfile.TemporaryDirectory() as tmp_dir:
                file1_list, file2_list = [], []
                for i in range(num_files):
                    a = os.path.join(tmp_dir, f"set1_{i}.bed")
                    b = os.path.join(tmp_dir, f"set2_{i}.bed")
                    generate_bed_file(num_intervals, a)
                    generate_bed_file(num_intervals, b)
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
                print("  bedtools...", end=" ", flush=True)
                bt = run_bedtools(file1_list_path, file2_list_path, num_threads)
                bt["sort_time"] = sort_time
                bedtools_runs.append(bt)
                rss = bt["max_rss_mb"]
                rss_str = f"{rss:.1f} MB" if rss is not None else "n/a"
                print(f"{bt['wall_time']:.2f}s wall, {rss_str} max RSS")

                for sub_b in sub_b_list:
                    label = f"subB={sub_b:g}" if sub_b < 1.0 else "subB=1.0"
                    print(f"  hammock-cpp Mode B {label}...", end=" ", flush=True)
                    hm = run_hammock(binary, file1_list_path, file2_list_path,
                                     precision, num_threads, sub_b=sub_b)
                    hm["sort_time"] = sort_time
                    hammock_runs_by_subb[sub_b].append(hm)
                    rss = hm["max_rss_mb"]
                    rss_str = f"{rss:.1f} MB" if rss is not None else "n/a"
                    print(f"{hm['wall_time']:.2f}s wall, {rss_str} max RSS")

        entry: Dict[str, Any] = {
            "num_files": num_files,
            "num_intervals_per_file": num_intervals,
            "num_threads": num_threads,
            "precision": precision,
            "sub_b_list": sub_b_list,
            "bedtools": aggregate(bedtools_runs, metric_keys),
        }
        for sub_b, runs_for in hammock_runs_by_subb.items():
            entry[tool_name_for_subb(sub_b)] = aggregate(runs_for, hammock_keys)
        results.append(entry)

    return results


def write_text_report(results: List[Dict[str, Any]], path: str) -> None:
    with open(path, "w") as f:
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
            sub_b_list = r.get("sub_b_list", [1.0])
            f.write(f"num_files={r['num_files']}\n")
            f.write("-" * 80 + "\n")
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
            for sub_b in sub_b_list:
                tool = tool_name_for_subb(sub_b)
                hm = r[tool]
                f.write(f"hammock-cpp Mode B [subB={sub_b:g}, mixed-stride]:\n")
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


def write_csv(results: List[Dict[str, Any]], path: str) -> None:
    cols = [
        "num_files", "num_threads", "precision", "sub_b", "tool",
        "mean_wall_time", "std_wall_time", "min_wall_time", "max_wall_time",
        "mean_cpu_time", "std_cpu_time", "min_cpu_time", "max_cpu_time",
        "mean_max_rss_mb", "std_max_rss_mb", "min_max_rss_mb", "max_max_rss_mb",
        "mean_sort_time", "std_sort_time", "min_sort_time", "max_sort_time",
        "mean_sketch_creation_time", "std_sketch_creation_time",
        "mean_comparison_time", "std_comparison_time",
    ]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in results:
            sub_b_list = r.get("sub_b_list", [1.0])
            tools = [("bedtools", "bedtools", "")]
            for sub_b in sub_b_list:
                tools.append((tool_name_for_subb(sub_b), tool_name_for_subb(sub_b), f"{sub_b:g}"))
            for tool, key, sub_b_str in tools:
                d = r[key]
                row = [r["num_files"], r["num_threads"], r["precision"], sub_b_str, tool]
                for c in cols[5:]:
                    v = d.get(c)
                    row.append(f"{v:.6f}" if isinstance(v, float) else ("" if v is None else v))
                w.writerow(row)


def plot(results: List[Dict[str, Any]], png_path: str) -> None:
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        print("matplotlib not available, skipping plot")
        return

    nfiles = [r["num_files"] for r in results]
    bt_wall = [r["bedtools"]["mean_wall_time"] for r in results]
    bt_cpu = [r["bedtools"]["mean_cpu_time"] for r in results]
    bt_rss = [r["bedtools"]["mean_max_rss_mb"] for r in results]
    threads = results[0]["num_threads"] if results else 1
    sub_b_list = results[0].get("sub_b_list", [1.0]) if results else [1.0]

    hm_colors = ["#ff7f0e", "#1f77b4", "#2ca02c", "#d62728"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].plot(nfiles, bt_wall, "ko-", label=f"bedtools (t={threads})", linewidth=2)
    for i, sub_b in enumerate(sub_b_list):
        tool = tool_name_for_subb(sub_b)
        ys = [r[tool]["mean_wall_time"] for r in results]
        axes[0].plot(nfiles, ys, "s--", color=hm_colors[i % len(hm_colors)],
                     label=f"hammock-cpp subB={sub_b:g}", linewidth=2)
    axes[0].set_xlabel("Number of files")
    axes[0].set_ylabel("Wall time (s)")
    axes[0].set_title("Wall time")
    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(nfiles, bt_cpu, "ko-", label="bedtools", linewidth=2)
    for i, sub_b in enumerate(sub_b_list):
        tool = tool_name_for_subb(sub_b)
        ys = [r[tool]["mean_cpu_time"] for r in results]
        axes[1].plot(nfiles, ys, "s--", color=hm_colors[i % len(hm_colors)],
                     label=f"hammock-cpp subB={sub_b:g}", linewidth=2)
    axes[1].set_xlabel("Number of files")
    axes[1].set_ylabel("CPU time (s)")
    axes[1].set_title("CPU time")
    axes[1].set_xscale("log", base=2)
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    hm_rss_any = any(r[tool_name_for_subb(s)]["mean_max_rss_mb"] is not None
                     for r in results for s in sub_b_list)
    if any(v is not None for v in bt_rss) or hm_rss_any:
        axes[2].plot(nfiles, bt_rss, "ko-", label="bedtools", linewidth=2)
        for i, sub_b in enumerate(sub_b_list):
            tool = tool_name_for_subb(sub_b)
            ys = [r[tool]["mean_max_rss_mb"] for r in results]
            axes[2].plot(nfiles, ys, "s--", color=hm_colors[i % len(hm_colors)],
                         label=f"hammock-cpp subB={sub_b:g}", linewidth=2)
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
    parser.add_argument("--test", action="store_true", help="Quick test (small files, few runs)")
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

    results = run_benchmark(
        binary=binary,
        num_files_list=num_files_list,
        num_runs=num_runs,
        num_intervals=num_intervals,
        num_threads=args.threads,
        precision=args.precision,
        sub_b_list=sub_b_list,
    )

    write_text_report(results, txt_path)
    write_csv(results, csv_path)
    plot(results, png_path)
    print(f"\nReports: {txt_path}\n         {csv_path}\nFigure:  {png_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
