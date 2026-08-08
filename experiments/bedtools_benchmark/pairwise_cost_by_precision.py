#!/usr/bin/env python3
"""What the similarity block costs, as a function of HyperLogLog precision.

Since v0.7.0 `hammock-cpp` emits `jaccard_similarity_ie` and the
containment/cosketch columns by default, and `--no-metrics` is the opt-out used
by every timed benchmark. This measures the difference between those two arms so
Methods can state the cost of the configuration the paper recommends rather than
assert it.

Precision is the axis, not file count. The extra work is a union sketch plus two
cardinality estimates per pair, all of it inside the pairwise phase, so the
per-pair cost is independent of N once the timer resolves it (measured: µs/pair
flattens above N≈16) while it grows steeply with 2^p. Fixing N and sweeping p
therefore isolates the quantity of interest at a fraction of the runtime -- a
p=24 point at N=512 spends ~13 min per run in *sketching* alone.

No bedtools here: this is a hammock-vs-hammock comparison and the bedtools leg
would be ~83% of the wall time while contributing nothing.

Two measurement notes that matter for how the numbers are read:

  * `us_per_pair` is `pair_time / (n*m)` at a FIXED thread count. It is a
    throughput figure, not a serial per-pair cost, and it does not transfer
    across thread counts. The thread count is recorded in every row; state it
    wherever the number is quoted.
  * `pair_time` excludes the serial `%.17g` write loop, which `write_time`
    carries separately. That split is the point: six extra output fields per row
    are precision-independent, so at low p they are most of the difference
    between the arms and would otherwise be charged to the estimator.

Arms alternate by run index so arm order is not confounded with position.

Usage:
    python3 pairwise_cost_by_precision.py --runs 5
    python3 pairwise_cost_by_precision.py --precisions 14,18,24 --runs 2   # quick
"""

import argparse
import csv
import os
import sys
import tempfile
from datetime import datetime

from benchmark_cpp_vs_bedtools import (
    check_binary_version,
    find_hammock_cpp,
    generate_bed_file,
    get_system_info,
    run_hammock,
    RESULTS_DIR,
)

DEFAULT_PRECISIONS = [12, 14, 16, 18, 20, 22, 24]
DEFAULT_NUM_FILES = 64        # per side, so 2N sketches and N*N pairs
DEFAULT_NUM_INTERVALS = 10000  # matches the synthetic scaling benchmark
DEFAULT_THREADS = 16           # pinned: the login node's OpenMP default is 48
DEFAULT_RUNS = 5

# Provenance travels in the rows, not just on stdout. The 20260804 CSV recorded
# none of it, so it cannot be checked for comparability after the fact: with
# `-march=native` baked in, a timing is only comparable to one from the same CPU
# model, and a run outside a SLURM allocation shares its cores with whatever else
# is on the node. Repeating five constant columns across 70 rows is a trivial
# price for a self-describing file that still joins cleanly.
PROVENANCE_COLS = [
    "hostname", "cpu_model", "cpu_count", "slurm_job_id", "binary_version",
]

ROW_COLS = [
    "precision", "num_files", "arm", "run_index", "threads", "num_intervals",
    "n_pairs", "sketch_time", "comparison_time", "pair_time", "write_time",
    "us_per_pair",
] + PROVENANCE_COLS


def make_corpus(num_files: int, num_intervals: int, tmp_dir: str):
    """Two disjoint lists of num_files BEDs each, as the scaling benchmark builds.

    No pre-sort: hammock is indifferent to input order and there is no bedtools
    leg here to require it.
    """
    file1, file2 = [], []
    for i in range(num_files):
        a = os.path.join(tmp_dir, f"set1_{i}.bed")
        b = os.path.join(tmp_dir, f"set2_{i}.bed")
        generate_bed_file(num_intervals, a)
        generate_bed_file(num_intervals, b)
        file1.append(a)
        file2.append(b)
    f1 = os.path.join(tmp_dir, "file1_list.txt")
    f2 = os.path.join(tmp_dir, "file2_list.txt")
    with open(f1, "w") as f:
        f.write("\n".join(file1) + "\n")
    with open(f2, "w") as f:
        f.write("\n".join(file2) + "\n")
    return f1, f2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--precisions", default=",".join(map(str, DEFAULT_PRECISIONS)))
    ap.add_argument("--num-files", type=int, default=DEFAULT_NUM_FILES)
    ap.add_argument("--num-intervals", type=int, default=DEFAULT_NUM_INTERVALS)
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    ap.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    ap.add_argument("--binary", default=None)
    ap.add_argument("--output-dir", default=RESULTS_DIR)
    args = ap.parse_args()

    precisions = [int(x) for x in args.precisions.split(",")]
    binary = args.binary or find_hammock_cpp()
    version = check_binary_version(binary)

    n_pairs = args.num_files * args.num_files
    print(f"hammock-cpp:   {binary} ({version})")
    print(f"precisions:    {precisions}")
    print(f"files/side:    {args.num_files}  ({2 * args.num_files} sketches, {n_pairs:,} pairs)")
    print(f"intervals:     {args.num_intervals}")
    print(f"threads:       {args.threads}")
    print(f"runs:          {args.runs}")
    sysinfo = get_system_info()
    print(f"system:        {sysinfo}")
    provenance = {
        "hostname": sysinfo.get("hostname", "unknown"),
        "cpu_model": sysinfo.get("cpu_model", "unknown"),
        "cpu_count": sysinfo.get("cpu_count", "unknown"),
        "slurm_job_id": sysinfo.get("slurm_job_id", "none"),
        "binary_version": version,
    }
    if provenance["slurm_job_id"] == "none":
        print("  NOTE: no SLURM allocation -- cores are shared with anything else "
              "on this node, which inflates wall times unpredictably.")

    os.makedirs(args.output_dir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(args.output_dir, f"pairwise_cost_by_precision_{stamp}.csv")

    rows = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        print(f"\ngenerating {2 * args.num_files} BED files...", flush=True)
        f1, f2 = make_corpus(args.num_files, args.num_intervals, tmp_dir)

        for p in precisions:
            for run_i in range(args.runs):
                # Alternate, so "metrics" is not always the second (warmer) arm.
                arms = [("no_metrics", False), ("metrics", True)]
                if run_i % 2:
                    arms.reverse()
                for arm, use_metrics in arms:
                    r = run_hammock(binary, f1, f2, p, args.threads,
                                    metrics=use_metrics)
                    pair_time = r["pair_time"]
                    us_per_pair = (pair_time * 1e6 / n_pairs) if pair_time else None
                    rows.append({
                        **provenance,
                        "precision": p,
                        "num_files": args.num_files,
                        "arm": arm,
                        "run_index": run_i,
                        "threads": args.threads,
                        "num_intervals": args.num_intervals,
                        "n_pairs": n_pairs,
                        "sketch_time": r["sketch_creation_time"],
                        "comparison_time": r["comparison_time"],
                        "pair_time": pair_time,
                        "write_time": r["write_time"],
                        "us_per_pair": us_per_pair,
                    })
                    print(f"  p={p:2d} run={run_i} {arm:10s} "
                          f"sketch={r['sketch_creation_time']:8.3f}s "
                          f"pair={pair_time:8.4f}s write={r['write_time']:7.4f}s "
                          f"{us_per_pair:9.2f} us/pair", flush=True)

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_path}")

    # Median per (precision, arm), so the summary is robust to one slow run.
    print("\nMedian us/pair by precision:")
    print(f"  {'p':>3}  {'no_metrics':>12}  {'metrics':>12}  {'ratio':>7}")
    for p in precisions:
        med = {}
        for arm in ("no_metrics", "metrics"):
            vals = sorted(r["us_per_pair"] for r in rows
                          if r["precision"] == p and r["arm"] == arm
                          and r["us_per_pair"] is not None)
            med[arm] = vals[len(vals) // 2] if vals else None
        if med["no_metrics"] and med["metrics"]:
            print(f"  {p:>3}  {med['no_metrics']:>12.2f}  {med['metrics']:>12.2f}  "
                  f"{med['metrics'] / med['no_metrics']:>7.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
