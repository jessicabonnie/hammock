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
emits the full 9-column block, matching hammock-cpp's default since v0.7.0),
so hammock-cpp here is run WITH --metrics for an apples-to-apples comparison
of the same output shape.

Small and quick by design (order-of-magnitude characterization for a
limitations sentence, not a headline figure): a handful of N points, one
subB, one precision, 3 replicates, tool order alternated.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from benchmark_cpp_vs_bedtools import (  # noqa: E402
    generate_bed_file, derive_seed, get_system_info, check_binary_version,
)


def run_with_wall(cmd: List[str]) -> Dict[str, Any]:
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed (exit {r.returncode}): {r.stderr[-2000:]}")
    return {"wall_time": wall}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--num-files", default="2,8,32,128")
    ap.add_argument("--num-intervals", type=int, default=10000)
    ap.add_argument("--precision", type=int, default=18)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--corpus-seed", type=int, default=20260811)
    ap.add_argument("--hammock-cpp-bin", default=os.environ.get(
        "HAMMOCK_CPP_BIN",
        "/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp"))
    ap.add_argument("--hammock-cli-bin", default=
        "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    for b in (args.hammock_cpp_bin, args.hammock_cli_bin):
        if not os.path.exists(b):
            print(f"not found: {b}", file=sys.stderr)
            return 1
    check_binary_version(args.hammock_cpp_bin)

    print(f"hammock-cpp: {args.hammock_cpp_bin}")
    print(f"hammock CLI: {args.hammock_cli_bin}")
    print(f"system: {get_system_info()}")

    num_files_list = [int(x) for x in args.num_files.split(",")]
    rows = []

    for num_files in num_files_list:
        print(f"\n{'=' * 60}\n{num_files} files x {args.num_intervals} intervals\n{'=' * 60}")
        for run_i in range(args.runs):
            with tempfile.TemporaryDirectory() as tmp_dir:
                file1_list, file2_list = [], []
                for i in range(num_files):
                    a = os.path.join(tmp_dir, f"set1_{i}.bed")
                    b = os.path.join(tmp_dir, f"set2_{i}.bed")
                    generate_bed_file(args.num_intervals, a, derive_seed(args.corpus_seed, run_i, i, 0))
                    generate_bed_file(args.num_intervals, b, derive_seed(args.corpus_seed, run_i, i, 1))
                    file1_list.append(a)
                    file2_list.append(b)
                # sort (both tools require it identically; not timed)
                for p in file1_list + file2_list:
                    subprocess.run(["sort", "-k1,1", "-k2,2n", "-o", p, p], check=True)
                l1 = os.path.join(tmp_dir, "l1.txt")
                l2 = os.path.join(tmp_dir, "l2.txt")
                with open(l1, "w") as f:
                    f.write("\n".join(file1_list) + "\n")
                with open(l2, "w") as f:
                    f.write("\n".join(file2_list) + "\n")

                cpp_out = os.path.join(tmp_dir, "cpp_out")
                cli_out = os.path.join(tmp_dir, "cli_out")
                cpp_cmd = [args.hammock_cpp_bin, l1, l2, "--mode", "B",
                           "-p", str(args.precision), "-o", cpp_out,
                           "--threads", str(args.threads), "--metrics"]
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
                                 "tool": label, "wall_time": r["wall_time"]})

    out_path = args.output or os.path.join(
        SCRIPT_DIR, "results", f"cli_overhead_{int(time.time())}.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    import csv
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["num_files", "run", "tool", "wall_time"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out_path}")

    # quick summary
    import statistics
    from collections import defaultdict
    by_cell = defaultdict(list)
    for r in rows:
        by_cell[(r["num_files"], r["tool"])].append(r["wall_time"])
    print(f"\n{'N':>6} {'cpp_median':>12} {'cli_median':>12} {'cli/cpp':>10}")
    for n in num_files_list:
        cpp = statistics.median(by_cell[(n, "hammock_cpp")])
        cli = statistics.median(by_cell[(n, "hammock_cli")])
        print(f"{n:>6} {cpp:>12.4f} {cli:>12.4f} {cli / cpp:>10.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
