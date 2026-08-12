#!/usr/bin/env python3
"""Isolation test for docs/seed-hammock-cpp-file-dispatch.md Part 1.

Single file per side -> no file-dispatch (ThreadPoolExecutor) parallelism is
possible for either front-end, so any CPU usage beyond ~100% * threads can
only come from the sketching phase's own OpenMP team ignoring --threads.

Writes a CSV so this claim has a committed artifact, matching the standard
the seed doc itself demands of every other quantitative claim in it.
"""
import argparse
import csv
import os
import random
import re
import subprocess
import sys
import tempfile
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from benchmark_cpp_vs_bedtools import get_system_info, check_binary_version  # noqa: E402

CPU_RE = re.compile(r"Percent of CPU this job got:\s+(\d+)%")
WALL_RE = re.compile(r"Elapsed \(wall clock\) time.*:\s+([\d:.]+)")
USER_RE = re.compile(r"User time \(seconds\):\s+([\d.]+)")


def parse_time_v(stderr: str):
    cpu = CPU_RE.search(stderr)
    user = USER_RE.search(stderr)
    return {
        "cpu_pct": int(cpu.group(1)) if cpu else None,
        "user_s": float(user.group(1)) if user else None,
    }


def run_timed(cmd):
    r = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed: {r.stderr[-1500:]}")
    return parse_time_v(r.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threads-list", default="2,4,16")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--num-intervals", type=int, default=200000)
    ap.add_argument("--precision", type=int, default=18)
    ap.add_argument("--hammock-cpp-bin", default=os.environ.get(
        "HAMMOCK_CPP_BIN",
        "/home/jbonnie1/interval_sketch/hammock_claude/build/cp310-cp310-linux_x86_64/hammock-cpp"))
    ap.add_argument("--hammock-cli-bin", default=
        "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    check_binary_version(args.hammock_cpp_bin)
    sysinfo = get_system_info()
    print(f"system: {sysinfo}")
    print(f"hammock-cpp: {args.hammock_cpp_bin}")
    print(f"hammock CLI: {args.hammock_cli_bin}")

    threads_list = [int(x) for x in args.threads_list.split(",")]
    random.seed(20260811)
    rows = []

    with tempfile.TemporaryDirectory() as tmp:
        a = os.path.join(tmp, "a.bed")
        b = os.path.join(tmp, "b.bed")
        for name in (a, b):
            with open(name, "w") as f:
                lines = []
                for _ in range(args.num_intervals):
                    chrom = f"chr{random.randint(1, 5)}"
                    s = random.randint(0, 100_000_000)
                    lines.append(f"{chrom}\t{s}\t{s + 150}\n")
                f.writelines(lines)
            subprocess.run(["sort", "-k1,1", "-k2,2n", "-o", name, name], check=True)
        l1 = os.path.join(tmp, "l1.txt")
        l2 = os.path.join(tmp, "l2.txt")
        with open(l1, "w") as f:
            f.write(a + "\n")
        with open(l2, "w") as f:
            f.write(b + "\n")

        for t in threads_list:
            for rep in range(args.reps):
                cpp_out = os.path.join(tmp, f"cpp_t{t}_r{rep}")
                cli_out = os.path.join(tmp, f"cli_t{t}_r{rep}")
                # PART9 (docs/seed-metrics-column-restructure.md): --metrics
                # here vs. cli_cmd's bare invocation below is a shape
                # mismatch (_full vs _ie) under the new three-shape contract.
                # Neither arm reads a CSV column here (only /usr/bin/time
                # output is parsed), so the fix is to drop --metrics from
                # cpp_cmd to match cli_cmd's bare shape; not retargeted here
                # per "mark PART9 only, don't retarget."
                cpp_cmd = [args.hammock_cpp_bin, l1, l2, "--mode", "B",
                           "-p", str(args.precision), "-o", cpp_out,
                           "--threads", str(t), "--metrics"]
                cli_cmd = [args.hammock_cli_bin, l1, l2, "--mode", "B",
                           "--precision", str(args.precision), "-o", cli_out,
                           "--threads", str(t)]
                # alternate order across reps to avoid a systematic first-vs-second confound
                pair = [("hammock_cpp", cpp_cmd), ("hammock_cli", cli_cmd)]
                if rep % 2:
                    pair = pair[::-1]
                for label, cmd in pair:
                    m = run_timed(cmd)
                    print(f"  threads={t} rep={rep} {label}: cpu={m['cpu_pct']}% user={m['user_s']}s")
                    rows.append({"threads": t, "rep": rep, "tool": label,
                                 "cpu_pct": m["cpu_pct"], "user_s": m["user_s"]})

    out_path = args.output or os.path.join(
        SCRIPT_DIR, "results", f"threads_isolation_{int(time.time())}.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["threads", "rep", "tool", "cpu_pct", "user_s"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out_path}")

    import statistics
    from collections import defaultdict
    by_cell = defaultdict(list)
    for r in rows:
        by_cell[(r["threads"], r["tool"])].append(r["cpu_pct"])
    print(f"\n{'threads':>8} {'cpp_median_cpu%':>16} {'cli_median_cpu%':>16}")
    for t in threads_list:
        cpp = statistics.median(by_cell[(t, "hammock_cpp")])
        cli = statistics.median(by_cell[(t, "hammock_cli")])
        print(f"{t:>8} {cpp:>16} {cli:>16}")


if __name__ == "__main__":
    main()
