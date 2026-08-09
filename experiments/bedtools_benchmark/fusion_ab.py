#!/usr/bin/env python3
"""Paired A/B of two hammock-cpp builds x two output shapes, in one allocation.

Why this exists. The cost of the similarity block was previously inferred by
comparing a 2026-08-04 run to a 2026-08-08 run -- different binary, different
node, different corpus, and the older run had no SLURM allocation at all. That
comparison is not sound: the `--no-metrics` arm runs byte-identical code in both,
and it measured 1.27-1.53x SLOWER on the old run, so any absolute cross-run
number carries an unknown machine factor of that size.

The design here follows standard practice for comparing two builds:

  * PAIRED. All four arms run against the same corpus in the same job on the
    same node, so a machine factor is common-mode and cancels in the ratio.
    Reported quantities are within-replicate ratios, never cross-run absolutes.
  * INTERLEAVED AND ORDER-RANDOMIZED. The four arms are permuted per replicate
    rather than run in blocks. Blocked execution confounds arm with position --
    which is exactly what the Maurano leg did (two invocations 12 minutes apart)
    and why its wall-time deltas came out with the wrong sign at 9 sigma.
  * SEEDED CORPUS, so a re-run measures the same inputs.
  * EFFECT SIZE + INTERVAL, not a point estimate: the per-replicate ratios are
    reported with their spread, so a difference smaller than the noise is
    visible as such.

Not duet-style (both binaries running simultaneously). Duet equalizes
interference by construction and is the better choice when spare cores exist,
but both arms here are OpenMP jobs sized to the whole allocation, so running
them concurrently would halve each and measure contention instead of the change.

The four arms:
    pre  x --no-metrics     <- unchanged code in both builds; the control
    pre  x --metrics        <- the old union-allocating path
    post x --no-metrics     <- control again
    post x --metrics        <- the fused path

The control pair is the point of the design: pre/post on --no-metrics must come
out at 1.00, because that code is identical. Any departure is the measurement
error of the whole experiment, and it bounds how much of the --metrics effect
can be believed.

Usage:
    python3 fusion_ab.py --pre-binary <path> --post-binary <path> --runs 5
"""

import argparse
import csv
import hashlib
import itertools
import os
import statistics as st
import sys
import tempfile
from datetime import datetime

from benchmark_cpp_vs_bedtools import (
    derive_seed,
    generate_bed_file,
    get_system_info,
    run_hammock,
    RESULTS_DIR,
)

DEFAULT_PRECISIONS = [12, 14, 16, 18, 20, 22, 24]
DEFAULT_NUM_FILES = 64
DEFAULT_NUM_INTERVALS = 10000
DEFAULT_THREADS = 16
DEFAULT_RUNS = 5
DEFAULT_SEED = 20260808

ROW_COLS = [
    "precision", "num_files", "build", "arm", "run_index", "order_position",
    "threads", "num_intervals", "n_pairs",
    "sketch_time", "comparison_time", "pair_time", "write_time",
    "wall_time", "cpu_time", "max_rss_mb", "us_per_pair",
    "binary_md5", "hostname", "cpu_model", "cpu_count", "slurm_job_id",
    "git_sha", "corpus_seed",
]


def md5_of(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def make_corpus(num_files: int, num_intervals: int, tmp_dir: str, seed: int):
    file1, file2 = [], []
    for i in range(num_files):
        a = os.path.join(tmp_dir, f"set1_{i}.bed")
        b = os.path.join(tmp_dir, f"set2_{i}.bed")
        generate_bed_file(num_intervals, a, derive_seed(seed, i, 0))
        generate_bed_file(num_intervals, b, derive_seed(seed, i, 1))
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
    ap.add_argument("--pre-binary", required=True)
    ap.add_argument("--post-binary", required=True)
    ap.add_argument("--precisions", default=",".join(map(str, DEFAULT_PRECISIONS)))
    ap.add_argument("--num-files", type=int, default=DEFAULT_NUM_FILES)
    ap.add_argument("--num-intervals", type=int, default=DEFAULT_NUM_INTERVALS)
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    ap.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    ap.add_argument("--corpus-seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--output-dir", default=RESULTS_DIR)
    args = ap.parse_args()

    precisions = [int(x) for x in args.precisions.split(",")]
    builds = {"pre": args.pre_binary, "post": args.post_binary}
    for name, path in builds.items():
        if not os.path.exists(path):
            print(f"{name} binary not found: {path}", file=sys.stderr)
            return 1
    md5s = {k: md5_of(v) for k, v in builds.items()}
    if md5s["pre"] == md5s["post"]:
        # Would silently produce a null result that looks like "no effect".
        print("pre and post binaries are byte-identical; nothing to compare",
              file=sys.stderr)
        return 1

    n_pairs = args.num_files * args.num_files
    sysinfo = get_system_info()
    print(f"pre  : {args.pre_binary}  (md5 {md5s['pre']})")
    print(f"post : {args.post_binary} (md5 {md5s['post']})")
    print(f"precisions: {precisions}   files/side: {args.num_files}   "
          f"threads: {args.threads}   runs: {args.runs}")
    print(f"corpus seed: {args.corpus_seed}")
    print(f"system: {sysinfo}")
    if sysinfo.get("slurm_job_id", "none") == "none":
        print("  WARNING: no SLURM allocation; cores are shared and the control "
              "arm will show it.")

    # All four arms, permuted per replicate. itertools.permutations is used
    # rather than random.shuffle so the schedule is deterministic and printed --
    # setup randomization only helps if it is reproducible.
    arms = [(b, m) for b in ("pre", "post") for m in (False, True)]
    perms = list(itertools.permutations(range(len(arms))))

    os.makedirs(args.output_dir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(args.output_dir, f"fusion_ab_{stamp}.csv")

    rows = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        print(f"\ngenerating {2 * args.num_files} BED files (seed {args.corpus_seed})...",
              flush=True)
        f1, f2 = make_corpus(args.num_files, args.num_intervals, tmp_dir,
                             args.corpus_seed)

        for p in precisions:
            for run_i in range(args.runs):
                order = perms[run_i % len(perms)]
                for pos, arm_idx in enumerate(order):
                    build, use_metrics = arms[arm_idx]
                    r = run_hammock(builds[build], f1, f2, p, args.threads,
                                    metrics=use_metrics)
                    pt = r["pair_time"]
                    rows.append({
                        "precision": p, "num_files": args.num_files,
                        "build": build,
                        "arm": "metrics" if use_metrics else "no_metrics",
                        "run_index": run_i, "order_position": pos,
                        "threads": args.threads,
                        "num_intervals": args.num_intervals, "n_pairs": n_pairs,
                        "sketch_time": r["sketch_creation_time"],
                        "comparison_time": r["comparison_time"],
                        "pair_time": pt, "write_time": r["write_time"],
                        "wall_time": r["wall_time"], "cpu_time": r["cpu_time"],
                        "max_rss_mb": r["max_rss_mb"],
                        "us_per_pair": (pt * 1e6 / n_pairs) if pt is not None else None,
                        "binary_md5": md5s[build],
                        "hostname": sysinfo.get("hostname", "unknown"),
                        "cpu_model": sysinfo.get("cpu_model", "unknown"),
                        "cpu_count": sysinfo.get("cpu_count", "unknown"),
                        "slurm_job_id": sysinfo.get("slurm_job_id", "none"),
                        "git_sha": sysinfo.get("git_sha", "unknown"),
                        "corpus_seed": args.corpus_seed,
                    })
                    print(f"  p={p:2d} run={run_i} pos={pos} {build:4s} "
                          f"{'metrics   ' if use_metrics else 'no_metrics'} "
                          f"pair={pt:8.4f}s wall={r['wall_time']:7.2f}s",
                          flush=True)
                # Flush per replicate: a crash at p=24 would otherwise discard
                # every completed measurement in the job.
                with open(out_path, "w", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=ROW_COLS)
                    w.writeheader()
                    w.writerows(rows)

    print(f"\nWrote {out_path}")

    def cell(p, build, arm, run_i, col="pair_time"):
        for r in rows:
            if (r["precision"] == p and r["build"] == build and r["arm"] == arm
                    and r["run_index"] == run_i and r[col] is not None):
                return r[col]
        return None

    def paired(p, col, num, den):
        """Ratios formed WITHIN a replicate, then summarised across replicates."""
        out = []
        for run_i in range(args.runs):
            a = cell(p, num[0], num[1], run_i, col)
            b = cell(p, den[0], den[1], run_i, col)
            if a and b:
                out.append(a / b)
        return out

    def fmt(v):
        if not v:
            return "     n/a"
        m = st.median(v)
        lo, hi = min(v), max(v)
        return f"{m:5.2f} [{lo:.2f}-{hi:.2f}]"

    print("\nCONTROL -- post/pre on --no-metrics. Identical code: must be ~1.00.")
    print("Its departure from 1.00 is this experiment's measurement error.")
    print(f"  {'p':>3}  {'pair_time':>20}  {'wall_time':>20}")
    for p in precisions:
        print(f"  {p:>3}  "
              f"{fmt(paired(p,'pair_time',('post','no_metrics'),('pre','no_metrics'))):>20}  "
              f"{fmt(paired(p,'wall_time',('post','no_metrics'),('pre','no_metrics'))):>20}")

    print("\nEFFECT -- post/pre on --metrics (the fusion speedup of the block):")
    print(f"  {'p':>3}  {'pair_time':>20}  {'comparison_time':>20}")
    for p in precisions:
        print(f"  {p:>3}  "
              f"{fmt(paired(p,'pair_time',('post','metrics'),('pre','metrics'))):>20}  "
              f"{fmt(paired(p,'comparison_time',('post','metrics'),('pre','metrics'))):>20}")

    print("\nBLOCK COST -- metrics/no_metrics within each build:")
    print(f"  {'p':>3}  {'pre pair_time':>20}  {'post pair_time':>20}")
    for p in precisions:
        print(f"  {p:>3}  "
              f"{fmt(paired(p,'pair_time',('pre','metrics'),('pre','no_metrics'))):>20}  "
              f"{fmt(paired(p,'pair_time',('post','metrics'),('post','no_metrics'))):>20}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
