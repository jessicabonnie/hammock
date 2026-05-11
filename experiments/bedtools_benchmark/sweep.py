#!/usr/bin/env python3
"""
Multi-axis sweep driver for hammock-cpp Mode B vs bedtools.

One axis per invocation (--axis precision|threads|intervals). Generates BED data
once per (num_files, num_intervals) combo within a run and reuses it across the
swept variable, so noise from data regeneration doesn't pollute the comparison.

Outputs a long-format CSV plus a two-panel PNG (wall time + max RSS) into
results/ and figures/.
"""

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from datetime import datetime

import numpy as np  # type: ignore

from benchmark_cpp_vs_bedtools import (
    BEDTOOLS_SCRIPT,
    FIGURES_DIR,
    RESULTS_DIR,
    find_hammock_cpp,
    generate_bed_file,
    get_system_info,
    run_bedtools,
    run_hammock,
)


METRIC_KEYS = ["wall_time", "cpu_time", "max_rss_mb"]
HAMMOCK_KEYS = METRIC_KEYS + ["sketch_creation_time", "comparison_time"]
ACCURACY_KEYS = [
    "jaccard_n_pairs",
    "jaccard_mae_vs_bt", "jaccard_max_err_vs_bt",
    "jaccard_mae_vs_hll", "jaccard_max_err_vs_hll",
]
# sort_time: wall time to pre-sort BED files for the (run_id, num_files, num_intervals)
# realization. Pre-sort happens outside per-tool timing (bedtools needs sorted input,
# hammock is indifferent). Same sort_time is attached to every row from the same
# realization. See docs/bedtools-parallelism-caveat.md.
ROW_COLS = [
    "axis", "tool", "precision", "threads", "num_files", "num_intervals", "run_id",
    "wall_time", "cpu_time", "max_rss_mb", "sort_time",
    "sketch_creation_time", "comparison_time",
] + ACCURACY_KEYS


def parse_bedtools_jaccards(stdout: str):
    """bedtools.sh emits one line per pair: <file_a>\\t<file_b>\\t<jaccard>.

    Returns {(basename_a, basename_b): jaccard_float}.
    """
    out = {}
    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        a, b, j = parts
        try:
            out[(os.path.basename(a), os.path.basename(b))] = float(j)
        except ValueError:
            continue
    return out


def parse_hammock_csv(path: str):
    """Hammock CSV is query\\treference\\tjaccard_similarity (basenames).

    Returns {(query, reference): jaccard_float}.
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as f:
        header = f.readline()  # skip
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                out[(parts[0], parts[1])] = float(parts[2])
            except ValueError:
                continue
    return out


def jaccard_error_stats(est: dict, *, vs_bt: dict, vs_hll: dict):
    """Compare hammock estimates against TWO ground truths.

    vs_bt: bedtools per-pair jaccards (set-jaccard ground truth — definitional gap)
    vs_hll: hammock@p_max per-pair jaccards (HLL ground truth — actual precision/accuracy)
    """
    out = {k: None for k in ACCURACY_KEYS}
    common_bt = set(est.keys()) & set(vs_bt.keys())
    common_hll = set(est.keys()) & set(vs_hll.keys()) if vs_hll else set()
    out["jaccard_n_pairs"] = len(common_bt) if common_bt else (len(common_hll) or None)
    if common_bt:
        d = np.array([abs(est[k] - vs_bt[k]) for k in common_bt])
        out["jaccard_mae_vs_bt"] = float(np.mean(d))
        out["jaccard_max_err_vs_bt"] = float(np.max(d))
    if common_hll:
        d = np.array([abs(est[k] - vs_hll[k]) for k in common_hll])
        out["jaccard_mae_vs_hll"] = float(np.mean(d))
        out["jaccard_max_err_vs_hll"] = float(np.max(d))
    return out


def _sort_one(path: str) -> None:
    """Sort a single BED file in-place by chrom,start."""
    sorted_path = path + ".sorted"
    with open(sorted_path, "w") as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", path], stdout=out, check=True)
    os.rename(sorted_path, path)


def make_data(num_files: int, num_intervals: int, tmp_dir: str, num_sort_workers: int = 8):
    """Generate and pre-sort 2*num_files BED files.

    Returns (file1_list_path, file2_list_path, sort_time). sort_time covers only
    the parallel sort step (not generation). Sort is fanned out across
    num_sort_workers threads — each thread blocks on its own external `sort`
    subprocess, so the GIL isn't in play. Pre-sort is outside per-tool timing
    because bedtools requires sorted input — see
    docs/bedtools-parallelism-caveat.md for the fairness framing.
    """
    file1, file2 = [], []
    for i in range(num_files):
        a = os.path.join(tmp_dir, f"set1_{i}.bed")
        b = os.path.join(tmp_dir, f"set2_{i}.bed")
        generate_bed_file(num_intervals, a)
        generate_bed_file(num_intervals, b)
        file1.append(a)
        file2.append(b)
    sort_start = time.time()
    all_paths = file1 + file2
    workers = max(1, min(num_sort_workers, len(all_paths)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        for _ in ex.map(_sort_one, all_paths):
            pass
    sort_time = time.time() - sort_start
    f1 = os.path.join(tmp_dir, "file1_list.txt")
    f2 = os.path.join(tmp_dir, "file2_list.txt")
    with open(f1, "w") as f:
        f.write("\n".join(file1) + "\n")
    with open(f2, "w") as f:
        f.write("\n".join(file2) + "\n")
    return f1, f2, sort_time


def _row(axis, tool, *, precision, threads, num_files, num_intervals, run_id, result, keys):
    r = {
        "axis": axis, "tool": tool,
        "precision": precision, "threads": threads,
        "num_files": num_files, "num_intervals": num_intervals,
        "run_id": run_id,
    }
    for k in keys:
        r[k] = result.get(k)
    return r


def sweep_precision(binary, precisions, num_files, num_intervals, num_threads, num_runs):
    """Vary precision; bedtools is precision-independent so it's run once per data realization.

    Computes jaccard accuracy with TWO ground truths:
      - bedtools (set-jaccard): exposes the definitional gap
      - hammock at max(precisions) (HLL ground truth): actual HLL precision/accuracy curve

    Also dumps every per-pair jaccard for downstream scatter plots.
    """
    import glob as _glob
    rows = []
    pair_rows = []  # for the scatter plot: one row per (run, query, reference, precision)
    p_max = max(precisions)

    for run_i in range(num_runs):
        with tempfile.TemporaryDirectory() as tmp_dir:
            print(f"\n[precision] run {run_i+1}/{num_runs}: gen {num_files} files × {num_intervals} intervals")
            f1, f2, sort_time = make_data(num_files, num_intervals, tmp_dir, num_sort_workers=num_threads)
            print(f"  sort:    {sort_time:.2f}s wall (parallel, {num_threads} workers; pre-sort, not in tool wall below)")

            print("  bedtools...", end=" ", flush=True)
            bt = run_bedtools(f1, f2, num_threads)
            bt["sort_time"] = sort_time
            rss = bt["max_rss_mb"]
            print(f"{bt['wall_time']:.2f}s wall, {rss:.1f} MB" if rss is not None else f"{bt['wall_time']:.2f}s wall")
            bt_truth = parse_bedtools_jaccards(bt.get("stdout", ""))
            print(f"    captured {len(bt_truth)} bedtools jaccards (set-jaccard ground truth)")
            rows.append(_row("precision", "bedtools",
                             precision=None, threads=num_threads,
                             num_files=num_files, num_intervals=num_intervals,
                             run_id=run_i, result=bt, keys=METRIC_KEYS + ["sort_time"]))

            # Run all precisions; collect per-pair estimates first so we can use
            # the largest-p result as the HLL ground truth for the others.
            estimates_by_p = {}
            hm_results_by_p = {}
            for p in precisions:
                print(f"  hammock-cpp p={p}...", end=" ", flush=True)
                hm = run_hammock(binary, f1, f2, p, num_threads, keep_output=True)
                hm["sort_time"] = sort_time
                rss = hm["max_rss_mb"]
                print(f"{hm['wall_time']:.2f}s wall, {rss:.1f} MB" if rss is not None else f"{hm['wall_time']:.2f}s wall")
                estimates_by_p[p] = parse_hammock_csv(hm.get("output_csv"))
                hm_results_by_p[p] = hm
                prefix = hm.get("_out_prefix")
                if prefix:
                    for f in _glob.glob(prefix + "*"):
                        try:
                            os.remove(f)
                        except OSError:
                            pass

            hll_truth = estimates_by_p.get(p_max, {})

            for p in precisions:
                est = estimates_by_p[p]
                err = jaccard_error_stats(est, vs_bt=bt_truth, vs_hll=hll_truth)
                print(f"    p={p}: n={err['jaccard_n_pairs']}, MAE_bt={err['jaccard_mae_vs_bt']:.4f}, MAE_hll={err['jaccard_mae_vs_hll']:.4f}")
                merged = dict(hm_results_by_p[p])
                merged.update(err)
                rows.append(_row("precision", "hammock_cpp_B",
                                 precision=p, threads=num_threads,
                                 num_files=num_files, num_intervals=num_intervals,
                                 run_id=run_i, result=merged, keys=HAMMOCK_KEYS + ["sort_time"] + ACCURACY_KEYS))
                for (q, r), j_hm in est.items():
                    pair_rows.append({
                        "run_id": run_i, "precision": p,
                        "query": q, "reference": r,
                        "hammock_jaccard": j_hm,
                        "bedtools_jaccard": bt_truth.get((q, r)),
                    })
    return rows, pair_rows


def sweep_threads(binary, thread_list, precision, num_files, num_intervals, num_runs):
    rows = []
    for run_i in range(num_runs):
        with tempfile.TemporaryDirectory() as tmp_dir:
            sort_workers = max(thread_list)
            print(f"\n[threads] run {run_i+1}/{num_runs}: gen {num_files} files × {num_intervals} intervals")
            f1, f2, sort_time = make_data(num_files, num_intervals, tmp_dir, num_sort_workers=sort_workers)
            print(f"  sort:    {sort_time:.2f}s wall (parallel, {sort_workers} workers; pre-sort, not in tool wall below)")

            for t in thread_list:
                print(f"  t={t} bedtools...", end=" ", flush=True)
                bt = run_bedtools(f1, f2, t)
                bt["sort_time"] = sort_time
                print(f"{bt['wall_time']:.2f}s wall")
                rows.append(_row("threads", "bedtools",
                                 precision=None, threads=t,
                                 num_files=num_files, num_intervals=num_intervals,
                                 run_id=run_i, result=bt, keys=METRIC_KEYS + ["sort_time"]))

                print(f"  t={t} hammock-cpp p={precision}...", end=" ", flush=True)
                hm = run_hammock(binary, f1, f2, precision, t)
                hm["sort_time"] = sort_time
                print(f"{hm['wall_time']:.2f}s wall")
                rows.append(_row("threads", "hammock_cpp_B",
                                 precision=precision, threads=t,
                                 num_files=num_files, num_intervals=num_intervals,
                                 run_id=run_i, result=hm, keys=HAMMOCK_KEYS + ["sort_time"]))
    return rows


def sweep_intervals(binary, intervals_list, precision, num_threads, num_files, num_runs):
    """Regenerates data per intervals point — file size is the swept variable."""
    rows = []
    for ni in intervals_list:
        for run_i in range(num_runs):
            with tempfile.TemporaryDirectory() as tmp_dir:
                print(f"\n[intervals] N={ni} × {num_files} files: run {run_i+1}/{num_runs}")
                f1, f2, sort_time = make_data(num_files, ni, tmp_dir, num_sort_workers=num_threads)
                print(f"  sort:    {sort_time:.2f}s wall (parallel, {num_threads} workers; pre-sort, not in tool wall below)")

                print("  bedtools...", end=" ", flush=True)
                bt = run_bedtools(f1, f2, num_threads)
                bt["sort_time"] = sort_time
                print(f"{bt['wall_time']:.2f}s wall")
                rows.append(_row("intervals", "bedtools",
                                 precision=None, threads=num_threads,
                                 num_files=num_files, num_intervals=ni,
                                 run_id=run_i, result=bt, keys=METRIC_KEYS + ["sort_time"]))

                print("  hammock-cpp...", end=" ", flush=True)
                hm = run_hammock(binary, f1, f2, precision, num_threads)
                hm["sort_time"] = sort_time
                print(f"{hm['wall_time']:.2f}s wall")
                rows.append(_row("intervals", "hammock_cpp_B",
                                 precision=precision, threads=num_threads,
                                 num_files=num_files, num_intervals=ni,
                                 run_id=run_i, result=hm, keys=HAMMOCK_KEYS + ["sort_time"]))
    return rows


def write_long_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            out = dict(r)
            for k, v in out.items():
                if isinstance(v, float):
                    out[k] = f"{v:.6f}"
                elif v is None:
                    out[k] = ""
            w.writerow(out)


def _agg(values_by_x, key):
    xs = sorted(values_by_x.keys())
    means, stds = [], []
    for x in xs:
        vals = [v[key] for v in values_by_x[x] if v.get(key) is not None]
        if vals:
            means.append(float(np.mean(vals)))
            stds.append(float(np.std(vals)) if len(vals) > 1 else 0.0)
        else:
            means.append(None)
            stds.append(None)
    return xs, means, stds


def plot_axis(rows, axis, png_path, title_suffix):
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        print("matplotlib not available, skipping plot", file=sys.stderr)
        return

    if axis == "precision":
        x_key, x_label, x_log = "precision", "HLL precision (p)", None
    elif axis == "threads":
        x_key, x_label, x_log = "threads", "Thread count", ("log", 2)
    else:
        x_key, x_label, x_log = "num_intervals", "Intervals per file", ("log", 10)

    bt = defaultdict(list)
    hm = defaultdict(list)
    for r in rows:
        x = r.get(x_key)
        if r["tool"] == "hammock_cpp_B" and x is not None:
            hm[x].append(r)
        elif r["tool"] == "bedtools" and (x is not None or axis == "precision"):
            bt[x if x is not None else "ref"].append(r)

    n_panels = 3 if axis == "precision" else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(6.5 * n_panels, 5))

    # Wall time
    if hm:
        xs, ys, errs = _agg(hm, "wall_time")
        axes[0].errorbar(xs, ys, yerr=errs, marker="s", linestyle="--",
                         color="#ff7f0e", label="hammock-cpp Mode B")
    if axis == "precision":
        bt_walls = [r["wall_time"] for r in rows if r["tool"] == "bedtools" and r.get("wall_time") is not None]
        if bt_walls:
            axes[0].axhline(np.mean(bt_walls), color="k", linestyle="-",
                            label=f"bedtools (avg)")
    else:
        if bt:
            xs, ys, errs = _agg(bt, "wall_time")
            axes[0].errorbar(xs, ys, yerr=errs, marker="o", linestyle="-",
                             color="k", label="bedtools")
    axes[0].set_xlabel(x_label)
    axes[0].set_ylabel("Wall time (s)")
    axes[0].set_title(f"Wall time {title_suffix}")
    if x_log:
        axes[0].set_xscale(x_log[0], base=x_log[1])
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    # Max RSS
    if hm:
        xs, ys, _ = _agg(hm, "max_rss_mb")
        ys_clean = [(x, y) for x, y in zip(xs, ys) if y is not None]
        if ys_clean:
            xs2, ys2 = zip(*ys_clean)
            axes[1].plot(xs2, ys2, "s--", color="#ff7f0e", label="hammock-cpp Mode B")
    if axis == "precision":
        bt_rss = [r["max_rss_mb"] for r in rows if r["tool"] == "bedtools" and r.get("max_rss_mb") is not None]
        if bt_rss:
            axes[1].axhline(np.mean(bt_rss), color="k", linestyle="-", label="bedtools (avg)")
    else:
        if bt:
            xs, ys, _ = _agg(bt, "max_rss_mb")
            ys_clean = [(x, y) for x, y in zip(xs, ys) if y is not None]
            if ys_clean:
                xs2, ys2 = zip(*ys_clean)
                axes[1].plot(xs2, ys2, "o-", color="k", label="bedtools")
    axes[1].set_xlabel(x_label)
    axes[1].set_ylabel("Max RSS (MB)")
    axes[1].set_title(f"Peak memory {title_suffix}")
    if x_log:
        axes[1].set_xscale(x_log[0], base=x_log[1])
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    # Accuracy panel (precision sweep only) — two ground truths.
    if axis == "precision" and hm:
        xs, ys_hll, errs_hll = _agg(hm, "jaccard_mae_vs_hll")
        valid_hll = [(x, y, e) for x, y, e in zip(xs, ys_hll, errs_hll) if y is not None]
        if valid_hll:
            xs2, ys2, errs2 = zip(*valid_hll)
            axes[2].errorbar(xs2, ys2, yerr=errs2, marker="s", linestyle="--",
                             color="#ff7f0e", label="MAE vs hammock@p_max (HLL noise)")
            xs_th = list(xs2)
            ys_th = [1.04 / np.sqrt(2 ** x) for x in xs_th]
            axes[2].plot(xs_th, ys_th, "k:", alpha=0.6, label="1.04/√(2^p) (theory)")

        xs, ys_bt, errs_bt = _agg(hm, "jaccard_mae_vs_bt")
        valid_bt = [(x, y, e) for x, y, e in zip(xs, ys_bt, errs_bt) if y is not None]
        if valid_bt:
            xs2, ys2, errs2 = zip(*valid_bt)
            axes[2].errorbar(xs2, ys2, yerr=errs2, marker="o", linestyle="-",
                             color="#1f77b4", label="MAE vs bedtools (definitional gap)")
        axes[2].set_xlabel(x_label)
        axes[2].set_ylabel("Jaccard MAE")
        axes[2].set_title(f"Accuracy {title_suffix}")
        axes[2].set_yscale("log")
        axes[2].grid(True, alpha=0.3)
        axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_pair_scatter(pair_rows, png_path, title_suffix):
    """Scatter of (bedtools_jaccard, hammock_jaccard), one point per pair, colored by precision.

    Shows the empirical bp-jaccard ↔ register-equality-jaccard mapping. y=x diagonal
    is reference; if hammock matched bedtools they'd lie on it. Realistic-input
    deviations from the diagonal are the definitional gap (see
    docs/jaccard-definitional-gap.md).
    """
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        return

    by_p = defaultdict(list)
    for r in pair_rows:
        if r.get("bedtools_jaccard") is None or r.get("hammock_jaccard") is None:
            continue
        by_p[r["precision"]].append((r["bedtools_jaccard"], r["hammock_jaccard"]))

    if not by_p:
        return

    fig, ax = plt.subplots(figsize=(7, 7))
    cmap = plt.get_cmap("viridis")
    ps = sorted(by_p.keys())
    for i, p in enumerate(ps):
        xs, ys = zip(*by_p[p])
        ax.scatter(xs, ys, s=8, alpha=0.4,
                   color=cmap(i / max(1, len(ps) - 1)),
                   label=f"p={p}")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5, label="y = x (set-jaccard equality)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("bedtools jaccard (set Jaccard, ground truth)")
    ax.set_ylabel("hammock jaccard (register-equality)")
    ax.set_title(f"Definitional gap: hammock vs bedtools jaccard\n{title_suffix}", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--axis", choices=["precision", "threads", "intervals"], required=True)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--binary", type=str, default=None)
    parser.add_argument("--output-dir", type=str, default=RESULTS_DIR)
    parser.add_argument("--figures-dir", type=str, default=FIGURES_DIR)

    parser.add_argument("--precisions", type=str, default="10,12,14,16,18")
    parser.add_argument("--thread-list", type=str, default="1,2,4,8,16")
    parser.add_argument("--intervals-list", type=str, default="1000,10000,100000,1000000")

    parser.add_argument("--precision", type=int, default=14)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--num-files", type=int, default=None,
                        help="Default: 64 (precision), 64 (threads), 16 (intervals)")
    parser.add_argument("--num-intervals", type=int, default=10000)

    args = parser.parse_args()
    binary = args.binary or find_hammock_cpp()
    if not os.path.exists(binary):
        print(f"hammock-cpp not found at {binary}", file=sys.stderr)
        return 1
    if not os.path.exists(BEDTOOLS_SCRIPT):
        print(f"bedtools.sh not found at {BEDTOOLS_SCRIPT}", file=sys.stderr)
        return 1

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.figures_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    print(f"system: {get_system_info()}")
    print(f"binary: {binary}")

    pair_rows = []
    if args.axis == "precision":
        precisions = [int(x) for x in args.precisions.split(",")]
        nf = args.num_files if args.num_files is not None else 64
        print(f"[precision] p ∈ {precisions}, t={args.threads}, files={nf}, intervals={args.num_intervals}")
        rows, pair_rows = sweep_precision(binary, precisions, nf, args.num_intervals, args.threads, args.runs)
        suffix = f"(t={args.threads}, files={nf}, intervals={args.num_intervals})"
    elif args.axis == "threads":
        thread_list = [int(x) for x in args.thread_list.split(",")]
        nf = args.num_files if args.num_files is not None else 64
        print(f"[threads] t ∈ {thread_list}, p={args.precision}, files={nf}, intervals={args.num_intervals}")
        rows = sweep_threads(binary, thread_list, args.precision, nf, args.num_intervals, args.runs)
        suffix = f"(p={args.precision}, files={nf}, intervals={args.num_intervals})"
    else:
        intervals_list = [int(x) for x in args.intervals_list.split(",")]
        nf = args.num_files if args.num_files is not None else 16
        print(f"[intervals] N ∈ {intervals_list}, p={args.precision}, t={args.threads}, files={nf}")
        rows = sweep_intervals(binary, intervals_list, args.precision, args.threads, nf, args.runs)
        suffix = f"(p={args.precision}, t={args.threads}, files={nf})"

    stem = f"sweep_{args.axis}_{timestamp}"
    csv_path = os.path.join(args.output_dir, stem + ".csv")
    png_path = os.path.join(args.figures_dir, stem + ".png")
    write_long_csv(rows, csv_path)
    plot_axis(rows, args.axis, png_path, suffix)
    print(f"\nResults: {csv_path}\nFigure:  {png_path}")

    if pair_rows:
        pairs_csv = os.path.join(args.output_dir, stem + "_pairs.csv")
        with open(pairs_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["run_id", "precision", "query", "reference",
                                              "bedtools_jaccard", "hammock_jaccard"],
                               extrasaction="ignore")
            w.writeheader()
            for r in pair_rows:
                w.writerow(r)
        scatter_png = os.path.join(args.figures_dir, stem + "_scatter.png")
        plot_pair_scatter(pair_rows, scatter_png, suffix)
        print(f"Pairs:   {pairs_csv}\nScatter: {scatter_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
