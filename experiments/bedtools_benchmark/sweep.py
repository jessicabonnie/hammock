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
    check_binary_version,
    run_bedtools,
    run_hammock,
    tool_name_for_subb,
)


METRIC_KEYS = ["wall_time", "cpu_time", "max_rss_mb"]
HAMMOCK_KEYS = METRIC_KEYS + ["sketch_creation_time", "comparison_time",
                              "pair_time", "write_time"]
ACCURACY_KEYS = [
    "jaccard_n_pairs",
    "jaccard_mae_vs_bt", "jaccard_max_err_vs_bt",
    "jaccard_mae_vs_hll", "jaccard_max_err_vs_hll",
    # Inclusion-exclusion twin. `jaccard_similarity` is register-equality,
    # which carries a chance-agreement floor and is NOT on bedtools' scale --
    # its _vs_bt columns measure the definitional gap, not error. The _ie_
    # columns are the ones comparable to bedtools. Left None when the sweep
    # ran with --no-ie. See docs/jaccard-definitional-gap.md.
    "jaccard_ie_mae_vs_bt", "jaccard_ie_max_err_vs_bt",
    "jaccard_ie_mae_vs_hll", "jaccard_ie_max_err_vs_hll",
]
# Timing of the --metrics invocation that supplies the jaccard_ie_* columns.
# Kept apart from METRIC_KEYS so it can never be mistaken for the timed arm.
METRICS_ARM_KEYS = ["metrics_wall_time", "metrics_cpu_time", "metrics_max_rss_mb"]
# Provenance travels in the rows, not just on stdout, mirroring
# pairwise_cost_by_precision.py's PROVENANCE_COLS. Before this the file recorded
# none of it: `get_system_info()` was printed once and then discarded, so a
# finished CSV could not be checked for comparability against any other. That
# matters for the p=18 figure specifically, because the two jobs behind it are
# two separate SLURM allocations and may land on different nodes -- and with
# `-march=native` baked into the binary, a different node can mean different
# code. Repeating six constant columns is a trivial price for a self-describing
# file.
PROVENANCE_COLS = [
    "hostname", "cpu_model", "cpu_count", "slurm_job_id", "binary_version",
    # binary_version cannot identify the code on its own: pyproject's version is
    # bumped per release, not per commit, so a pre- and post-change binary both
    # report the same string.
    "git_sha",
]

# sort_time: wall time to pre-sort BED files for the (run_id, num_files, num_intervals)
# realization. Pre-sort happens outside per-tool timing (bedtools needs sorted input,
# hammock is indifferent). Same sort_time is attached to every row from the same
# realization. See docs/bedtools-parallelism-caveat.md. On a fixed corpus that
# ships pre-sorted it is 0.0 -- "sorted before the experiment began", the same
# starting condition for both tools, not "sorting was free".
ROW_COLS = [
    "axis", "tool", "sub_b", "precision", "threads", "num_files", "num_intervals", "run_id",
    # Which BED corpus the row came from: "synthetic" (generated per replicate)
    # or "maurano" (20 real fetal DHS files, fixed across replicates). Speed and
    # accuracy are only jointly interpretable within one corpus.
    "corpus",
    "wall_time", "cpu_time", "max_rss_mb", "sort_time",
    # comparison_time keeps its historical meaning (pair loop + serial write);
    # pair_time/write_time decompose it and are blank on bedtools rows.
    "sketch_creation_time", "comparison_time", "pair_time", "write_time",
    # Wall time of the SECOND, --metrics invocation on the precision axis --
    # the one that produces the jaccard_ie_* columns. It used to be printed as
    # "(not recorded)" and thrown away, which meant a frontier plot would show
    # the *speed* of the 3-column arm against the *accuracy* of the 9-column
    # arm. Blank when the sweep ran --no-ie, and on every non-precision axis.
    # Not comparable to `wall_time` as a substitute: it is a different output
    # shape. It is there to answer "what does computing the x-axis cost?".
] + METRICS_ARM_KEYS + ACCURACY_KEYS + PROVENANCE_COLS


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


def parse_hammock_csv(path: str, column: str = "jaccard_similarity"):
    """Read one similarity column out of a hammock-cpp TSV, by header name.

    Returns {(query, reference): value}. Column lookup is by name, not
    position: a --no-metrics file has 3 columns and the default has 9, so an
    index would silently read containment_AB as if it were a Jaccard.

    Raises KeyError if the file exists but lacks `column` -- that means the
    binary ran with --no-metrics, and returning {} instead would show up
    downstream as "no pairs in common", i.e. a silently empty accuracy column
    rather than a failure.
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            idx = header.index(column)
        except ValueError:
            raise KeyError(
                f"{path} has no {column!r} column (header: {header}). "
                f"Re-run hammock-cpp without --no-metrics.") from None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= idx:
                continue
            try:
                out[(parts[0], parts[1])] = float(parts[idx])
            except ValueError:
                continue
    return out


def _is_self_pair(key) -> bool:
    """True when a (query, reference) key compares a file with itself.

    Only reachable when the same list is passed as both operands, which is the
    fixed-corpus convention (`--corpus maurano`, mirroring
    experiments/subB_mixed_stride/run_sweep.py). The synthetic path builds two
    disjoint lists, so nothing there matches.
    """
    q, r = key
    return os.path.basename(q) == os.path.basename(r)


def jaccard_error_stats(est: dict, *, vs_bt: dict, vs_hll: dict, prefix: str = "jaccard",
                        drop_self_pairs: bool = False):
    """Compare hammock estimates against TWO ground truths.

    vs_bt: bedtools per-pair jaccards (set-jaccard ground truth — definitional gap)
    vs_hll: hammock@p_max per-pair jaccards (HLL ground truth — actual precision/accuracy)

    prefix selects the key family: "jaccard" for register-equality,
    "jaccard_ie" for inclusion-exclusion. The vs_hll ground truth must be
    drawn from the SAME estimator as `est` -- comparing IE estimates against a
    register-equality p_max reference would report the definitional gap as if
    it were precision error.

    drop_self_pairs excludes file-vs-itself comparisons. With one 20-file list
    passed twice, 20 of the 400 ordered pairs are self-comparisons where both
    tools return ~1.0 at zero error -- free correctness that dilutes MAE by
    ~5% (1.152e-3 -> 1.094e-3 at p=18 on Maurano) without measuring anything.
    Off by default so the synthetic path, where the two lists are disjoint and
    no key can match, keeps producing byte-identical numbers to the archive.
    """
    out = {}
    if drop_self_pairs:
        est = {k: v for k, v in est.items() if not _is_self_pair(k)}
    common_bt = set(est.keys()) & set(vs_bt.keys())
    common_hll = set(est.keys()) & set(vs_hll.keys()) if vs_hll else set()
    if prefix == "jaccard":
        out["jaccard_n_pairs"] = len(common_bt) if common_bt else (len(common_hll) or None)
    if common_bt:
        d = np.array([abs(est[k] - vs_bt[k]) for k in common_bt])
        out[f"{prefix}_mae_vs_bt"] = float(np.mean(d))
        out[f"{prefix}_max_err_vs_bt"] = float(np.max(d))
    if common_hll:
        d = np.array([abs(est[k] - vs_hll[k]) for k in common_hll])
        out[f"{prefix}_mae_vs_hll"] = float(np.mean(d))
        out[f"{prefix}_max_err_vs_hll"] = float(np.max(d))
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


# Maurano fetal-tissue DHS corpus, already sorted on disk by
# experiments/subB_mixed_stride/prepare_maurano.sh. This file lives in
# <repo>/experiments/bedtools_benchmark/, so the sibling experiment is two
# levels up and back down.
_SWEEP_DIR = os.path.dirname(os.path.abspath(__file__))
MAURANO_DIR = os.path.join(os.path.dirname(_SWEEP_DIR),
                           "subB_mixed_stride", "data", "maurano_sorted")


def make_maurano_data(tmp_dir: str):
    """Point the sweep at the 20 real Maurano BEDs instead of generating a corpus.

    Returns the same (file1_list_path, file2_list_path, sort_time) triple as
    make_data, so sweep_precision does not care which it got.

    Two deliberate differences from the synthetic path:

      * **The same list is passed as both operands**, following
        experiments/subB_mixed_stride/run_sweep.py. That yields 400 ordered
        pairs of which 20 are self-comparisons -- see `drop_self_pairs` in
        jaccard_error_stats for why those are excluded from MAE but still
        timed. They are legitimately part of the timed workload; they are just
        not evidence about accuracy.
      * **sort_time is 0.0, not a measurement.** These files ship pre-sorted,
        so there is no sort to time. It is 0.0 rather than None so the column
        stays numeric, but do not read it as "sorting was free" -- it is "the
        corpus was sorted before the experiment began", which is the same
        starting condition for both tools.

    The corpus is NOT copied into tmp_dir: 20 files at ~100 Mbp each is ~2 GB,
    and copying would put the page-cache state of the copy, not of the shared
    read-only originals, in front of whichever tool runs first.
    """
    beds = sorted(f for f in os.listdir(MAURANO_DIR) if f.endswith(".bed"))
    if not beds:
        raise FileNotFoundError(
            f"no .bed files in {MAURANO_DIR}; run "
            f"experiments/subB_mixed_stride/prepare_maurano.sh first")
    paths = [os.path.join(MAURANO_DIR, b) for b in beds]
    f1 = os.path.join(tmp_dir, "file1_list.txt")
    with open(f1, "w") as f:
        f.write("\n".join(paths) + "\n")
    return f1, f1, 0.0


def _row(axis, tool, *, precision, threads, num_files, num_intervals, run_id, result, keys,
         sub_b=None, corpus="synthetic", provenance=None):
    r = {
        "axis": axis, "tool": tool, "sub_b": sub_b,
        "precision": precision, "threads": threads,
        "num_files": num_files, "num_intervals": num_intervals,
        "run_id": run_id, "corpus": corpus,
    }
    if provenance:
        r.update(provenance)
    for k in keys:
        r[k] = result.get(k)
    return r


def sweep_precision(binary, precisions, num_files, num_intervals, num_threads, num_runs,
                    sub_b_list, ie=True, corpus="synthetic", provenance=None):
    """Vary precision; bedtools is precision-independent so it's run once per data realization.

    Computes jaccard accuracy with TWO ground truths:
      - bedtools (set-jaccard): exposes the definitional gap
      - hammock at subB=1.0, max(precisions) (HLL ground truth): actual HLL
        precision/accuracy curve. For subB<1.0 runs, this same ground truth
        is used, so MAE_hll captures both HLL noise AND sampling noise.

    Also dumps every per-pair jaccard for downstream scatter plots.
    """
    import glob as _glob
    rows = []
    pair_rows = []  # for the scatter plot: one row per (run, query, reference, precision, sub_b)
    p_max = max(precisions)

    drop_self = corpus != "synthetic"

    for run_i in range(num_runs):
        with tempfile.TemporaryDirectory() as tmp_dir:
            if corpus == "maurano":
                print(f"\n[precision] run {run_i+1}/{num_runs}: fixed Maurano corpus "
                      f"({num_files} files, same list both sides)")
                f1, f2, sort_time = make_maurano_data(tmp_dir)
                print(f"  sort:    n/a (corpus ships pre-sorted)")
            else:
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
            if not bt_truth:
                # An empty parse is the silent-failure mode this whole sweep is
                # most exposed to: every MAE column would come out blank and the
                # run would still "succeed" after ~12 minutes of compute.
                raise RuntimeError(
                    "bedtools produced no parseable jaccards -- every accuracy "
                    "column would be empty. Check bedtools.sh and that the "
                    "inputs are sorted.")
            rows.append(_row("precision", "bedtools",
                             precision=None, threads=num_threads,
                             num_files=num_files, num_intervals=num_intervals,
                             run_id=run_i, result=bt, keys=METRIC_KEYS + ["sort_time"],
                             corpus=corpus, provenance=provenance))

            # Run all (precision, subB) combos; collect per-pair estimates first
            # so subB=1.0 @ p_max can serve as the HLL ground truth for all
            # other (p, subB) combos.
            estimates: dict = {}
            estimates_ie: dict = {}
            hm_results: dict = {}
            metrics_timing: dict = {}

            def _drop(hm):
                prefix = hm.get("_out_prefix")
                if not prefix:
                    return
                for f in _glob.glob(prefix + "*"):
                    try:
                        os.remove(f)
                    except OSError:
                        pass

            for sub_b in sub_b_list:
                for p in precisions:
                    print(f"  hammock-cpp p={p} subB={sub_b:g}...", end=" ", flush=True)
                    hm = run_hammock(binary, f1, f2, p, num_threads,
                                     keep_output=True, sub_b=sub_b)
                    hm["sort_time"] = sort_time
                    rss = hm["max_rss_mb"]
                    print(f"{hm['wall_time']:.2f}s wall, {rss:.1f} MB" if rss is not None else f"{hm['wall_time']:.2f}s wall")
                    estimates[(p, sub_b)] = parse_hammock_csv(hm.get("output_csv"))
                    hm_results[(p, sub_b)] = hm
                    _drop(hm)

                    if ie:
                        # SECOND pass, in its OWN columns. --metrics adds a
                        # cardinality estimate per pair, so folding it into the
                        # run above would inflate comparison_time and break
                        # comparability with the published RESULTS.md numbers.
                        # Since 0.7.0 the timed pass above is explicitly
                        # --no-metrics, so the split is belt-and-braces rather
                        # than the only thing keeping the timings honest.
                        #
                        # What changed: this used to print "(not recorded)" and
                        # throw the timing away. That made a precision-frontier
                        # plot dishonest by construction -- the x-axis (IE MAE)
                        # comes from THIS invocation while the y-axis (speed)
                        # came from the other one, so the figure would report
                        # the accuracy of an arm nobody timed at the speed of an
                        # arm that does not compute it. The numbers go to
                        # metrics_* columns, never into wall_time, so no
                        # archived series changes meaning.
                        print(f"    + --metrics pass...", end=" ", flush=True)
                        hm2 = run_hammock(binary, f1, f2, p, num_threads,
                                          keep_output=True, sub_b=sub_b, metrics=True)
                        estimates_ie[(p, sub_b)] = parse_hammock_csv(
                            hm2.get("output_csv"), column="jaccard_similarity_ie")
                        metrics_timing[(p, sub_b)] = {
                            "metrics_wall_time": hm2.get("wall_time"),
                            "metrics_cpu_time": hm2.get("cpu_time"),
                            "metrics_max_rss_mb": hm2.get("max_rss_mb"),
                        }
                        _drop(hm2)
                        ratio = (hm2["wall_time"] / hm["wall_time"]
                                 if hm.get("wall_time") else None)
                        print(f"{hm2['wall_time']:.2f}s wall"
                              + (f" ({ratio:.2f}x the --no-metrics arm)" if ratio else ""))

            hll_truth = estimates.get((p_max, 1.0), {})
            # IE ground truth must be the IE value at p_max, not the
            # register-equality one -- they estimate different quantities.
            hll_truth_ie = estimates_ie.get((p_max, 1.0), {})

            for sub_b in sub_b_list:
                for p in precisions:
                    est = estimates[(p, sub_b)]
                    est_ie = estimates_ie.get((p, sub_b), {})
                    err = {k: None for k in ACCURACY_KEYS}
                    err.update(jaccard_error_stats(est, vs_bt=bt_truth, vs_hll=hll_truth,
                                                   drop_self_pairs=drop_self))
                    if est_ie:
                        err.update(jaccard_error_stats(
                            est_ie, vs_bt=bt_truth, vs_hll=hll_truth_ie,
                            prefix="jaccard_ie", drop_self_pairs=drop_self))

                    def _f(k):
                        v = err[k]
                        return f"{v:.4f}" if v is not None else "n/a"
                    print(f"    p={p} subB={sub_b:g}: n={err['jaccard_n_pairs']}, "
                          f"MAE_bt={_f('jaccard_mae_vs_bt')} (register-equality, "
                          f"includes the definitional gap), "
                          f"MAE_bt_ie={_f('jaccard_ie_mae_vs_bt')}, "
                          f"MAE_hll={_f('jaccard_mae_vs_hll')}")
                    merged = dict(hm_results[(p, sub_b)])
                    merged.update(err)
                    merged.update(metrics_timing.get((p, sub_b), {}))
                    rows.append(_row("precision", tool_name_for_subb(sub_b),
                                     precision=p, threads=num_threads,
                                     num_files=num_files, num_intervals=num_intervals,
                                     run_id=run_i, result=merged,
                                     keys=HAMMOCK_KEYS + ["sort_time"] + ACCURACY_KEYS
                                          + METRICS_ARM_KEYS,
                                     sub_b=sub_b, corpus=corpus, provenance=provenance))
                    for (q, r), j_hm in est.items():
                        pair_rows.append({
                            "run_id": run_i, "precision": p, "sub_b": sub_b,
                            "query": q, "reference": r,
                            "hammock_jaccard": j_hm,
                            "hammock_jaccard_ie": est_ie.get((q, r)),
                            "bedtools_jaccard": bt_truth.get((q, r)),
                        })
    return rows, pair_rows


def sweep_threads(binary, thread_list, precision, num_files, num_intervals, num_runs,
                  sub_b_list, provenance=None):
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
                                 run_id=run_i, result=bt, keys=METRIC_KEYS + ["sort_time"],
                                 provenance=provenance))

                for sub_b in sub_b_list:
                    print(f"  t={t} hammock-cpp p={precision} subB={sub_b:g}...", end=" ", flush=True)
                    hm = run_hammock(binary, f1, f2, precision, t, sub_b=sub_b)
                    hm["sort_time"] = sort_time
                    print(f"{hm['wall_time']:.2f}s wall")
                    rows.append(_row("threads", tool_name_for_subb(sub_b),
                                     precision=precision, threads=t,
                                     num_files=num_files, num_intervals=num_intervals,
                                     run_id=run_i, result=hm,
                                     keys=HAMMOCK_KEYS + ["sort_time"], sub_b=sub_b,
                                     provenance=provenance))
    return rows


def sweep_intervals(binary, intervals_list, precision, num_threads, num_files, num_runs,
                    sub_b_list, provenance=None):
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
                                 run_id=run_i, result=bt, keys=METRIC_KEYS + ["sort_time"],
                                 provenance=provenance))

                for sub_b in sub_b_list:
                    print(f"  hammock-cpp subB={sub_b:g}...", end=" ", flush=True)
                    hm = run_hammock(binary, f1, f2, precision, num_threads, sub_b=sub_b)
                    hm["sort_time"] = sort_time
                    print(f"{hm['wall_time']:.2f}s wall")
                    rows.append(_row("intervals", tool_name_for_subb(sub_b),
                                     precision=precision, threads=num_threads,
                                     num_files=num_files, num_intervals=ni,
                                     run_id=run_i, result=hm,
                                     keys=HAMMOCK_KEYS + ["sort_time"], sub_b=sub_b,
                                     provenance=provenance))
    return rows


def write_long_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            out = dict(r)
            for k, v in out.items():
                if isinstance(v, float):
                    # %.9g, not %.6f. The old format gave six *decimal places*,
                    # not six significant figures, so an accuracy column was
                    # written at whatever precision its magnitude happened to
                    # allow: at p=24 an IE MAE near 1.4e-4 landed as 0.000144 --
                    # three significant figures for the quantity the precision
                    # frontier is plotted against, and anything under 5e-7
                    # rounded to a hard 0.000000, indistinguishable from a real
                    # zero. Timing columns are unaffected in practice (they are
                    # O(1-100 s), where the two formats agree to more digits
                    # than the measurement carries).
                    out[k] = f"{v:.9g}"
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
    hm_by_subb: dict = defaultdict(lambda: defaultdict(list))  # sub_b -> x -> [rows]
    for r in rows:
        x = r.get(x_key)
        if r["tool"].startswith("hammock_cpp_B") and x is not None:
            sub_b = r.get("sub_b") if r.get("sub_b") is not None else 1.0
            hm_by_subb[float(sub_b)][x].append(r)
        elif r["tool"] == "bedtools" and (x is not None or axis == "precision"):
            bt[x if x is not None else "ref"].append(r)
    sub_bs = sorted(hm_by_subb.keys(), reverse=True)  # 1.0 first
    hm_colors = ["#ff7f0e", "#1f77b4", "#2ca02c", "#d62728"]

    n_panels = 3 if axis == "precision" else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(6.5 * n_panels, 5))

    # Wall time
    for i, sub_b in enumerate(sub_bs):
        hm = hm_by_subb[sub_b]
        if not hm:
            continue
        xs, ys, errs = _agg(hm, "wall_time")
        axes[0].errorbar(xs, ys, yerr=errs, marker="s", linestyle="--",
                         color=hm_colors[i % len(hm_colors)],
                         label=f"hammock-cpp subB={sub_b:g}")
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
    for i, sub_b in enumerate(sub_bs):
        hm = hm_by_subb[sub_b]
        if not hm:
            continue
        xs, ys, _ = _agg(hm, "max_rss_mb")
        ys_clean = [(x, y) for x, y in zip(xs, ys) if y is not None]
        if ys_clean:
            xs2, ys2 = zip(*ys_clean)
            axes[1].plot(xs2, ys2, "s--", color=hm_colors[i % len(hm_colors)],
                         label=f"hammock-cpp subB={sub_b:g}")
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

    # Accuracy panel (precision sweep only) — two ground truths per subB.
    if axis == "precision" and hm_by_subb:
        for i, sub_b in enumerate(sub_bs):
            hm = hm_by_subb[sub_b]
            if not hm:
                continue
            color = hm_colors[i % len(hm_colors)]
            xs, ys_hll, errs_hll = _agg(hm, "jaccard_mae_vs_hll")
            valid_hll = [(x, y, e) for x, y, e in zip(xs, ys_hll, errs_hll) if y is not None]
            if valid_hll:
                xs2, ys2, errs2 = zip(*valid_hll)
                axes[2].errorbar(xs2, ys2, yerr=errs2, marker="s", linestyle="--",
                                 color=color, label=f"MAE vs hammock@p_max (subB={sub_b:g})")
            xs, ys_bt, errs_bt = _agg(hm, "jaccard_mae_vs_bt")
            valid_bt = [(x, y, e) for x, y, e in zip(xs, ys_bt, errs_bt) if y is not None]
            if valid_bt:
                xs2, ys2, errs2 = zip(*valid_bt)
                axes[2].errorbar(xs2, ys2, yerr=errs2, marker="o", linestyle="-",
                                 color=color, label=f"MAE vs bedtools (subB={sub_b:g})")
        # Theory curve (HLL 1.04/√(2^p)) — only meaningful for subB=1.0.
        if 1.0 in hm_by_subb:
            xs_th = sorted(hm_by_subb[1.0].keys())
            ys_th = [1.04 / np.sqrt(2 ** x) for x in xs_th]
            axes[2].plot(xs_th, ys_th, "k:", alpha=0.6, label="1.04/√(2^p) (theory)")
        axes[2].set_xlabel(x_label)
        axes[2].set_ylabel("Jaccard MAE")
        axes[2].set_title(f"Accuracy {title_suffix}")
        axes[2].set_yscale("log")
        axes[2].grid(True, alpha=0.3)
        axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_pair_scatter(pair_rows, png_path, title_suffix,
                      value_key="hammock_jaccard",
                      y_label="hammock jaccard (register-equality)"):
    """Scatter of (bedtools_jaccard, hammock estimate), one point per pair, colored by precision.

    With multiple subB values, plots one panel per subB so each panel keeps a
    clean p-colored scatter.

    value_key selects the estimator. The default register-equality column is
    NOT on bedtools' scale -- its offset from y=x is the definitional gap, not
    error -- so the same plot is also emitted for jaccard_similarity_ie, which
    is. Returns False if there is nothing to plot.
    """
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        return

    by_subb_p = defaultdict(lambda: defaultdict(list))
    for r in pair_rows:
        if r.get("bedtools_jaccard") is None or r.get(value_key) is None:
            continue
        sub_b = r.get("sub_b") if r.get("sub_b") is not None else 1.0
        by_subb_p[float(sub_b)][r["precision"]].append((r["bedtools_jaccard"], r[value_key]))

    if not by_subb_p:
        return False

    sub_bs = sorted(by_subb_p.keys(), reverse=True)
    n = len(sub_bs)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7), squeeze=False)
    cmap = plt.get_cmap("viridis")
    for j, sub_b in enumerate(sub_bs):
        ax = axes[0, j]
        ps = sorted(by_subb_p[sub_b].keys())
        for i, p in enumerate(ps):
            xs, ys = zip(*by_subb_p[sub_b][p])
            ax.scatter(xs, ys, s=8, alpha=0.4,
                       color=cmap(i / max(1, len(ps) - 1)),
                       label=f"p={p}")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.5, label="y = x (set-jaccard equality)")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("bedtools jaccard (set Jaccard, ground truth)")
        ax.set_ylabel(y_label)
        ax.set_title(f"subB={sub_b:g}", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=9)
    fig.suptitle(f"{y_label} vs bedtools jaccard\n{title_suffix}", fontsize=10)
    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--axis", choices=["precision", "threads", "intervals"], required=True)
    parser.add_argument("--corpus", choices=["synthetic", "maurano"], default="synthetic",
                        help="BED corpus for the precision axis. 'synthetic' "
                             "regenerates per replicate (the default, and what "
                             "every archived sweep used). 'maurano' uses the 20 "
                             "real fetal DHS files in "
                             "experiments/subB_mixed_stride/data/maurano_sorted/, "
                             "passed as BOTH operands. Use maurano when the run "
                             "has to carry an accuracy axis: the synthetic "
                             "generator draws every file from the same 10 Mbp "
                             "span over 24 chromosomes, so true J is ~0.10 for "
                             "essentially every pair -- a delta function, with "
                             "no range for a MAE-vs-precision curve to be read "
                             "against. Maurano spans J = 0.135-0.627. Ignored "
                             "on the threads and intervals axes, which sweep a "
                             "property of the corpus itself.")
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
    parser.add_argument("--subB-list", dest="sub_b_list", type=str, default="1.0",
                        help="Comma-separated subB values (default '1.0'). "
                             "Each adds a hammock variant invoked with "
                             "--subB <val> --subB-method mixed-stride. subB=1.0 "
                             "emits 'hammock_cpp_B' for backwards compat; other "
                             "values emit 'hammock_cpp_B_subB<val>'.")

    parser.add_argument("--no-ie", dest="ie", action="store_false", default=True,
                        help="Skip the second, untimed --metrics pass on the "
                             "precision axis. That pass is what supplies the "
                             "bedtools-comparable jaccard_similarity_ie "
                             "columns; it re-sketches every input, so it "
                             "roughly DOUBLES the axis's wall time. Recorded "
                             "timings come only from the first pass either way.")

    args = parser.parse_args()
    sub_b_list = [float(x) for x in args.sub_b_list.split(",") if x.strip()]
    if not sub_b_list:
        print("--subB-list must contain at least one value", file=sys.stderr)
        return 1
    for s in sub_b_list:
        if not (0.0 < s <= 1.0):
            print(f"subB values must lie in (0, 1.0]; got {s}", file=sys.stderr)
            return 1
    if args.corpus != "synthetic" and args.axis != "precision":
        # Fail loudly rather than silently ignoring the flag: the threads and
        # intervals axes sweep properties of the generated corpus itself, so
        # there is no coherent fixed-corpus meaning for them.
        print(f"--corpus {args.corpus} is only supported on --axis precision "
              f"(got --axis {args.axis})", file=sys.stderr)
        return 1
    binary = args.binary or find_hammock_cpp()
    if not os.path.exists(binary):
        print(f"hammock-cpp not found at {binary}", file=sys.stderr)
        return 1
    try:
        binary_version = check_binary_version(binary)
        print(f"hammock-cpp version: {binary_version}")
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 1
    if not os.path.exists(BEDTOOLS_SCRIPT):
        print(f"bedtools.sh not found at {BEDTOOLS_SCRIPT}", file=sys.stderr)
        return 1

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.figures_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    sysinfo = get_system_info()
    print(f"system: {sysinfo}")
    print(f"binary: {binary}")
    provenance = {
        "hostname": sysinfo.get("hostname", "unknown"),
        "cpu_model": sysinfo.get("cpu_model", "unknown"),
        "cpu_count": sysinfo.get("cpu_count", "unknown"),
        "slurm_job_id": sysinfo.get("slurm_job_id", "none"),
        "binary_version": binary_version,
        "git_sha": sysinfo.get("git_sha", "unknown"),
    }
    if provenance["slurm_job_id"] == "none":
        print("  NOTE: no SLURM allocation -- cores are shared with anything else "
              "on this node, which inflates wall times unpredictably.")

    pair_rows = []
    subb_tag = ",".join(f"{s:g}" for s in sub_b_list)
    if args.axis == "precision":
        precisions = [int(x) for x in args.precisions.split(",")]
        if args.corpus == "maurano":
            # num_files/num_intervals are inputs to the synthetic generator but
            # only labels on a fixed corpus -- and they reach the CSV and the
            # plot suffix, so they have to describe what actually ran rather
            # than keep whatever default the generator would have used.
            beds = sorted(f for f in os.listdir(MAURANO_DIR) if f.endswith(".bed"))
            nf = len(beds)
            # Coarse interval count from file size (~30 B/line), matching
            # subB_mixed_stride/run_sweep.py's convention. A label, not a
            # measurement.
            avg_bytes = sum(os.path.getsize(os.path.join(MAURANO_DIR, b))
                            for b in beds) / max(1, nf)
            num_intervals = int(avg_bytes / 30)
        else:
            nf = args.num_files if args.num_files is not None else 64
            num_intervals = args.num_intervals
        print(f"[precision] corpus={args.corpus}, p ∈ {precisions}, t={args.threads}, "
              f"files={nf}, intervals≈{num_intervals}, subB ∈ {sub_b_list}")
        rows, pair_rows = sweep_precision(binary, precisions, nf, num_intervals,
                                          args.threads, args.runs, sub_b_list,
                                          ie=args.ie, corpus=args.corpus,
                                          provenance=provenance)
        suffix = (f"({args.corpus}, t={args.threads}, files={nf}, "
                  f"intervals={num_intervals}, subB={subb_tag})")
    elif args.axis == "threads":
        thread_list = [int(x) for x in args.thread_list.split(",")]
        nf = args.num_files if args.num_files is not None else 64
        print(f"[threads] t ∈ {thread_list}, p={args.precision}, files={nf}, "
              f"intervals={args.num_intervals}, subB ∈ {sub_b_list}")
        rows = sweep_threads(binary, thread_list, args.precision, nf,
                             args.num_intervals, args.runs, sub_b_list,
                             provenance=provenance)
        suffix = f"(p={args.precision}, files={nf}, intervals={args.num_intervals}, subB={subb_tag})"
    else:
        intervals_list = [int(x) for x in args.intervals_list.split(",")]
        nf = args.num_files if args.num_files is not None else 16
        print(f"[intervals] N ∈ {intervals_list}, p={args.precision}, t={args.threads}, "
              f"files={nf}, subB ∈ {sub_b_list}")
        rows = sweep_intervals(binary, intervals_list, args.precision, args.threads,
                               nf, args.runs, sub_b_list, provenance=provenance)
        suffix = f"(p={args.precision}, t={args.threads}, files={nf}, subB={subb_tag})"

    stem = f"sweep_{args.axis}_{timestamp}"
    csv_path = os.path.join(args.output_dir, stem + ".csv")
    png_path = os.path.join(args.figures_dir, stem + ".png")
    write_long_csv(rows, csv_path)
    plot_axis(rows, args.axis, png_path, suffix)
    print(f"\nResults: {csv_path}\nFigure:  {png_path}")

    if pair_rows:
        pairs_csv = os.path.join(args.output_dir, stem + "_pairs.csv")
        with open(pairs_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["run_id", "precision", "sub_b",
                                              "query", "reference",
                                              "bedtools_jaccard", "hammock_jaccard",
                                              "hammock_jaccard_ie"],
                               extrasaction="ignore")
            w.writeheader()
            for r in pair_rows:
                w.writerow(r)
        scatter_png = os.path.join(args.figures_dir, stem + "_scatter.png")
        plot_pair_scatter(pair_rows, scatter_png, suffix)
        print(f"Pairs:   {pairs_csv}\nScatter: {scatter_png}")
        ie_png = os.path.join(args.figures_dir, stem + "_scatter_ie.png")
        if plot_pair_scatter(pair_rows, ie_png, suffix,
                             value_key="hammock_jaccard_ie",
                             y_label="hammock jaccard_similarity_ie (inclusion-exclusion)"):
            print(f"Scatter (IE): {ie_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
