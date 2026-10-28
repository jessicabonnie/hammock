#!/usr/bin/env python3
"""Side-by-side summaries: archived 0.7.0 (or 0.5.0) numbers vs the 0.12.0 single-build rerun.

Reads the archived CSVs read-only and the new CSVs under results/, and writes
only into summary/finalize_<tag>/ (a fresh directory; refuses to reuse one).
Any benchmark whose new output is missing is reported and skipped, so a partial
rerun still yields summaries for what finished.
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

EXP = os.path.dirname(os.path.abspath(__file__))
REPO = "/home/jbonnie1/interval_sketch/hammock_claude"
DOCS = os.path.join(REPO, "docs", "data")
OLD_RESULTS = "/vast/blangme2/jbonnie/hammock_claude_experiments/bedtools_benchmark/results"

OLD = {
    "fig3a": os.path.join(DOCS, "cpp_vs_bedtools_t16_p18.csv"),  # job 29671317, 0.7.0
    "threads": os.path.join(DOCS, "sweep_threads_p18.csv"),  # job 29670792, 0.7.0
    "precision": os.path.join(DOCS, "sweep_precision_maurano_p18_t16.csv"),  # job 29670793, 0.7.0
    "precision_pairs": os.path.join(OLD_RESULTS, "sweep_precision_20260809_230812_pairs.csv"),
    "acc18": os.path.join(DOCS, "hammock_hll_p18_jaccB_full.csv"),  # hammock 0.5.0 (Python CLI)
    "acc21": os.path.join(DOCS, "hammock_hll_p21_jaccB_full.csv"),  # hammock 0.5.0 (Python CLI)
    "truth": os.path.join(DOCS, "maurano_bedtools_ref.tsv"),
    "pcost": os.path.join(OLD_RESULTS, "pairwise_cost_by_precision_20260808_183458.csv"),  # job 29620391, 0.7.0
}

notes = []


def note(msg):
    print(msg)
    notes.append(msg)


def newest(pattern, exclude=()):
    files = [f for f in sorted(glob.glob(pattern)) if not any(e in f for e in exclude)]
    return files


def fig3a(out):
    files = newest(os.path.join(EXP, "results", "fig3a", "cpp_vs_bedtools_t16_*.csv"),
                   exclude=("_runs.csv",))
    if not files:
        note("fig3a: no new CSV; skipped")
        return
    new = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    if new.duplicated(["num_files", "tool"]).any():
        note(f"fig3a: duplicate (N, tool) rows across {files}; keeping last")
        new = new.drop_duplicates(["num_files", "tool"], keep="last")
    old = pd.read_csv(OLD["fig3a"])
    cols = ["mean_wall_time", "std_wall_time", "mean_sketch_creation_time",
            "mean_comparison_time", "mean_pair_time", "mean_bedtools_parallel_eff"]
    m = old[["num_files", "tool"] + cols].merge(
        new[["num_files", "tool"] + cols + ["hostname", "slurm_job_id", "binary_version"]],
        on=["num_files", "tool"], how="outer", suffixes=("_070", "_0120"))
    m["wall_ratio_0120_over_070"] = m["mean_wall_time_0120"] / m["mean_wall_time_070"]
    for v in ("070", "0120"):
        m[f"sketch_frac_{v}"] = m[f"mean_sketch_creation_time_{v}"] / m[f"mean_wall_time_{v}"]
        m[f"comparison_frac_{v}"] = m[f"mean_comparison_time_{v}"] / m[f"mean_wall_time_{v}"]
    m = m.sort_values(["num_files", "tool"])
    m.to_csv(os.path.join(out, "fig3a_side_by_side.csv"), index=False)

    # Speedup table: bedtools / hammock IE arm (the plotted headline curve)
    rows = []
    for n in sorted(set(m["num_files"])):
        r = {"num_files": n}
        for v in ("070", "0120"):
            sub = m[m["num_files"] == n].set_index("tool")[f"mean_wall_time_{v}"]
            bt, ie, re_ = sub.get("bedtools"), sub.get("hammock_ie_B"), sub.get("hammock_cpp_B")
            r[f"bedtools_over_hammock_ie_{v}"] = bt / ie if bt and ie else np.nan
            r[f"bedtools_over_hammock_re_{v}"] = bt / re_ if bt and re_ else np.nan
            r[f"re_wall_reduction_vs_ie_{v}"] = 1 - re_ / ie if re_ and ie else np.nan
        rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(out, "fig3a_speedup.csv"), index=False)
    note(f"fig3a: new inputs {files}")


def threads(out):
    files = newest(os.path.join(EXP, "results", "threads", "sweep_threads_*.csv"))
    if not files:
        note("threads: no new CSV; skipped")
        return
    new = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    old = pd.read_csv(OLD["threads"])
    res = []
    for label, df in (("070", old), ("0120", new)):
        g = df.groupby(["tool", "threads"])["wall_time"].agg(["mean", "median", "std", "count"]).reset_index()
        g.columns = ["tool", "threads"] + [f"{c}_wall_{label}" for c in ("mean", "median", "std", "n")]
        t1 = g[g["threads"] == 1].set_index("tool")[f"mean_wall_{label}"]
        g[f"speedup_vs_t1_{label}"] = g.apply(lambda r: t1[r["tool"]] / r[f"mean_wall_{label}"], axis=1)
        res.append(g)
    m = res[0].merge(res[1], on=["tool", "threads"], how="outer")
    m["wall_ratio_0120_over_070"] = m["mean_wall_0120"] / m["mean_wall_070"]
    m.sort_values(["tool", "threads"]).to_csv(os.path.join(out, "threads_side_by_side.csv"), index=False)
    note(f"threads: new inputs {files}")


def _unordered(df, a, b):
    """First occurrence of each unordered pair, self pairs dropped (R distinct())."""
    key = df.apply(lambda r: tuple(sorted((r[a], r[b]))), axis=1)
    df = df.assign(_key=key)
    df = df[df[a] != df[b]]
    return df.drop_duplicates("_key", keep="first")


def _agree(est, truth):
    est, truth = np.asarray(est, float), np.asarray(truth, float)
    gap = est - truth
    return {
        "n": len(est),
        "mae": np.mean(np.abs(gap)),
        "rmse": np.sqrt(np.mean(gap ** 2)),
        "pearson": stats.pearsonr(est, truth)[0],
        "spearman": stats.spearmanr(est, truth)[0],
        "kendall_tau_b": stats.kendalltau(est, truth)[0],
    }


def precision(out):
    files = newest(os.path.join(EXP, "results", "precision", "sweep_precision_*.csv"),
                   exclude=("_pairs.csv",))
    if not files:
        note("precision: no new CSV; skipped")
        return
    new = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    old = pd.read_csv(OLD["precision"])
    res = []
    for label, df in (("070", old), ("0120", new)):
        bt_med = df[df["tool"] == "bedtools"]["wall_time"].median()
        hm = df[df["tool"] == "hammock_cpp_B"]
        g = hm.groupby("precision").agg(
            median_wall=("wall_time", "median"), min_wall=("wall_time", "min"),
            max_wall=("wall_time", "max"), median_ie_pass_wall=("metrics_wall_time", "median"),
            median_sketch=("sketch_creation_time", "median"),
            median_pair=("pair_time", "median"),
            ie_mae_vs_bt=("jaccard_ie_mae_vs_bt", "mean"),
            ie_mae_vs_bt_sd=("jaccard_ie_mae_vs_bt", "std"),
            re_mae_vs_bt=("jaccard_mae_vs_bt", "mean"),
            n=("wall_time", "count")).reset_index()
        g["median_bedtools_wall"] = bt_med
        g["rel_wall_hammock_over_bedtools"] = g["median_wall"] / bt_med
        g["ie_pass_vs_re_pct"] = 100 * (g["median_ie_pass_wall"] / g["median_wall"] - 1)
        g = g.rename(columns={c: f"{c}_{label}" for c in g.columns if c != "precision"})
        res.append(g)
    m = res[0].merge(res[1], on="precision", how="outer")
    m["rel_wall_ratio_0120_over_070"] = (m["rel_wall_hammock_over_bedtools_0120"]
                                         / m["rel_wall_hammock_over_bedtools_070"])

    # Rank/linear agreement from the per-pair dumps (replicate 0; estimates are
    # deterministic, bedtools is exact).
    pfiles = newest(os.path.join(EXP, "results", "precision", "sweep_precision_*_pairs.csv"))
    for label, path in (("070", OLD["precision_pairs"]), ("0120", pfiles[-1] if pfiles else None)):
        if not path or not os.path.exists(path):
            note(f"precision: pairs file for {label} missing")
            continue
        pr = pd.read_csv(path)
        pr = pr[pr["run_id"] == 0]
        recs = []
        for p, d in pr.groupby("precision"):
            d = _unordered(d, "query", "reference")
            a = _agree(d["hammock_jaccard_ie"], d["bedtools_jaccard"])
            recs.append({"precision": p, f"ie_pearson_{label}": a["pearson"],
                         f"ie_kendall_{label}": a["kendall_tau_b"],
                         f"ie_mae_unordered_{label}": a["mae"], f"n_unordered_{label}": a["n"]})
        m = m.merge(pd.DataFrame(recs), on="precision", how="left")
    m.to_csv(os.path.join(out, "precision_side_by_side.csv"), index=False)
    note(f"precision: new inputs {files} + {pfiles}")


def accuracy(out):
    truth = pd.read_csv(OLD["truth"], sep="\t")
    truth = _unordered(truth, "file1", "file2")[["_key", "jaccard"]]
    recs, pair_tables = [], []
    for p in (18, 21):
        newf = newest(os.path.join(EXP, "results", "accuracy", f"maurano_job*_hll_p{p}_jaccB_full.csv"))
        if not newf:
            note(f"accuracy p={p}: no new CSV; skipped")
            continue
        old = pd.read_csv(OLD[f"acc{p}"])
        new = pd.read_csv(newf[-1], sep="\t").rename(columns={"query": "file1", "reference": "file2"})
        o = _unordered(old, "file1", "file2")
        n = _unordered(new, "file1", "file2")
        j = (o[["_key", "file1", "file2", "reg_eq_similarity", "jaccard_similarity_ie"]]
             .merge(n[["_key", "reg_eq_similarity", "jaccard_similarity_ie"]], on="_key",
                    suffixes=("_050", "_0120"))
             .merge(truth, on="_key"))
        if len(j) != 190:
            note(f"accuracy p={p}: joined {len(j)} unordered pairs, expected 190")
        for est in ("jaccard_similarity_ie", "reg_eq_similarity"):
            for v in ("050", "0120"):
                a = _agree(j[f"{est}_{v}"], j["jaccard"])
                recs.append({"precision": p, "estimator": est, "build": v, **a})
            d = (j[f"{est}_0120"] - j[f"{est}_050"]).abs()
            recs.append({"precision": p, "estimator": est, "build": "0120_minus_050",
                         "n": len(d), "max_abs_diff": d.max(),
                         "n_identical": int((d == 0).sum()),
                         "n_within_1e-12": int((d <= 1e-12).sum())})
        j.insert(0, "precision", p)
        pair_tables.append(j.drop(columns="_key").rename(columns={"jaccard": "bedtools_jaccard"}))
        note(f"accuracy p={p}: new input {newf[-1]}")
    if recs:
        pd.DataFrame(recs).to_csv(os.path.join(out, "accuracy_side_by_side.csv"), index=False)
        pd.concat(pair_tables).to_csv(os.path.join(out, "accuracy_pairs_050_vs_0120.csv"), index=False)


def pairwise_cost(out):
    srcs = {"070_6448y": OLD["pcost"]}
    for tag in ("6448y", "6248r"):
        f = newest(os.path.join(EXP, "results", f"pairwise_cost_{tag}", "pairwise_cost_by_precision_*.csv"))
        if f:
            srcs[f"0120_{tag}"] = f[-1]
            note(f"pairwise_cost {tag}: new input {f[-1]}")
        else:
            note(f"pairwise_cost {tag}: no new CSV; skipped")
    if len(srcs) == 1:
        return
    tabs = []
    for label, path in srcs.items():
        df = pd.read_csv(path)
        g = df.groupby(["precision", "arm"])["us_per_pair"].median().unstack("arm")
        g = g.rename(columns={"no_metrics": f"re_us_per_pair_{label}",
                              "metrics": f"ie_full_us_per_pair_{label}"})
        g[f"ie_over_re_{label}"] = g[f"ie_full_us_per_pair_{label}"] / g[f"re_us_per_pair_{label}"]
        tabs.append(g)
    m = pd.concat(tabs, axis=1).reset_index()
    m.to_csv(os.path.join(out, "pairwise_cost_side_by_side.csv"), index=False)


def main():
    tag = os.environ.get("SLURM_JOB_ID") or pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    out = os.path.join(EXP, "summary", f"finalize_{tag}")
    os.makedirs(out, exist_ok=False)
    for fn in (fig3a, threads, precision, accuracy, pairwise_cost):
        try:
            fn(out)
        except Exception as e:  # keep going; one failed benchmark must not hide the rest
            note(f"{fn.__name__}: FAILED with {type(e).__name__}: {e}")
    with open(os.path.join(out, "NOTES.txt"), "w") as f:
        f.write("\n".join(notes) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
