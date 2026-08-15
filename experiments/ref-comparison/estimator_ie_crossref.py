#!/usr/bin/env python3
"""Does the cross-reference tissue result hold on `jaccard_similarity_ie`?

Exp A established that peak sets cluster by tissue rather than by the reference
genome they were called against, and picked k=15, w=15 as the lead cell. Every
number behind that was computed on `jaccard_similarity`, the register-equality
column, because that is what the sweep emitted. `jaccard_similarity_ie` is the
bedtools-comparable estimator (CLAUDE.md divergence #2), so this asks the narrow
question: does the ranking survive the column change?

No sketching. `jaccard_similarity_ie` is recovered exactly from the shipped
`containment_AB`/`containment_BA` columns, using the canonical helper.

The statistic is the one `docs/exp_a_results.md` endorses:

    delta = median(same-tissue cross-reference J) - median(different-tissue J)

Deliberately *not* the full-separation criterion `min(xref) > max(diff)`. That
doc already records (2026-08-06) that on `jaccard_similarity` all cells separate
completely, so the criterion is satisfied everywhere and carries no information;
separation is reported here only as a diagnostic column. Rank on delta.

No p-values. The archived Wilcoxon p was computed on n=18/54 *ordered* pairs,
which double-counts every unordered comparison -- the independent units are 9
and 27, over only 9 files. Emitting a fresh p here would re-introduce a number
the experiment already knows is wrong.

Usage:
    python experiments/ref-comparison/estimator_ie_crossref.py
    python experiments/ref-comparison/estimator_ie_crossref.py --csv out.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics as st
import sys

import numpy as np

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "python"))
from hammock.runner import _jaccard_ie_from_containments as _ie  # noqa: E402

RESULTS = os.path.join(_REPO, "experiments", "ref-comparison", "results", "exp_a")
META = os.path.join(_REPO, "docs", "data", "exp_a_metadata.tsv")

# Set once, the first time a row is read with the pre-rename
# `jaccard_similarity` column instead of `reg_eq_similarity`
# (docs/seed-jaccard-reg-eq-rename.md Step 2). `cell_stats` is called once
# per CSV via an outer glob over potentially many "cell" directories, so this
# flag is module-level -- logged once for the whole run, not once per cell.
_REG_EQ_FALLBACK_LOGGED = False


def _reg_eq_value(row: dict, path: str) -> float:
    """Read the register-equality column, preferring the post-rename name.

    Falls back to the legacy `jaccard_similarity` name for archived CSVs
    written before the rename; logs the fallback once for the whole run.
    """
    global _REG_EQ_FALLBACK_LOGGED
    val = row.get("reg_eq_similarity")
    if val is not None:
        return float(val)
    if not _REG_EQ_FALLBACK_LOGGED:
        print("estimator_ie_crossref.py: 'reg_eq_similarity' column not "
              f"found (first seen in {path}); falling back to legacy "
              "'jaccard_similarity' column name (pre-rename hammock "
              "output).", file=sys.stderr)
        _REG_EQ_FALLBACK_LOGGED = True
    return float(row["jaccard_similarity"])


def load_tissue_map(path: str) -> dict[str, str]:
    with open(path) as fh:
        return {r["sample_id"]: r["tissue"] for r in csv.DictReader(fh, delimiter="\t")}


def parse_label(fasta_path: str) -> tuple[str, str]:
    """`.../fastas/<peak_type>/<REF>/<sample_id>.fa` -> (sample_id, ref).

    The reference is the *parent directory*, not part of the filename, so a
    basename-only parse would silently collapse the three references onto one
    label and make every cross-reference pair look like a self-comparison.
    """
    sample = os.path.basename(fasta_path)
    for ext in (".fa", ".fasta", ".fna"):
        if sample.endswith(ext):
            sample = sample[: -len(ext)]
            break
    return sample, os.path.basename(os.path.dirname(fasta_path))


def cell_stats(path: str, tissue_of: dict[str, str]) -> dict | None:
    """Split one cell's pairs into same-tissue-cross-reference vs different-tissue."""
    xref = {"jaccard_similarity": [], "jaccard_similarity_ie": []}
    diff = {"jaccard_similarity": [], "jaccard_similarity_ie": []}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            s1, r1 = parse_label(row["file1"])
            s2, r2 = parse_label(row["file2"])
            if s1 == s2 and r1 == r2:
                continue  # self-comparison
            t1, t2 = tissue_of.get(s1), tissue_of.get(s2)
            if t1 is None or t2 is None:
                return None  # unannotated sample: refuse to guess
            reg = _reg_eq_value(row, path)
            ie = float(_ie(np.array([[float(row["containment_AB"])]]),
                           np.array([[float(row["containment_BA"])]]))[0, 0])
            bucket = xref if t1 == t2 else diff
            bucket["jaccard_similarity"].append(reg)
            bucket["jaccard_similarity_ie"].append(ie)
    if not xref["jaccard_similarity"] or not diff["jaccard_similarity"]:
        return None
    out = {"n_xref": len(xref["jaccard_similarity"]),
           "n_diff": len(diff["jaccard_similarity"])}
    for col in ("jaccard_similarity", "jaccard_similarity_ie"):
        tag = "reg" if col == "jaccard_similarity" else "ie"
        out[f"median_xref_{tag}"] = st.median(xref[col])
        out[f"median_diff_{tag}"] = st.median(diff[col])
        out[f"delta_{tag}"] = st.median(xref[col]) - st.median(diff[col])
        out[f"separated_{tag}"] = min(xref[col]) > max(diff[col])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=RESULTS)
    ap.add_argument("--meta", default=META)
    ap.add_argument("--csv", help="write the per-cell table here")
    args = ap.parse_args()

    tissue_of = load_tissue_map(args.meta)
    rows = []
    for peak_type in sorted(
            os.path.basename(d) for d in glob.glob(os.path.join(args.results, "*"))
            if os.path.isdir(d)):
        for cell_dir in sorted(glob.glob(os.path.join(args.results, peak_type, "k*_w*"))):
            # cell_stats needs containment_AB/containment_BA, present only in
            # the "_full" (--metrics) shape -- narrow the glob so a "_ie"/"_re"
            # file from the same directory can never be picked instead (was
            # `sorted(glob.glob(*.csv))[0]`, order-dependent across shapes).
            csvs = sorted(glob.glob(os.path.join(cell_dir, "*_full.csv")))
            if not csvs:
                continue
            stats = cell_stats(csvs[0], tissue_of)
            if stats is None:
                continue
            k, w = os.path.basename(cell_dir).split("_")
            rows.append({"peak_type": peak_type, "cell": os.path.basename(cell_dir),
                         "k": int(k[1:]), "w": int(w[1:]), **stats})

    if not rows:
        print(f"no usable cells found under {args.results}", file=sys.stderr)
        return 1

    for peak_type in sorted({r["peak_type"] for r in rows}):
        sub = sorted((r for r in rows if r["peak_type"] == peak_type),
                     key=lambda r: -r["delta_reg"])
        print(f"\n=== {peak_type} peaks: cells ranked by delta on "
              f"jaccard_similarity (n_xref={sub[0]['n_xref']}, "
              f"n_diff={sub[0]['n_diff']}; ordered pairs) ===")
        print(f"{'cell':>10} {'delta_reg':>10} {'delta_ie':>10} "
              f"{'rank_reg':>9} {'rank_ie':>8} {'sep_reg':>8} {'sep_ie':>7}")
        by_ie = {r["cell"]: i + 1 for i, r in
                 enumerate(sorted(sub, key=lambda r: -r["delta_ie"]))}
        for i, r in enumerate(sub):
            print(f"{r['cell']:>10} {r['delta_reg']:>10.4f} {r['delta_ie']:>10.4f} "
                  f"{i + 1:>9} {by_ie[r['cell']]:>8} "
                  f"{str(r['separated_reg']):>8} {str(r['separated_ie']):>7}")
        # Report the ranking's stability, not just whether the argmax moved.
        # The top cells sit within a few thousandths of each other, so which one
        # is nominally first is not a meaningful difference; the rank
        # correlation and the lead margin are what say whether the ordering is
        # the same ordering.
        n = len(sub)
        d2 = sum((i + 1 - by_ie[r["cell"]]) ** 2 for i, r in enumerate(sub))
        rho = 1 - 6 * d2 / (n * (n * n - 1))
        lead_reg = sub[0]["cell"]
        lead_ie = min(sub, key=lambda r: by_ie[r["cell"]])["cell"]
        margin_reg = sub[0]["delta_reg"] - sub[1]["delta_reg"]
        ie_sorted = sorted(sub, key=lambda r: -r["delta_ie"])
        margin_ie = ie_sorted[0]["delta_ie"] - ie_sorted[1]["delta_ie"]
        print(f"  Spearman rho between the two rankings: {rho:.4f} "
              f"({sum(1 for i, r in enumerate(sub) if i + 1 != by_ie[r['cell']])}"
              f"/{n} cells change rank)")
        print(f"  top cell: {lead_reg} (register-equality, margin "
              f"{margin_reg:+.4f}) vs {lead_ie} (IE, margin {margin_ie:+.4f})")
        # The claim Exp A actually rests on is the k>=15 / k<=10 split, not the
        # identity of the single best cell.
        hi = [r for r in sub if r["k"] >= 15]
        lo = [r for r in sub if r["k"] <= 10]
        if hi and lo:
            for tag in ("reg", "ie"):
                gap = min(r[f"delta_{tag}"] for r in hi) - max(r[f"delta_{tag}"] for r in lo)
                print(f"  k>=15 vs k<=10 delta gap ({tag}): {gap:+.4f} "
                      f"({'disjoint' if gap > 0 else 'OVERLAPPING'})")
        else:
            print("  k>=15 vs k<=10 delta gap: not reported "
                  "(input does not contain both groups)")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
        print(f"\nWrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
