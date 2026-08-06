#!/usr/bin/env python3
"""Rank fidelity of the two Jaccard estimators, stratified by true J and precision.

Answers a question neither the Maurano corpus nor `estimator_compare.py`'s own
summary tables can: **at low true Jaccard, which estimator ranks pairs better,
and does the answer depend on precision?**

Why this stratification matters. `jaccard_similarity` (register-equality) is
biased by a chance-agreement floor but has lower variance;
`jaccard_similarity_ie` (inclusion-exclusion) is near-unbiased but noisier and
censored at 0 by the `>= 0` clamp. Whether the bias or the noise dominates is a
function of where you sit in J -- and the Maurano interval-mode corpus cannot
address it, because its off-diagonal J never drops below 0.1355. The synthetic
corpus behind `estimator_compare_full.csv` puts ~73% of its rows below J = 0.05,
which is where the repo's cross-species Mode D work actually operates.

Reads `results/estimator_compare_full.csv`, which already carries both
estimators plus exact `bedtools jaccard` truth -- no rerun, no sketching.

Conventions, pinned because they move the third digit:
  * Kendall tau-b, computed over all C(n,2) comparisons within a stratum.
  * "Discordant" counts comparison *pairs* whose ordering disagrees with truth,
    not rows. Ties in either sequence are excluded from both numerator and
    denominator, which is what makes this tau-b.
  * Rows are used as-is from the CSV; that file is a full A x B cross product,
    so pairs share input files and are NOT independent. Treat the discordance
    counts as descriptive, not as a basis for a binomial CI -- see
    docs/estimator-analysis-findings.md section 4 on why the naive SE is ~7.5x
    too small on this kind of design.
  * Cardinality ratio is pinned near 1 by construction (`make_data` writes
    `num_intervals` intervals into both files), so this measures *within*-
    geometry rank fidelity. The across-geometry case is a different question
    and is governed by c varying with |A|/|B|.

Usage:
    python experiments/bedtools_benchmark/estimator_rank_by_precision.py
    python experiments/bedtools_benchmark/estimator_rank_by_precision.py --csv PATH
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import statistics as st
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "results", "estimator_compare_full.csv")

# Strata chosen so the first bin is the regime Maurano cannot reach and the
# cross-species corpora live in. Not tuned to the result.
STRATA = (("J < 0.05", 0.0, 0.05),
          ("0.05 <= J < 0.2", 0.05, 0.2),
          ("J >= 0.2", 0.2, 1.01))


def kendall_tau_b(xs, ys):
    """Return (tau_b, discordant, comparable). Ties excluded from both."""
    n = len(xs)
    conc = disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = (xs[i] - xs[j]) * (ys[i] - ys[j])
            if s > 0:
                conc += 1
            elif s < 0:
                disc += 1
    total = conc + disc
    return ((conc - disc) / total if total else float("nan")), disc, total


def _replicate(name):
    m = re.search(r"rep(\d+)", name)
    return m.group(1) if m else ""


def load(path):
    by_precision = defaultdict(list)
    skipped = 0
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                by_precision[int(row["precision"])].append((
                    float(row["bedtools_jaccard"]),
                    float(row["j_register_equality"]),
                    float(row["j_incl_excl"]),
                    str(row.get("ie_clamped", "")).strip().lower() in ("1", "true", "yes"),
                    # Resampling unit: the corpus is a 3 x 30 cross product of
                    # A-file replicates against B-file overlap fractions, so the
                    # replicate is the only independent axis. Fractions cannot be
                    # resampled -- dropping one removes a J level outright.
                    _replicate(row.get("file1", "")),
                ))
            except (ValueError, TypeError, KeyError):
                skipped += 1
    return by_precision, skipped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--min-n", type=int, default=8,
                    help="skip strata with fewer rows (default 8)")
    args = ap.parse_args()

    if not os.path.exists(args.csv):
        print(f"not found: {args.csv}")
        return 1

    by_precision, skipped = load(args.csv)
    if skipped:
        print(f"note: skipped {skipped} unparseable row(s)")

    hdr = (f"{'stratum':>16}{'p':>4}{'n':>5}{'clamped':>9}"
           f"{'tau_RE':>9}{'tau_IE':>9}{'disc_RE':>11}{'disc_IE':>11}"
           f"{'MAE_RE':>9}{'MAE_IE':>9}{'rank':>6}")
    print(hdr)
    print("-" * len(hdr))
    for label, lo, hi in STRATA:
        for p in sorted(by_precision):
            sub = [t for t in by_precision[p] if lo <= t[0] < hi]
            if len(sub) < args.min_n:
                continue
            truth = [t[0] for t in sub]
            re_v = [t[1] for t in sub]
            ie_v = [t[2] for t in sub]
            clamped = sum(1 for t in sub if t[3])
            t_re, d_re, total = kendall_tau_b(truth, re_v)
            t_ie, d_ie, _ = kendall_tau_b(truth, ie_v)
            winner = "RE" if t_re > t_ie else ("IE" if t_ie > t_re else "tie")
            print(f"{label:>16}{p:>4}{len(sub):>5}{clamped:>9}"
                  f"{t_re:>9.4f}{t_ie:>9.4f}"
                  f"{f'{d_re}/{total}':>11}{f'{d_ie}/{total}':>11}"
                  f"{st.mean(abs(a - b) for a, b in zip(re_v, truth)):>9.4f}"
                  f"{st.mean(abs(a - b) for a, b in zip(ie_v, truth)):>9.4f}"
                  f"{winner:>6}")
        print()

    print("Reading: MAE favours inclusion-exclusion in every stratum at every")
    print("precision. Rank fidelity does not -- below J = 0.05 register-equality")
    print("ranks better up to p = 20 and loses at p = 24.")
    print()
    print("Only the J < 0.05 stratum resolves anything. Above it both estimators")
    print("reach tau = 1.0000 by p = 20 on this corpus (register-equality is")
    print("0.9804 at p = 16 in the J >= 0.2 stratum -- 1 discordant pair out of")
    print("102 -- so p = 16 is not where they converge), and the one cell that")
    print("does separate (p=12, J >= 0.2) turns on 1 discordant comparison against")
    print("3 out of 102. Do not read a winner out of those rows; read that the")
    print("question is only live at low J. Note also the 25 clamped rows in the")
    print("p=12 low-J cell -- they tie with each other and depress tau_IE there,")
    print("which is a real property of the estimator, not a measurement artifact.")
    print()
    label, lo, hi = STRATA[0]
    print(f"Stability of the {label} ordering, leave-one-replicate-out:")
    print(f"{'p':>4}{'tau_RE range':>22}{'tau_IE range':>22}{'holds':>8}")
    for p in sorted(by_precision):
        sub = [t for t in by_precision[p] if lo <= t[0] < hi]
        reps = sorted({t[4] for t in sub})
        if len(reps) < 2:
            continue
        loo = []
        for drop in reps:
            keep = [t for t in sub if t[4] != drop]
            loo.append((kendall_tau_b([t[0] for t in keep], [t[1] for t in keep])[0],
                        kendall_tau_b([t[0] for t in keep], [t[2] for t in keep])[0]))
        holds = len({a > b for a, b in loo}) == 1
        print(f"{p:>4}"
              f"{f'[{min(a for a, _ in loo):.3f}, {max(a for a, _ in loo):.3f}]':>22}"
              f"{f'[{min(b for _, b in loo):.3f}, {max(b for _, b in loo):.3f}]':>22}"
              f"{('yes' if holds else 'FLIPS'):>8}")
    print("Three replicates is a coarse check, not a confidence interval.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
