#!/usr/bin/env python3
"""Does --subB degrade `jaccard_similarity_ie`? Analysis for the Maurano grid.

Primary statistic is accuracy: |IE - bedtools| per pair, since IE is the
bedtools-calibrated column. Register-equality is reported only as *drift* vs its
own subB=1.0 baseline -- it is not on bedtools' scale (chance-agreement floor,
CLAUDE.md divergence #2) and its level must never be compared to IE's.

Pairs are deduplicated to the 190 unordered pairs. The archived analysis counted
190 pairs x 5 byte-identical replicates as 950 observations; that is
pseudo-replication and is not repeated here.

Uncertainty is leave-one-file-out jackknife over the 20 files, which is what
docs/estimator-analysis-findings.md section 4 endorses at this n -- the 190 pairs
are not independent (each file appears in 19 of them).
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

REPO = "/home/jbonnie1/interval_sketch/hammock_claude"
REF = os.path.join(REPO, "docs/data/maurano_bedtools_ref.tsv")
GRID = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else None


def load_truth():
    j, cov = {}, {}
    with open(REF) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["file1"] == r["file2"]:
                cov[r["file1"]] = float(r["union"])   # self-union == covered bp
            else:
                j[frozenset((r["file1"], r["file2"]))] = float(r["jaccard"])
    return j, cov


def jackknife(pairs, files, stat):
    """pairs: list of (a, b, value). stat: callable over a list of values."""
    reps = []
    for f in files:
        sub = [v for a, b, v in pairs if a != f and b != f]
        if sub:
            reps.append(stat(sub))
    n = len(reps)
    if n < 2:
        return float("nan")
    m = sum(reps) / n
    return math.sqrt((n - 1) / n * sum((x - m) ** 2 for x in reps))


def main() -> int:
    truth, cov = load_truth()
    median_cov = sorted(cov.values())[len(cov) // 2]

    # cell -> unordered pair -> row
    cells = defaultdict(dict)
    with open(GRID) as fh:
        for r in csv.DictReader(fh):
            key = (int(r["precision"]), float(r["subB"]), r["method"])
            cells[key].setdefault(frozenset((r["file_a"], r["file_b"])), r)

    # --- containment self-check: clamp FIRST, then derive; 2 ulp tolerance -----
    # hammock_cli.cpp clamps containments to 1.0 before deriving _ie but emits
    # them unclamped, so a strict-equality check false-alarms by construction.
    worst, n_checked, n_clamped = 0.0, 0, 0
    for rows in cells.values():
        for r in rows.values():
            ca, cb = float(r["containment_AB"]), float(r["containment_BA"])
            if ca > 1.0 or cb > 1.0:
                n_clamped += 1
            ca, cb = min(ca, 1.0), min(cb, 1.0)
            want = 0.0 if (ca <= 0 or cb <= 0) else 1.0 / (1.0 / ca + 1.0 / cb - 1.0)
            got = float(r["jaccard_similarity_ie"])
            denom = max(abs(want), 1e-300)
            worst = max(worst, abs(want - got) / denom)
            n_checked += 1
    print(f"containment self-check: {n_checked} rows, max relative deviation "
          f"{worst:.2e} ({n_clamped} rows had a containment > 1.0)")
    print("  (expected ~1e-16; a large value means the parse or the clamp is wrong)\n")

    files = sorted(cov)
    out_rows = []
    for (p, sub_b, method) in sorted(cells):
        rows = cells[(p, sub_b, method)]
        ie = [(a, b, float(r["jaccard_similarity_ie"]) - truth[k])
              for k, r in rows.items()
              if k in truth for a, b in [tuple(k)]]
        n = len(ie)
        mae = sum(abs(v) for _, _, v in ie) / n
        bias = sum(v for _, _, v in ie) / n
        se_mae = jackknife([(a, b, abs(v)) for a, b, v in ie], files,
                           lambda xs: sum(xs) / len(xs))
        se_bias = jackknife(ie, files, lambda xs: sum(xs) / len(xs))
        lam = median_cov * sub_b / (2 ** p)
        n_zero = sum(1 for r in rows.values()
                     if float(r["jaccard_similarity_ie"]) == 0.0)
        out_rows.append({
            "precision": p, "subB": sub_b, "method": method, "n_pairs": n,
            "lambda_median_file": lam,
            "ie_mae_vs_bedtools": mae, "ie_mae_se": se_mae,
            "ie_signed_bias": bias, "ie_bias_se": se_bias,
            "frac_ie_zero": n_zero / len(rows),
        })

    # --- register-equality drift vs its OWN subB=1.0 baseline (secondary) -----
    for p in sorted({r["precision"] for r in out_rows}):
        base = cells.get((p, 1.0, "none"))
        if not base:
            continue
        bl = {k: float(v["jaccard_similarity"]) for k, v in base.items()}
        for row in out_rows:
            if row["precision"] != p:
                continue
            rows = cells[(p, row["subB"], row["method"])]
            d = [abs(float(v["jaccard_similarity"]) - bl[k])
                 for k, v in rows.items() if k in bl]
            row["re_drift_vs_nosub"] = sum(d) / len(d) if d else float("nan")

    hdr = ["precision", "subB", "method", "n_pairs", "lambda_median_file",
           "ie_mae_vs_bedtools", "ie_mae_se", "ie_signed_bias", "ie_bias_se",
           "frac_ie_zero", "re_drift_vs_nosub"]
    for p in sorted({r["precision"] for r in out_rows}):
        floor = next(r for r in out_rows if r["precision"] == p and r["subB"] == 1.0)
        print(f"=== p={p}  (subB=1.0 floor: MAE={floor['ie_mae_vs_bedtools']:.3e} "
              f"+/- {floor['ie_mae_se']:.1e}) ===")
        print(f"  {'subB':>5} {'method':<15} {'lambda':>8} {'IE MAE':>10} {'SE':>9} "
              f"{'x floor':>8} {'IE bias':>10} {'RE drift':>9}")
        for r in sorted((r for r in out_rows if r["precision"] == p),
                        key=lambda r: (-r["subB"], r["method"])):
            print(f"  {r['subB']:>5g} {r['method']:<15} {r['lambda_median_file']:>8.2f} "
                  f"{r['ie_mae_vs_bedtools']:>10.3e} {r['ie_mae_se']:>9.1e} "
                  f"{r['ie_mae_vs_bedtools']/floor['ie_mae_vs_bedtools']:>8.2f} "
                  f"{r['ie_signed_bias']:>10.2e} {r.get('re_drift_vs_nosub',0):>9.2e}")
        print()

    if OUT:
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=hdr)
            w.writeheader()
            w.writerows(out_rows)
        print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
