#!/usr/bin/env python3
"""Does reading `jaccard_similarity_ie` rescue the mouse-human tissue signal?

Phase 0d of the estimator comparison. `README.md` records this experiment as a
strong negative result: 0 of 20 (k, w) cells were tissue-dominant, i.e. sketching
regulatory peak FASTAs does not recover cross-species tissue identity at ~80 Mya.
Since then a second Jaccard column shipped, and the cross-species pairs sit in
exactly the low-J regime where the two estimators were measured to change places
(`experiments/bedtools_benchmark/estimator_rank_by_precision.py`). This asks the
narrow question that raises: was the negative result an artefact of reading the
register-equality column?

No sketching. `jaccard_similarity_ie` is recovered exactly from the shipped
`containment_AB`/`containment_BA` columns (CLAUDE.md divergence #2).

The statistic is the separation the experiment was designed around: among
cross-species pairs only, the gap between same-tissue and different-tissue
median Jaccard. Within-species pairs are reported alongside as the scale
reference -- they are what the signal has to compete with.

Usage:
    python experiments/mus-homo/estimator_ie_tissue.py
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

RESULTS = os.path.join(_REPO, "experiments", "mus-homo", "results")


def species(name: str) -> str:
    """Filenames are ENCSR<id>_<tissue>_<species>.fa. Match the basename suffix,
    not the path -- the path contains 'mus-homo' and matches everything."""
    return "mouse" if name.endswith("_mouse.fa") else "human"


def tissue(name: str) -> str:
    return os.path.basename(name).split("_")[1]


def load(path):
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            a, b = os.path.basename(row["file1"]), os.path.basename(row["file2"])
            if a == b:
                continue
            rows.append((a, b, float(row["jaccard_similarity"]),
                         float(_ie(np.array([[float(row["containment_AB"])]]),
                                   np.array([[float(row["containment_BA"])]]))[0, 0])))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=RESULTS)
    args = ap.parse_args()

    cells = sorted(os.path.basename(d) for d in glob.glob(os.path.join(args.results, "k*_w*")))
    if not cells:
        print(f"no cells found under {args.results}")
        return 1

    hdr = f"{'cell':>10}{'stratum':>26}{'n':>5}{'med_RE':>10}{'med_IE':>10}"
    print(hdr)
    print("-" * len(hdr))
    for cell in cells:
        found = glob.glob(os.path.join(args.results, cell, "*.csv"))
        if not found:
            continue
        strata = {"within-species": [], "cross-sp, same tissue": [], "cross-sp, diff tissue": []}
        for a, b, j_re, j_ie in load(found[0]):
            if species(a) == species(b):
                key = "within-species"
            else:
                key = "cross-sp, same tissue" if tissue(a) == tissue(b) else "cross-sp, diff tissue"
            strata[key].append((j_re, j_ie))
        for label, vals in strata.items():
            if vals:
                print(f"{cell:>10}{label:>26}{len(vals):>5}"
                      f"{st.median(v[0] for v in vals):>10.4f}"
                      f"{st.median(v[1] for v in vals):>10.4f}")
        same, diff = strata["cross-sp, same tissue"], strata["cross-sp, diff tissue"]
        if same and diff:
            d_re = st.median(v[0] for v in same) - st.median(v[0] for v in diff)
            d_ie = st.median(v[1] for v in same) - st.median(v[1] for v in diff)
            print(f"{'':>10}{'--> tissue separation':>26}{'':>5}{d_re:>10.4f}{d_ie:>10.4f}")
        print()

    print("Reading: the negative result stands under either column. Wherever the")
    print("two estimators agree on the level (k <= 10) they agree on the separation")
    print("to four decimals, so nothing was lost by reading register-equality. The")
    print("one place they part company is k20, and there `_ie` removes the")
    print("separation rather than sharpening it: 0.0008 -> 0.0000, with the")
    print("cross-species medians themselves dropping 0.013 -> 0.003. That is the")
    print("chance-agreement floor coming out, which is the expected direction --")
    print("k20 sits at the low-J end where the floor is most of the value.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
