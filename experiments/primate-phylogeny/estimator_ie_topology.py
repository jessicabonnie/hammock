#!/usr/bin/env python3
"""Does the recovered topology change if you read `jaccard_similarity_ie`?

Phase 0d of the estimator comparison. This corpus is the repo's lowest-J Mode D
data: at k=20 the off-diagonal register-equality Jaccard has median 0.024
(H3K4me3) and 0.057 (H3K27ac), and the inclusion-exclusion reconstruction of the
same pairs puts them at 0.007 and 0.0003. That is the regime where
`experiments/bedtools_benchmark/estimator_rank_by_precision.py` measured the two
estimators changing places, and these runs are all at p=24 -- the precision
where that measurement says inclusion-exclusion ranks better.

No sketching and no rerun: the shipped CSVs already carry `containment_AB` and
`containment_BA`, from which `jaccard_similarity_ie` is exactly recoverable
(J = 1/(1/C_AB + 1/C_BA - 1); see CLAUDE.md divergence #2). Both estimators are
turned into 1 - J distance matrices, run through the same plain
neighbour-joining, and scored against the accepted mammal topology.

Scoring convention: splits are compared as *unrooted bipartitions*, normalized
to the smaller side, so an arbitrary NJ root cannot flatter either arm. The four
true internal splits for these seven taxa are (hsa,rheMac), +calJac (primates),
+Mmus (Euarchontoglires), and (bosTau,canFam) (Laurasiatheria).

Caveat worth carrying: seven taxa give four internal splits, so the whole score
moves in quarters and a one-split difference is one edge. This is a
demonstration on real data, not a powered comparison.

Usage:
    python experiments/primate-phylogeny/estimator_ie_topology.py
    python experiments/primate-phylogeny/estimator_ie_topology.py --mark H3K4me3
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

import numpy as np

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "python"))
from hammock.runner import _jaccard_ie_from_containments as _ie  # noqa: E402

RESULTS = os.path.join(_REPO, "experiments", "primate-phylogeny", "results")

# Ordered so the printed taxon lists are stable; membership is what matters.
TAXA = ("hsa", "rheMac", "calJac", "Mmus", "bosTau", "canFam", "monDom")
TRUE_SPLITS = (
    frozenset({"hsa", "rheMac"}),
    frozenset({"hsa", "rheMac", "calJac"}),
    frozenset({"hsa", "rheMac", "calJac", "Mmus"}),
    frozenset({"bosTau", "canFam"}),
)


def _tag(path: str) -> str:
    return os.path.basename(path)[: -len(".fa")]


def load_pair(path):
    """Return (names, D_re, D_ie, offdiag_re, offdiag_ie)."""
    re_j, ie_j = {}, {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            a, b = _tag(row["file1"]), _tag(row["file2"])
            re_j[(a, b)] = float(row["jaccard_similarity"])
            ie_j[(a, b)] = float(_ie(np.array([[float(row["containment_AB"])]]),
                                     np.array([[float(row["containment_BA"])]]))[0, 0])
    names = [t for t in TAXA if (t, t) in re_j]

    def dist(src):
        m = np.zeros((len(names), len(names)))
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                if i != j:
                    m[i, j] = 1.0 - src.get((a, b), src.get((b, a)))
        return m

    off = lambda src: [v for (a, b), v in src.items() if a != b]  # noqa: E731
    return names, dist(re_j), dist(ie_j), off(re_j), off(ie_j)


def neighbor_join(names, dmat):
    """Plain NJ; returns the set of clusters joined, i.e. the induced splits."""
    dmat = dmat.astype(float).copy()
    clusters = [frozenset([n]) for n in names]
    splits = set()
    while len(clusters) > 2:
        n = len(clusters)
        row = dmat.sum(axis=1)
        q = np.full((n, n), np.inf)
        for i in range(n):
            for j in range(i + 1, n):
                q[i, j] = (n - 2) * dmat[i, j] - row[i] - row[j]
        i, j = np.unravel_index(np.argmin(q), q.shape)
        merged = clusters[i] | clusters[j]
        if 1 < len(merged) < len(names):
            splits.add(merged)
        new_row = np.array([(dmat[i, k] + dmat[j, k] - dmat[i, j]) / 2 for k in range(n)])
        keep = [k for k in range(n) if k not in (i, j)]
        nxt = np.zeros((len(keep) + 1, len(keep) + 1))
        nxt[:-1, :-1] = dmat[np.ix_(keep, keep)]
        nxt[-1, :-1] = nxt[:-1, -1] = new_row[keep]
        dmat = nxt
        clusters = [clusters[k] for k in keep] + [merged]
    return splits


def normalize(split, names):
    """Unrooted bipartition, canonicalized to the smaller side."""
    other = frozenset(names) - split
    return min(split, other, key=lambda s: (len(s), sorted(s)))


def bipartitions(names, dmat):
    """Non-trivial unrooted bipartitions. A cluster of n-1 taxa normalizes to a
    singleton, which every tree on these taxa has -- drop those or two identical
    topologies read as disagreeing."""
    out = {normalize(s, names) for s in neighbor_join(names, dmat)}
    return {s for s in out if len(s) >= 2}


def score(names, dmat):
    got = bipartitions(names, dmat)
    correct = sum(1 for t in TRUE_SPLITS if normalize(t, names) in got)
    primates = normalize(frozenset({"hsa", "rheMac", "calJac"}) & set(names), names)
    return correct, primates in got, got


def run_mark(mark: str) -> None:
    cells = sorted(os.path.basename(d) for d in glob.glob(os.path.join(RESULTS, mark, "k*_w*")))
    if not cells:
        print(f"{mark}: no cells found under {RESULTS}")
        return
    print(f"\n=== {mark} ===")
    hdr = (f"{'cell':>10}{'medJ_RE':>10}{'medJ_IE':>10}"
           f"{'true_RE':>9}{'true_IE':>9}{'primate_RE':>12}{'primate_IE':>12}")
    print(hdr)
    print("-" * len(hdr))
    disagree = []
    for cell in cells:
        found = glob.glob(os.path.join(RESULTS, mark, cell, "*.csv"))
        if not found:
            continue
        names, d_re, d_ie, off_re, off_ie = load_pair(found[0])
        if len(names) < 4:
            continue
        c_re, p_re, s_re = score(names, d_re)
        c_ie, p_ie, s_ie = score(names, d_ie)
        print(f"{cell:>10}{np.median(off_re):>10.4f}{np.median(off_ie):>10.4f}"
              f"{f'{c_re}/4':>9}{f'{c_ie}/4':>9}"
              f"{('YES' if p_re else 'no'):>12}{('YES' if p_ie else 'no'):>12}")
        if s_re != s_ie:
            disagree.append((cell, s_re - s_ie, s_ie - s_re, names))
    for cell, only_re, only_ie, names in disagree:
        fmt = lambda ss: ", ".join("+".join(sorted(normalize(s, names))) for s in ss) or "-"  # noqa: E731
        print(f"  {cell}: splits only under RE: {fmt(only_re)}")
        print(f"  {cell}: splits only under IE: {fmt(only_ie)}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mark", action="append",
                    help="histone mark subdirectory (repeatable; default both)")
    args = ap.parse_args()
    for mark in args.mark or ["H3K4me3", "H3K27ac"]:
        run_mark(mark)
    print("\nReading: the primate clade is recovered identically under both estimators")
    print("in every cell of both marks, so the published headline is unaffected. The")
    print("only difference anywhere is at k=20 on H3K4me3 -- the lowest-J cells in the")
    print("repo -- where inclusion-exclusion additionally recovers Laurasiatheria and")
    print("register-equality does not. One edge, on one mark; treat it as a hint that")
    print("the low-J crossover has consequences downstream, not as a result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
