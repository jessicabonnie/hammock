#!/usr/bin/env python3
"""Check whether average-linkage clustering recovers one clade per sample in every sweep matrix."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import squareform

from analyze import PEAK_TYPES, SEEDS, parse_label, read_cells, validate_matrix


MIN_K_CLAIM = 8  # The manuscript claims clade recovery at every configuration with k >= 8.


def clade_stats(rows: list[dict[str, str]]) -> dict[str, float | bool]:
    labels = sorted({parse_label(row["file1"]) for row in rows})
    if len(labels) != 9 or len({sample for sample, _ in labels}) != 3:
        raise ValueError(f"expected 3 samples on 3 references, found {labels}")
    index = {label: i for i, label in enumerate(labels)}
    sim = np.empty((9, 9))
    for row in rows:
        sim[index[parse_label(row["file1"])], index[parse_label(row["file2"])]] = float(
            row["jaccard_similarity_ie"])
    # validate_matrix allows asymmetry up to 1e-12, which squareform would reject.
    dist = 1 - (sim + sim.T) / 2
    np.fill_diagonal(dist, 0)
    tree = linkage(squareform(dist), method="average")
    clusters = cut_tree(tree, n_clusters=3).ravel()
    partition = {frozenset(i for i in range(9) if clusters[i] == c) for c in set(clusters)}
    expected = {
        frozenset(i for label, i in index.items() if label[0] == sample)
        for sample, _ in labels
    }
    same = [sim[i, j] for a, i in index.items() for b, j in index.items()
            if a[0] == b[0] and a[1] != b[1]]
    diff = [sim[i, j] for a, i in index.items() for b, j in index.items() if a[0] != b[0]]
    heights = tree[:, 2]
    return {
        "clades_recovered": partition == expected,
        "separated_ie": min(same) > max(diff),
        "min_same_sample": min(same),
        "max_different_sample": max(diff),
        # Tied merges can make the recovered tree depend on the implementation.
        "tied_merge_heights": any(
            math.isclose(a, b, abs_tol=1e-12) for a, b in zip(heights, heights[1:])),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--precision", type=int, required=True)
    args = parser.parse_args()

    rows_out: list[dict[str, object]] = []
    for seed in SEEDS:
        for peak_type in PEAK_TYPES:
            for k, w in read_cells(args.cells):
                path = (args.results / f"seed{seed}" / peak_type / f"k{k}_w{w}" /
                        f"exp_a_mnmzr_p{args.precision}_jaccD_k{k}_w{w}_rehash-selector64_full.csv")
                if not path.is_file():
                    raise FileNotFoundError(path)
                rows = validate_matrix(path, k, w, args.precision)
                rows_out.append({"seed": seed, "peak_type": peak_type, "cell": f"k{k}_w{w}",
                                 "k": k, "w": w, **clade_stats(rows)})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_out[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)

    claim = [r for r in rows_out if int(r["k"]) >= MIN_K_CLAIM]
    rest = [r for r in rows_out if int(r["k"]) < MIN_K_CLAIM]
    n_claim = sum(bool(r["clades_recovered"]) for r in claim)
    print(f"k>={MIN_K_CLAIM}: clades recovered in {n_claim}/{len(claim)} matrices")
    print(f"k<{MIN_K_CLAIM}: clades recovered in "
          f"{sum(bool(r['clades_recovered']) for r in rest)}/{len(rest)} matrices; "
          f"{sum(bool(r['tied_merge_heights']) for r in rest)} have tied merge heights")
    print(f"wrote {args.output}")
    return 0 if n_claim == len(claim) else 1


if __name__ == "__main__":
    raise SystemExit(main())
