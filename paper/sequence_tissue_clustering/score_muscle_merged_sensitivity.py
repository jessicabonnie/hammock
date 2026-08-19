#!/usr/bin/env python3
"""Rescore the Mode D sweep after merging arm/back/leg into one muscle class."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import squareform


FILE_RE = re.compile(
    r"hammock_mnmzr_p(?P<p>\d+)_jaccD_k(?P<k>\d+)_w(?P<w>\d+)"
    r"(?P<full>_full)?\.csv$"
)


def strip_ext(value: str) -> str:
    name = Path(value).name
    return re.sub(r"\.(bed|fa|fasta|fna)$", "", name, flags=re.IGNORECASE)


def choose2(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sum(values * (values - 1) / 2))


def adjusted_rand(labels_a: list[str], labels_b: np.ndarray) -> float:
    table = pd.crosstab(pd.Series(labels_a), pd.Series(labels_b)).to_numpy()
    n = table.sum()
    index = choose2(table)
    rows = choose2(table.sum(axis=1))
    cols = choose2(table.sum(axis=0))
    expected = rows * cols / math.comb(int(n), 2)
    maximum = (rows + cols) / 2
    return 0.0 if math.isclose(maximum, expected) else (index - expected) / (maximum - expected)


def normalized_mi(labels_a: list[str], labels_b: np.ndarray) -> float:
    table = pd.crosstab(pd.Series(labels_a), pd.Series(labels_b)).to_numpy(dtype=float)
    n = table.sum()
    p_i = table.sum(axis=1) / n
    p_j = table.sum(axis=0) / n
    p_ij = table / n
    outer = np.outer(p_i, p_j)
    mask = p_ij > 0
    mutual_info = float(np.sum(p_ij[mask] * np.log(p_ij[mask] / outer[mask])))
    h_i = float(-np.sum(p_i[p_i > 0] * np.log(p_i[p_i > 0])))
    h_j = float(-np.sum(p_j[p_j > 0] * np.log(p_j[p_j > 0])))
    return 0.0 if h_i == 0 or h_j == 0 else 2 * mutual_info / (h_i + h_j)


def add_ie(frame: pd.DataFrame) -> pd.DataFrame:
    if "jaccard_similarity_ie" in frame:
        return frame
    required = {"containment_AB", "containment_BA"}
    if not required.issubset(frame.columns):
        raise ValueError("input lacks both native IE and directional containments")
    c_ab = np.minimum(frame["containment_AB"].astype(float), 1.0)
    c_ba = np.minimum(frame["containment_BA"].astype(float), 1.0)
    frame = frame.copy()
    frame["jaccard_similarity_ie"] = np.where(
        (c_ab <= 0) | (c_ba <= 0), 0.0, 1.0 / (1.0 / c_ab + 1.0 / c_ba - 1.0)
    )
    return frame


def similarity_matrix(frame: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    stems1 = frame["file1"].map(strip_ext)
    stems2 = frame["file2"].map(strip_ext)
    stems = sorted(set(stems1) | set(stems2))
    index = {stem: i for i, stem in enumerate(stems)}
    matrix = np.full((len(stems), len(stems)), np.nan)
    for stem1, stem2, value in zip(stems1, stems2, frame["jaccard_similarity_ie"]):
        matrix[index[stem1], index[stem2]] = value
    matrix = np.where(np.isnan(matrix), matrix.T, matrix)
    if np.isnan(matrix).any():
        raise ValueError("similarity matrix is incomplete")
    np.fill_diagonal(matrix, 1.0)
    return stems, np.clip(matrix, 0.0, 1.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=18)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--key", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    experiment = repo / "experiments" / "maurano_dhs_validation"
    raw_dir = args.raw_dir or experiment / "results" / "raw_d"
    key_path = args.key or experiment / "data" / "maurano_filenames_key.tsv"
    output = args.output or Path(__file__).with_name(
        f"muscle_merged_p{args.precision}_sensitivity.csv"
    )

    key = pd.read_csv(key_path, sep="\t")
    tissue_for = dict(zip(key["File"].map(strip_ext), key["Biosample_term_name"]))

    # Prefer native, current `_full` outputs when both archived and current
    # files exist for a cell.
    chosen: dict[tuple[int, int], tuple[bool, Path]] = {}
    for path in raw_dir.glob(f"hammock_mnmzr_p{args.precision}_jaccD_k*_w*.csv"):
        match = FILE_RE.fullmatch(path.name)
        if not match or int(match["p"]) != args.precision:
            continue
        cell = (int(match["k"]), int(match["w"]))
        is_full = bool(match["full"])
        if cell not in chosen or (is_full and not chosen[cell][0]):
            chosen[cell] = (is_full, path)

    rows = []
    for (k, w), (_, path) in sorted(chosen.items()):
        frame = add_ie(pd.read_csv(path))
        stems, matrix = similarity_matrix(frame)
        hierarchy = linkage(squareform(1.0 - matrix, checks=False), method="average")

        truth_10 = [tissue_for[stem] for stem in stems]
        truth_8 = ["fMuscle" if label.startswith("fMuscle_") else label for label in truth_10]
        prediction_10 = cut_tree(hierarchy, n_clusters=10).ravel()
        prediction_8 = cut_tree(hierarchy, n_clusters=8).ravel()
        rows.append(
            {
                "precision": args.precision,
                "k": k,
                "w": w,
                "input": path.name,
                "ari_10class": adjusted_rand(truth_10, prediction_10),
                "nmi_10class": normalized_mi(truth_10, prediction_10),
                "ari_8class": adjusted_rand(truth_8, prediction_8),
                "nmi_8class": normalized_mi(truth_8, prediction_8),
            }
        )

    result = pd.DataFrame(rows).sort_values(
        ["ari_8class", "nmi_8class", "k", "w"],
        ascending=[False, False, True, True],
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)

    best_ari = result[result["ari_8class"] == result["ari_8class"].max()]
    best_nmi = result[result["nmi_8class"] == result["nmi_8class"].max()]
    print(f"Scored {len(result)} cells at p={args.precision}")
    print("Maximum eight-class ARI:")
    print(best_ari[["k", "w", "ari_8class", "nmi_8class"]].to_string(index=False))
    print("Maximum eight-class NMI:")
    print(best_nmi[["k", "w", "ari_8class", "nmi_8class"]].to_string(index=False))
    print(f"Wrote: {output}")


if __name__ == "__main__":
    main()
