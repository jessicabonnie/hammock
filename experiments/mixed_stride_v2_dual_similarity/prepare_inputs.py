#!/usr/bin/env python3
"""Summarize register-equality and join it to the existing +IE results."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


ALL_RATES = (1.0, 0.3, 0.1, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001)
PLOT_RATES = (1.0, 0.1, 0.01, 0.001)
FIELDS = (
    "similarity", "subB", "n_reps", "threads", "wall_median", "wall_mean",
    "wall_sd", "n_pairs", "mae_vs_bedtools", "max_err_vs_bedtools",
)


def load_truth(path: Path) -> dict[frozenset[str], float]:
    truth = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["file1"] != row["file2"]:
                truth[frozenset((row["file1"], row["file2"]))] = float(row["jaccard"])
    if len(truth) != 190:
        raise RuntimeError(f"expected 190 BEDTools truth pairs, got {len(truth)}")
    return truth


def summarize_re(raw_path: Path, truth: dict[frozenset[str], float]) -> list[dict]:
    walls = defaultdict(dict)
    estimates = defaultdict(dict)
    with raw_path.open() as handle:
        for row in csv.DictReader(handle):
            if row["method"] != "mixed-stride" or row["metrics"] != "False":
                raise RuntimeError("register-equality input has the wrong method or shape")
            rate = float(row["subB"])
            walls[rate][row["run_id"]] = float(row["wall_time"])
            key = frozenset((row["file_a"], row["file_b"]))
            if key not in truth:
                continue
            value = float(row["jaccard"])
            previous = estimates[rate].setdefault(key, value)
            if previous != value:
                raise RuntimeError(f"non-deterministic RE estimate at subB={rate}, pair={key}")

    output = []
    for rate in ALL_RATES:
        cell_walls = list(walls[rate].values())
        cell = estimates[rate]
        if len(cell_walls) != 5 or cell.keys() != truth.keys():
            raise RuntimeError(
                f"incomplete RE cell subB={rate}: runs={len(cell_walls)}, pairs={len(cell)}"
            )
        errors = [abs(cell[key] - truth[key]) for key in truth]
        output.append({
            "similarity": "register-equality",
            "subB": rate,
            "n_reps": len(cell_walls),
            "threads": 8,
            "wall_median": statistics.median(cell_walls),
            "wall_mean": statistics.mean(cell_walls),
            "wall_sd": statistics.stdev(cell_walls),
            "n_pairs": len(cell),
            "mae_vs_bedtools": statistics.mean(errors),
            "max_err_vs_bedtools": max(errors),
        })
    return output


def load_ie(path: Path) -> list[dict]:
    by_rate = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            rate = float(row["subB"])
            by_rate[rate] = {
                "similarity": "inclusion-exclusion",
                "subB": rate,
                "n_reps": int(row["n_reps"]),
                "threads": int(row["threads"]),
                "wall_median": float(row["wall_median"]),
                "wall_mean": float(row["wall_mean"]),
                "wall_sd": float(row["wall_sd"]),
                "n_pairs": int(row["n_pairs"]),
                "mae_vs_bedtools": float(row["mae_ie_vs_bedtools"]),
                "max_err_vs_bedtools": float(row["max_err_ie_vs_bedtools"]),
            }
    if set(by_rate) != set(ALL_RATES):
        raise RuntimeError(f"unexpected +IE rate grid: {sorted(by_rate)}")
    for rate, row in by_rate.items():
        if row["n_reps"] != 5 or row["threads"] != 8 or row["n_pairs"] != 190:
            raise RuntimeError(f"incompatible +IE cell at subB={rate}: {row}")
    return [by_rate[rate] for rate in ALL_RATES]


def write_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--register-equality", type=Path, required=True)
    parser.add_argument("--inclusion-exclusion", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    truth = load_truth(args.truth)
    re_rows = summarize_re(args.register_equality, truth)
    ie_rows = load_ie(args.inclusion_exclusion)
    write_rows(args.output_dir / "register_equality_summary.csv", re_rows)
    selected = [row for rate in PLOT_RATES for row in re_rows + ie_rows
                if row["subB"] == rate]
    write_rows(args.output_dir / "dual_similarity_plot.csv", selected)
    print(f"wrote summaries under {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
