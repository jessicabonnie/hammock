#!/usr/bin/env python3
"""Prepare Figure 3 plotting inputs from the reproduction's raw outputs."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def load_truth(path: Path) -> dict[frozenset[str], float]:
    truth = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["file1"] != row["file2"]:
                truth[frozenset((row["file1"], row["file2"]))] = float(row["jaccard"])
    if len(truth) != 190:
        raise RuntimeError(f"expected 190 bedtools truth pairs, got {len(truth)}")
    return truth


def summarize_panel_b(hammock_csv: Path, truth_path: Path,
                      output: Path, figure3_output: Path) -> None:
    truth = load_truth(truth_path)
    walls = defaultdict(dict)
    estimates = defaultdict(dict)
    with hammock_csv.open() as handle:
        for row in csv.DictReader(handle):
            if row["method"] != "mixed-stride":
                raise RuntimeError(f"unexpected panel B method: {row['method']}")
            rate = float(row["subB"])
            walls[rate][row["run_id"]] = float(row["wall_time"])
            key = frozenset((row["file_a"], row["file_b"]))
            if key not in truth:
                continue
            value = float(row["jaccard_ie"])
            previous = estimates[rate].setdefault(key, value)
            if previous != value:
                raise RuntimeError(f"non-deterministic estimate at subB={rate}, pair={key}")

    output_rows = []
    rates = (1.0, 0.3, 0.1, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001)
    for rate in rates:
        cell_walls = list(walls[rate].values())
        cell = estimates[rate]
        if len(cell_walls) != 5 or cell.keys() != truth.keys():
            raise RuntimeError(
                f"incomplete panel B cell subB={rate}: "
                f"runs={len(cell_walls)}, pairs={len(cell)}")
        errors = [abs(cell[key] - truth[key]) for key in truth]
        output_rows.append({
            "method": "mixed-stride", "size_class": "maurano", "subB": rate,
            "n_reps": len(cell_walls), "threads": 8,
            "wall_median": statistics.median(cell_walls),
            "wall_mean": statistics.mean(cell_walls),
            "wall_sd": statistics.stdev(cell_walls), "n_pairs": len(cell),
            "mae_ie_vs_bedtools": statistics.mean(errors),
            "max_err_ie_vs_bedtools": max(errors),
        })
    fields = list(output_rows[0])
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    with figure3_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(row for row in output_rows
                         if row["subB"] in {1.0, 0.1, 0.01, 0.001})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-dir", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    args = parser.parse_args()
    results = args.experiment_dir / "results"
    prepared = results / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    summarize_panel_b(
        results / "panel_b" / "hammock.csv",
        args.truth,
        prepared / "panel_b_summary.csv",
        prepared / "panel_b_figure3_summary.csv",
    )
    print(f"wrote prepared inputs under {prepared}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
