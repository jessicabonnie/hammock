#!/usr/bin/env python3
"""Validate and aggregate the synthetic subB scaling replicate files."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

N_VALUES = (2, 4, 8, 16, 32, 64, 128, 256, 512)
SUBB_VALUES = (1.0, 0.1, 0.01, 0.001)


def expected_reps(n: int, small_runs: int, large_runs: int) -> int:
    return small_runs if n <= 32 else large_runs


def mean(values: list[float]) -> float:
    return statistics.mean(values)


def sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--num-files", default=",".join(map(str, N_VALUES)))
    parser.add_argument("--small-runs", type=int, default=20)
    parser.add_argument("--large-runs", type=int, default=3)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    n_values = tuple(int(value) for value in args.num_files.split(","))
    if any(value not in N_VALUES for value in n_values):
        raise ValueError(f"N must be selected from {N_VALUES}")

    rows = []
    for n in n_values:
        for rep in range(expected_reps(n, args.small_runs, args.large_runs)):
            path = args.raw_dir / f"N{n}_rep{rep:02d}.csv"
            if not path.exists():
                raise RuntimeError(f"missing replicate: {path}")
            with path.open() as handle:
                current = list(csv.DictReader(handle))
            if len(current) != 9:
                raise RuntimeError(f"expected nine arms in {path}, got {len(current)}")
            rows.extend(current)

    expected_total = sum(expected_reps(n, args.small_runs, args.large_runs)
                         for n in n_values) * 9
    if len(rows) != expected_total:
        raise RuntimeError(f"expected {expected_total} rows, got {len(rows)}")

    all_fields = list(rows[0])
    with (args.output_dir / "all_runs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    cells = defaultdict(list)
    for row in rows:
        rate = "" if row["tool"] == "bedtools" else float(row["subB"])
        cells[(int(row["num_files"]), row["tool"], row["similarity"], rate)].append(row)

    summary = []
    metrics = ("wall_time", "cpu_time", "max_rss_mb", "sketch_creation_time",
               "comparison_time", "pair_time", "write_time",
               "mae_vs_bedtools", "max_err_vs_bedtools")
    for (n, tool, similarity, rate), cell in sorted(cells.items()):
        if len(cell) != expected_reps(n, args.small_runs, args.large_runs):
            raise RuntimeError(f"incomplete cell N={n}, {similarity}, subB={rate}")
        entry = {
            "num_files": n,
            "num_pairs": n * n,
            "tool": tool,
            "similarity": similarity,
            "subB": rate,
            "n_reps": len(cell),
            "num_intervals_per_file": int(cell[0]["num_intervals_per_file"]),
            "precision": int(cell[0]["precision"]),
            "threads": int(cell[0]["threads"]),
        }
        for metric in metrics:
            values = [float(row[metric]) for row in cell if row[metric] != ""]
            entry[f"mean_{metric}"] = mean(values) if values else ""
            entry[f"median_{metric}"] = statistics.median(values) if values else ""
            entry[f"sd_{metric}"] = sd(values) if values else ""
        summary.append(entry)

    summary_fields = list(summary[0])
    with (args.output_dir / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    print(f"validated {len(rows)} run rows and wrote {len(summary)} summary cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
