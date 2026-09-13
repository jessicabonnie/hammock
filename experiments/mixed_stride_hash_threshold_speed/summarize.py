#!/usr/bin/env python3
"""Validate and summarize a mixed-stride/hash-threshold timing run."""

from __future__ import annotations

import argparse
import csv
import random
import statistics
from collections import defaultdict
from pathlib import Path


METHODS = ("mixed-stride", "hash-threshold")
METRICS = ("sketch_seconds", "wall_seconds", "cpu_seconds")


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def bootstrap_ci(values: list[float], seed: int, draws: int = 10000) -> tuple[float, float]:
    rng = random.Random(seed)
    medians = [
        statistics.median(rng.choices(values, k=len(values)))
        for _ in range(draws)
    ]
    return percentile(medians, 0.025), percentile(medians, 0.975)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    raw_path = args.run_dir / "raw_runs.csv"
    with raw_path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"no rows in {raw_path}")

    by_cell: dict[tuple[float, str], list[dict[str, str]]] = defaultdict(list)
    by_round: dict[tuple[float, int], dict[str, dict[str, str]]] = defaultdict(dict)
    total_points = {int(row["total_points"]) for row in rows}
    if len(total_points) != 1:
        raise RuntimeError(f"total-point count changed across runs: {sorted(total_points)}")
    for row in rows:
        rate = float(row["subB"])
        method = row["method"]
        round_index = int(row["round"])
        if method not in METHODS:
            raise RuntimeError(f"unexpected method {method}")
        by_cell[(rate, method)].append(row)
        if method in by_round[(rate, round_index)]:
            raise RuntimeError(f"duplicate rate/round/method: {rate}/{round_index}/{method}")
        by_round[(rate, round_index)][method] = row

    rates = sorted({rate for rate, _ in by_cell}, reverse=True)
    rounds = sorted({int(row["round"]) for row in rows})
    output_rows = []
    for rate_index, rate in enumerate(rates):
        for method in METHODS:
            cell = by_cell[(rate, method)]
            if len(cell) != len(rounds):
                raise RuntimeError(
                    f"incomplete cell subB={rate}, method={method}: "
                    f"{len(cell)} rows for {len(rounds)} rounds"
                )
            digests = {row["output_sha256"] for row in cell}
            if len(digests) != 1:
                raise RuntimeError(f"non-deterministic output for subB={rate}, method={method}")
        paired = [by_round[(rate, round_index)] for round_index in rounds]
        if any(set(pair) != set(METHODS) for pair in paired):
            raise RuntimeError(f"unpaired round at subB={rate}")

        row: dict[str, object] = {
            "subB": rate,
            "rounds": len(rounds),
            "mixed_stride_realized_fraction": statistics.median(
                float(item["realized_fraction"])
                for item in by_cell[(rate, "mixed-stride")]
            ),
            "hash_threshold_realized_fraction": statistics.median(
                float(item["realized_fraction"])
                for item in by_cell[(rate, "hash-threshold")]
            ),
        }
        for metric_index, metric in enumerate(METRICS):
            mixed = [float(item[metric]) for item in by_cell[(rate, "mixed-stride")]]
            hashed = [float(item[metric]) for item in by_cell[(rate, "hash-threshold")]]
            ratios = [
                float(pair["hash-threshold"][metric]) /
                float(pair["mixed-stride"][metric])
                for pair in paired
            ]
            low, high = bootstrap_ci(
                ratios, seed=20260901 + rate_index * len(METRICS) + metric_index
            )
            stem = metric.removesuffix("_seconds")
            row[f"mixed_stride_{stem}_median_s"] = statistics.median(mixed)
            row[f"hash_threshold_{stem}_median_s"] = statistics.median(hashed)
            row[f"paired_{stem}_speedup_median"] = statistics.median(ratios)
            row[f"paired_{stem}_speedup_ci95_low"] = low
            row[f"paired_{stem}_speedup_ci95_high"] = high
        output_rows.append(row)

    summary_path = args.run_dir / "summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    print("subB  sketch speedup [95% CI]  wall speedup [95% CI]")
    for row in output_rows:
        print(
            f"{row['subB']:<5g} "
            f"{row['paired_sketch_speedup_median']:.3f} "
            f"[{row['paired_sketch_speedup_ci95_low']:.3f}, "
            f"{row['paired_sketch_speedup_ci95_high']:.3f}]  "
            f"{row['paired_wall_speedup_median']:.3f} "
            f"[{row['paired_wall_speedup_ci95_low']:.3f}, "
            f"{row['paired_wall_speedup_ci95_high']:.3f}]"
        )
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
