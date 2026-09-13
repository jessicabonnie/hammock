#!/usr/bin/env python3
"""Plot the archived mixed-stride versus hash-threshold speed result."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SUMMARY = SCRIPT_DIR / "results" / "job_30459651" / "summary.csv"
DEFAULT_OUTPUT = SCRIPT_DIR / "figures" / "mixed_stride_vs_hash_threshold_speed"


def read_column(rows: list[dict[str, str]], name: str) -> list[float]:
    return [float(row[name]) for row in rows]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    with args.summary.open() as handle:
        rows = sorted(csv.DictReader(handle), key=lambda row: float(row["subB"]))
    if not rows:
        raise RuntimeError(f"no results in {args.summary}")

    rates = read_column(rows, "subB")
    series = (
        (
            "Sketch construction",
            "paired_sketch_speedup_median",
            "paired_sketch_speedup_ci95_low",
            "paired_sketch_speedup_ci95_high",
            "#0072B2",
            "o",
        ),
        (
            "End-to-end wall time",
            "paired_wall_speedup_median",
            "paired_wall_speedup_ci95_low",
            "paired_wall_speedup_ci95_high",
            "#D55E00",
            "s",
        ),
    )

    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8.5,
        "svg.fonttype": "none",
    })
    figure, (time_axis, speed_axis) = plt.subplots(
        1, 2, figsize=(8.0, 3.45), sharex=True, constrained_layout=True
    )
    time_axis.plot(
        rates,
        read_column(rows, "mixed_stride_sketch_median_s"),
        color="#0072B2",
        marker="o",
        markersize=5,
        linewidth=1.6,
        label="Mixed-stride",
    )
    time_axis.plot(
        rates,
        read_column(rows, "hash_threshold_sketch_median_s"),
        color="#CC79A7",
        marker="^",
        markersize=5,
        linewidth=1.6,
        label="Hash-threshold",
    )
    for label, middle, low, high, color, marker in series:
        values = read_column(rows, middle)
        lows = read_column(rows, low)
        highs = read_column(rows, high)
        errors = (
            [value - bound for value, bound in zip(values, lows)],
            [bound - value for value, bound in zip(values, highs)],
        )
        speed_axis.errorbar(
            rates,
            values,
            yerr=errors,
            color=color,
            marker=marker,
            markersize=5,
            linewidth=1.6,
            capsize=2.5,
            label=label,
        )

    time_axis.set_ylabel("Median sketch-construction time (s)")
    time_axis.set_ylim(0, 9)
    time_axis.set_title("A", loc="left", fontweight="bold")
    time_axis.legend(frameon=False, loc="upper left")

    speed_axis.axhline(1, color="#666666", linewidth=1, linestyle="--", zorder=0)
    speed_axis.set_ylabel("Speedup over hash-threshold (fold)")
    speed_axis.set_ylim(0.95, 2.4)
    speed_axis.set_title("B", loc="left", fontweight="bold")
    speed_axis.legend(frameon=False, loc="upper right")

    for axis in (time_axis, speed_axis):
        axis.set_xscale("log")
        axis.set_xticks(rates)
        axis.set_xticklabels([f"{rate:g}" for rate in rates], rotation=35, ha="right")
        axis.set_xlabel("Sampling fraction (subB)")
        axis.grid(axis="y", color="#dddddd", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg_path = args.output_prefix.with_suffix(".svg")
    figure.savefig(svg_path)
    figure.savefig(args.output_prefix.with_suffix(".png"), dpi=300)
    plt.close(figure)
    # Matplotlib writes trailing spaces in SVG path data. Normalize them so
    # repository-wide whitespace validation remains useful.
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n"
    )
    print(svg_path)
    print(args.output_prefix.with_suffix(".png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
