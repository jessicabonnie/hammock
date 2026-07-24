#!/usr/bin/env python3
"""Reproduce the ChIP-Atlas cumulative experiment growth figure.

Run from the repository root:

    python paper/database_counts/scripts/plot_chip_atlas_growth.py

The script reads the transcribed source table and writes both SVG and PNG
outputs under paper/database_counts/results/.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import StrMethodFormatter


HERE = Path(__file__).resolve().parent
BASE = HERE.parent
INPUT = BASE / "source_data" / "chip_atlas_annual_experiments.tsv"
OUTPUT_DIR = BASE / "results"


def load_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t")
    required = {"year", "cumulative_experiments"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if data["year"].duplicated().any():
        raise ValueError("Years must be unique")
    if not data["year"].is_monotonic_increasing:
        raise ValueError("Years must be sorted in increasing order")
    if (data["cumulative_experiments"] < 0).any():
        raise ValueError("Experiment counts cannot be negative")
    if not data["cumulative_experiments"].is_monotonic_increasing:
        raise ValueError("Cumulative experiment counts must not decrease")

    endpoints = data.set_index("year")["cumulative_experiments"]
    expected = {2015: 37_720, 2025: 464_655}
    for year, value in expected.items():
        if year not in endpoints or int(endpoints.loc[year]) != value:
            raise ValueError(
                f"Unexpected source value for {year}; expected {value:,}. "
                "Verify against ChIP-Atlas 2025 update, Table 1."
            )

    return data


def make_plot(data: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(
        data["year"],
        data["cumulative_experiments"],
        marker="o",
        linewidth=2,
    )

    ax.set_title("Growth of ChIP-Atlas, 2015–2025")
    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative experiments")
    ax.set_xticks(data["year"])
    ax.tick_params(axis="x", rotation=45)
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax.grid(axis="y", alpha=0.25)

    events = data.set_index("year")
    ax.annotate(
        "hg38, mm10, dm6, ce11 added",
        xy=(2020, events.loc[2020, "cumulative_experiments"]),
        xytext=(2015.2, 355_000),
        arrowprops={"arrowstyle": "->", "linewidth": 1},
        fontsize=9,
    )
    ax.annotate(
        "ATAC-seq and Bisulfite-seq added",
        xy=(2021, events.loc[2021, "cumulative_experiments"]),
        xytext=(2017.7, 280_000),
        arrowprops={"arrowstyle": "->", "linewidth": 1},
        fontsize=9,
    )

    fig.text(
        0.01,
        0.01,
        "Source: ChIP-Atlas 2025 update, Table 1 (doi:10.1093/nar/gkag378).",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


def main() -> None:
    data = load_data(INPUT)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = make_plot(data)
    fig.savefig(OUTPUT_DIR / "chip_atlas_growth.svg")
    fig.savefig(OUTPUT_DIR / "chip_atlas_growth.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
