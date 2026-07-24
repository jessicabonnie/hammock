#!/usr/bin/env python3
"""Plot annual growth of NCBI GEO from official fourth-quarter snapshots."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "source_data" / "geo_annual_q4_counts.tsv"
OUTPUT_SVG = ROOT / "results" / "geo_growth.svg"
OUTPUT_PNG = ROOT / "results" / "geo_growth.png"

EXPECTED_COLUMNS = ["year", "series", "platforms", "samples"]
EXPECTED_FIRST = {"year": 2001, "series": 13, "platforms": 19, "samples": 670}
EXPECTED_LAST = {
    "year": 2025,
    "series": 270_498,
    "platforms": 27_918,
    "samples": 8_216_725,
}


def validate(df: pd.DataFrame) -> None:
    """Validate the transcribed cumulative time series before plotting."""
    if list(df.columns) != EXPECTED_COLUMNS:
        raise ValueError(f"Expected columns {EXPECTED_COLUMNS}; found {list(df.columns)}")
    if df.empty:
        raise ValueError("GEO source table is empty")
    if df["year"].duplicated().any() or not df["year"].is_monotonic_increasing:
        raise ValueError("Years must be unique and strictly increasing")
    if int(df["year"].max()) > 2025:
        raise ValueError("Partial 2026 observations must not be mixed with complete years")

    for column in ["series", "platforms", "samples"]:
        if not df[column].is_monotonic_increasing:
            raise ValueError(f"Cumulative {column} counts decrease in at least one year")

    first = {key: int(df.iloc[0][key]) for key in EXPECTED_COLUMNS}
    last = {key: int(df.iloc[-1][key]) for key in EXPECTED_COLUMNS}
    if first != EXPECTED_FIRST:
        raise ValueError(f"Unexpected first row: {first}")
    if last != EXPECTED_LAST:
        raise ValueError(f"Unexpected final row: {last}")


def main() -> None:
    df = pd.read_csv(INPUT, sep="\t")
    validate(df)
    OUTPUT_SVG.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    axes[0].plot(df["year"], df["series"], marker="o", markersize=3, linewidth=2)
    axes[0].set_title("Growth of NCBI GEO, 2001–2025")
    axes[0].set_ylabel("Cumulative series")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].ticklabel_format(style="plain", axis="y")

    axes[1].plot(df["year"], df["samples"], marker="o", markersize=3, linewidth=2)
    axes[1].set_xlabel("Year (fourth-quarter snapshot)")
    axes[1].set_ylabel("Cumulative samples")
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].ticklabel_format(style="plain", axis="y")
    axes[1].set_xticks(range(2001, 2026, 2))
    axes[1].tick_params(axis="x", rotation=45)

    fig.tight_layout()
    fig.savefig(OUTPUT_SVG, bbox_inches="tight")
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    print(f"Wrote {OUTPUT_SVG}")
    print(f"Wrote {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
