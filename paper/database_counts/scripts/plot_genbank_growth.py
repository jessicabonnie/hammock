from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data" / "genbank_august_growth_2009_2025.tsv"
OUT_SVG = ROOT / "results" / "genbank_growth.svg"
OUT_PNG = ROOT / "results" / "genbank_growth.png"


def validate(df: pd.DataFrame) -> None:
    expected_columns = [
        "year",
        "release",
        "genbank_bases",
        "genbank_sequences",
        "wgs_bases",
        "wgs_sequences",
        "total_bases",
        "total_sequences",
    ]
    if list(df.columns) != expected_columns:
        raise ValueError(f"Unexpected columns: {list(df.columns)}")
    if not df["year"].is_monotonic_increasing:
        raise ValueError("Years must be increasing.")
    if int(df.iloc[0]["year"]) != 2009 or int(df.iloc[-1]["year"]) != 2025:
        raise ValueError("Expected August snapshots from 2009 through 2025.")
    if not (df["total_bases"] == df["genbank_bases"] + df["wgs_bases"]).all():
        raise ValueError("Stored total_bases do not equal GenBank plus WGS bases.")
    if not (
        df["total_sequences"]
        == df["genbank_sequences"] + df["wgs_sequences"]
    ).all():
        raise ValueError("Stored total_sequences do not equal GenBank plus WGS sequences.")
    for column in ["total_bases", "total_sequences"]:
        if (df[column].diff().fillna(0) < 0).any():
            raise ValueError(f"Derived cumulative series decreases in {column}.")


def make_plot(df: pd.DataFrame) -> None:
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2), sharex=True)

    axes[0].plot(df["year"], df["total_sequences"], marker="o", linewidth=2)
    axes[0].set_ylabel("Sequence records")
    axes[0].set_yscale("log")
    axes[0].grid(axis="y", alpha=0.25)

    axes[1].plot(df["year"], df["total_bases"], marker="o", linewidth=2)
    axes[1].set_ylabel("Base pairs")
    axes[1].set_xlabel("Year")
    axes[1].set_yscale("log")
    axes[1].grid(axis="y", alpha=0.25)

    axes[1].set_xticks(df["year"])
    axes[1].tick_params(axis="x", rotation=45)
    fig.suptitle("Growth of GenBank and WGS, August releases 2009–2025")
    fig.tight_layout()
    fig.savefig(OUT_SVG)
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)


def main() -> None:
    dataframe = pd.read_csv(SOURCE, sep="\t")
    validate(dataframe)
    make_plot(dataframe)
    print(f"Wrote {OUT_SVG}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
