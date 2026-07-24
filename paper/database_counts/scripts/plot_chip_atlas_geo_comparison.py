from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data" / "chip_atlas_geo_comparison_2015_2025.tsv"
OUT_SVG = ROOT / "results" / "chip_atlas_geo_comparison.svg"
OUT_PNG = ROOT / "results" / "chip_atlas_geo_comparison.png"


def validate(df: pd.DataFrame) -> None:
    expected_columns = ["year", "chip_atlas_experiments", "geo_series"]
    if list(df.columns) != expected_columns:
        raise ValueError(f"Unexpected columns: {list(df.columns)}")
    if not df["year"].is_monotonic_increasing:
        raise ValueError("Years must be increasing.")
    for column in ["chip_atlas_experiments", "geo_series"]:
        if (df[column].diff().fillna(0) < 0).any():
            raise ValueError(f"Cumulative series decreases in {column}.")
    if int(df.iloc[0]["year"]) != 2015 or int(df.iloc[-1]["year"]) != 2025:
        raise ValueError("Expected overlapping comparison window 2015–2025.")


def make_plot(df: pd.DataFrame) -> None:
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(
        df["year"],
        df["chip_atlas_experiments"],
        marker="o",
        linewidth=2,
        label="ChIP-Atlas experiments",
    )
    ax.plot(
        df["year"],
        df["geo_series"],
        marker="o",
        linewidth=2,
        label="GEO Series",
    )
    ax.set_title("Growth of ChIP-Atlas and GEO, 2015–2025")
    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative records")
    ax.set_xticks(df["year"])
    ax.tick_params(axis="x", rotation=45)
    ax.ticklabel_format(style="plain", axis="y")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    ax.annotate(
        "ATAC-seq added to ChIP-Atlas",
        xy=(2021, 196136),
        xytext=(2017.8, 275000),
        arrowprops={"arrowstyle": "->", "lw": 1},
        fontsize=9,
    )
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
