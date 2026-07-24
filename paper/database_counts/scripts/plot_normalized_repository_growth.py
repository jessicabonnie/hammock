from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CHIP_GEO_SOURCE = (
    ROOT / "source_data" / "chip_atlas_geo_comparison_2015_2025.tsv"
)
SRA_SOURCE = ROOT / "source_data" / "sra_annual_experiment_counts.tsv"
OUT_SVG = ROOT / "results" / "normalized_repository_growth.svg"
OUT_PNG = ROOT / "results" / "normalized_repository_growth.png"
BASE_YEAR = 2015


def load_and_join() -> pd.DataFrame:
    chip_geo = pd.read_csv(CHIP_GEO_SOURCE, sep="\t")
    sra = pd.read_csv(SRA_SOURCE, sep="\t")

    required_sra = {"year", "cumulative_experiments"}
    missing = required_sra.difference(sra.columns)
    if missing:
        raise ValueError(f"Missing SRA columns: {sorted(missing)}")

    sra = sra[["year", "cumulative_experiments"]].rename(
        columns={"cumulative_experiments": "sra_experiments"}
    )
    joined = chip_geo.merge(sra, on="year", how="inner")
    joined = joined[joined["year"] >= BASE_YEAR].copy()

    if joined.empty or int(joined.iloc[0]["year"]) != BASE_YEAR:
        raise ValueError(f"All three series must contain the base year {BASE_YEAR}.")

    columns = ["chip_atlas_experiments", "geo_series", "sra_experiments"]
    for column in columns:
        if (joined[column] <= 0).any():
            raise ValueError(f"Values in {column} must be positive.")
        if (joined[column].diff().fillna(0) < 0).any():
            raise ValueError(f"Cumulative series decreases in {column}.")
        joined[f"{column}_index"] = joined[column] / joined.iloc[0][column] * 100

    return joined


def make_plot(df: pd.DataFrame) -> None:
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.8, 5.2))

    ax.plot(
        df["year"],
        df["chip_atlas_experiments_index"],
        marker="o",
        linewidth=2,
        label="ChIP-Atlas experiments",
    )
    ax.plot(
        df["year"],
        df["geo_series_index"],
        marker="o",
        linewidth=2,
        label="GEO Series",
    )
    ax.plot(
        df["year"],
        df["sra_experiments_index"],
        marker="o",
        linewidth=2,
        label="SRA experiments",
    )

    ax.axhline(100, linewidth=1, linestyle="--")
    ax.set_title(f"Relative repository growth, {BASE_YEAR} = 100")
    ax.set_xlabel("Year")
    ax.set_ylabel(f"Growth index ({BASE_YEAR} = 100)")
    ax.set_xticks(df["year"])
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_SVG)
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)


def main() -> None:
    dataframe = load_and_join()
    make_plot(dataframe)
    print(f"Wrote {OUT_SVG}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
