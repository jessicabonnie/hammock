from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data" / "sra_annual_experiment_counts.tsv"
OUT_SVG = ROOT / "results" / "sra_experiment_growth.svg"
OUT_PNG = ROOT / "results" / "sra_experiment_growth.png"


def validate(df: pd.DataFrame) -> None:
    required = {
        "year",
        "cumulative_experiments",
        "query_term",
        "retrieval_date",
        "query_url",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if not df["year"].is_monotonic_increasing:
        raise ValueError("Years must be increasing.")
    if (df["cumulative_experiments"].diff().fillna(0) < 0).any():
        raise ValueError("Cumulative SRA experiment count decreases.")


def make_plot(df: pd.DataFrame) -> None:
    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(
        df["year"],
        df["cumulative_experiments"],
        marker="o",
        linewidth=2,
    )
    ax.set_title("Growth of SRA experiments")
    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative experiments")
    ax.set_xticks(df["year"])
    ax.tick_params(axis="x", rotation=45)
    ax.ticklabel_format(style="plain", axis="y")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_SVG)
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(SOURCE, sep="\t")
    validate(df)
    make_plot(df)
    print(f"Wrote {OUT_SVG}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
