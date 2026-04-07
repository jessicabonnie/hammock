#!/usr/bin/env python3
"""
scripts/plot_sweep.py
=====================
Visualize sweep_kw.py output: bias and variance of minimizer sketch
relative to true k-mer Jaccard, across (k, w) combinations.

Usage:
    python scripts/plot_sweep.py results/sweep_*.csv
    python scripts/plot_sweep.py results/sweep_*.csv --out results/figs/
"""

import argparse
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load(paths):
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)


def plot_bias_heatmap(df, col, ax, title):
    """Mean signed error (bias) as a heatmap over k × w."""
    pivot = df.groupby(["k", "w"])[col].mean().unstack("w")
    sns.heatmap(pivot, ax=ax, center=0, cmap="RdBu_r", annot=True, fmt=".3f",
                cbar_kws={"label": "mean error (sketch − true)"})
    ax.set_title(title)
    ax.set_xlabel("w")
    ax.set_ylabel("k")


def plot_mae_heatmap(df, col, ax, title):
    """Mean absolute error as a heatmap over k × w."""
    pivot = df.groupby(["k", "w"])[col].apply(lambda x: x.abs().mean()).unstack("w")
    sns.heatmap(pivot, ax=ax, cmap="YlOrRd", annot=True, fmt=".3f",
                cbar_kws={"label": "MAE"})
    ax.set_title(title)
    ax.set_xlabel("w")
    ax.set_ylabel("k")


def plot_by_mutation(df, col, out_dir):
    """Line plot: sketched vs true Jaccard at each mutation rate, per (k,w)."""
    df["kw"] = "k=" + df["k"].astype(str) + " w=" + df["w"].astype(str)
    mean = df.groupby(["kw", "mutation_rate"])[
        ["true_kmer_jaccard", col]
    ].mean().reset_index()

    fig, ax = plt.subplots(figsize=(10, 5))
    palette = sns.color_palette("tab20", n_colors=mean["kw"].nunique())
    for (kw, grp), color in zip(mean.groupby("kw"), palette):
        ax.plot(grp["true_kmer_jaccard"], grp[col], marker="o", label=kw, color=color)
    # Identity line
    lo = df["true_kmer_jaccard"].min()
    hi = df["true_kmer_jaccard"].max()
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="identity")
    ax.set_xlabel("true k-mer Jaccard")
    ax.set_ylabel(f"sketched ({col})")
    ax.set_title(f"Sketch vs true Jaccard by (k, w)\n[{col}]")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
    fig.tight_layout()
    path = os.path.join(out_dir, f"sketch_vs_true_{col}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  saved {path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csvs", nargs="+", help="sweep CSV file(s)")
    p.add_argument("--out", default="results/figs", help="output directory")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    df = load(args.csvs)

    # Drop rows with NaN errors (degenerate sketch cases)
    df = df.dropna(subset=["error_jaccard", "error_jaccard_with_ends"])

    for col, err_col in [
        ("sketched_jaccard",           "error_jaccard"),
        ("sketched_jaccard_with_ends", "error_jaccard_with_ends"),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        plot_bias_heatmap(df, err_col, axes[0], f"Bias [{col}]")
        plot_mae_heatmap( df, err_col, axes[1], f"MAE  [{col}]")
        fig.tight_layout()
        path = os.path.join(args.out, f"heatmap_{col}.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  saved {path}")

        plot_by_mutation(df, col, args.out)

    print("Done.")


if __name__ == "__main__":
    main()
