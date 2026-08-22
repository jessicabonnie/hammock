#!/usr/bin/env python3
"""Render experiment-only figures after the complete p=18 evidence gate."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import EXPERIMENT, grid, load_config


def grid_heatmap(frame: pd.DataFrame, value: str, title: str, label: str, output: Path,
                 config) -> None:
    ks = [int(item) for item in config["k_values"]]
    ws = [int(item) for item in config["w_values"]]
    matrix = np.full((len(ks), len(ws)), np.nan)
    for row in frame.itertuples(index=False):
        matrix[ks.index(int(row.k)), ws.index(int(row.w))] = float(getattr(row, value))
    masked = np.ma.masked_invalid(matrix)
    figure, axis = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
    image = axis.imshow(masked, aspect="auto", interpolation="none", cmap="viridis")
    axis.set_xticks(range(len(ws)), ws)
    axis.set_yticks(range(len(ks)), ks)
    axis.set_xlabel("window size w")
    axis.set_ylabel("k-mer size k")
    axis.set_title(title)
    figure.colorbar(image, ax=axis, label=label)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    summary = EXPERIMENT / "results" / "summaries"
    exact = pd.read_csv(summary / "exact_scores.csv")
    runs = pd.read_csv(summary / "rehash_run_scores.csv")
    aggregates = pd.read_csv(summary / "rehash_cell_aggregates.csv")
    primary = runs[runs["precision"] == int(config["primary_precision"])]
    if len(exact) != 37 or len(primary) != 296 or len(aggregates) != 37:
        raise SystemExit("refusing partial figures: expected 37 exact cells and 296 p=18 runs")
    output = EXPERIMENT / "figures"
    output.mkdir(parents=True, exist_ok=True)
    grid_heatmap(exact, "ari_10class", "Exact selected-feature tissue recovery",
                 "ten-class ARI", output / "exact_ari_heatmap.png", config)
    grid_heatmap(aggregates, "exact_partition_frequency",
                 "Exact-partition recovery over eight rehash seeds", "recovery frequency",
                 output / "rehash_exact_partition_recovery.png", config)
    grid_heatmap(aggregates, "median_exact_mae",
                 "Rehashed HLL error against exact selected features", "median MAE",
                 output / "rehash_exact_mae.png", config)

    ordered = grid(config)
    labels = [f"{k},{w}" for k, w in ordered]
    values = [primary[(primary.k == k) & (primary.w == w)]["ari_10class"].to_numpy()
              for k, w in ordered]
    figure, axis = plt.subplots(figsize=(14, 5), constrained_layout=True)
    axis.boxplot(values, tick_labels=labels, showmeans=True, meanline=True)
    axis.tick_params(axis="x", labelrotation=90, labelsize=7)
    axis.set_xlabel("(k,w)")
    axis.set_ylabel("ten-class ARI")
    axis.set_title("Tissue recovery across eight prespecified rehash seeds")
    figure.savefig(output / "rehash_ari_distributions.png", dpi=180)
    plt.close(figure)

    merged = exact[["k", "w", "cut_gap"]].merge(
        aggregates[["k", "w", "median_cut_gap"]], on=["k", "w"])
    figure, axis = plt.subplots(figsize=(6, 5), constrained_layout=True)
    axis.scatter(merged["cut_gap"], merged["median_cut_gap"], c=merged["k"], cmap="viridis")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("exact ten-cluster cut gap")
    axis.set_ylabel("median rehash cut gap")
    axis.set_title("Exact and rehashed clustering margins")
    figure.savefig(output / "cut_margin_comparison.png", dpi=180)
    plt.close(figure)
    print(f"wrote five figures to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
