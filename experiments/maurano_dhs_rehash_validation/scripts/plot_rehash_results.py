#!/usr/bin/env python3
"""Render experiment-only figures after the complete evidence gates."""
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
    precision = pd.read_csv(summary / "rehash_precision_aggregates.csv")
    precision_8 = pd.read_csv(summary / "rehash_precision_aggregates_8class.csv")
    primary = (runs[runs["phase"] == "primary"] if "phase" in runs else
               runs[runs["precision"] == int(config["primary_precision"])])
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

    historical = runs[(runs.k == 10) & (runs.w == 30)]
    historical_precision = precision[(precision.k == 10) & (precision.w == 30)].sort_values("precision")
    expected_seeds = 101
    plotted_precisions = sorted(set(map(int, config["precisions"])) |
                                set(map(int, config["interpolation"]["precisions"])))
    if (historical_precision["precision"].astype(int).tolist() == plotted_precisions and
            (historical_precision["seeds_observed"] == expected_seeds).all()):
        precisions = historical_precision["precision"].astype(int).tolist()
        ari_values = [historical[historical.precision == value]["ari_10class"].to_numpy()
                      for value in precisions]
        figure, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        axes[0].boxplot(ari_values, tick_labels=precisions, showmeans=True, meanline=True)
        rng = np.random.default_rng(0)
        for position, values in enumerate(ari_values, start=1):
            axes[0].scatter(position + rng.uniform(-0.12, 0.12, len(values)), values,
                            s=12, alpha=0.25, color="black", linewidths=0)
        axes[0].set_xlabel("HLL precision p")
        axes[0].set_ylabel("ten-class ARI")
        axes[0].set_title("Seed-dependent tissue recovery")
        frequency = historical_precision["exact_partition_frequency"].to_numpy()
        lower = historical_precision["exact_partition_wilson95_low"].to_numpy()
        upper = historical_precision["exact_partition_wilson95_high"].to_numpy()
        axes[1].errorbar(precisions, frequency, yerr=[frequency - lower, upper - frequency],
                         fmt="o-", capsize=4)
        axes[1].axhline(0.5, color="0.6", linestyle="--", linewidth=1)
        axes[1].set_xticks(precisions)
        axes[1].set_ylim(0, 1)
        axes[1].set_xlabel("HLL precision p")
        axes[1].set_ylabel("exact-partition recovery")
        axes[1].set_title("101 seeds; Wilson 95% intervals")
        figure.savefig(output / "historical_precision_extension.png", dpi=180)
        plt.close(figure)
        historical_precision_8 = precision_8[
            (precision_8.k == 10) & (precision_8.w == 30)].sort_values("precision")
        if (historical_precision_8["precision"].astype(int).tolist() == precisions and
                (historical_precision_8["seeds_observed"] == expected_seeds).all()):
            ari_values_8 = [historical[historical.precision == value]["ari_8class"].to_numpy()
                            for value in precisions]
            figure, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
            axes[0].boxplot(ari_values_8, tick_labels=precisions, showmeans=True, meanline=True)
            rng = np.random.default_rng(0)
            for position, values in enumerate(ari_values_8, start=1):
                axes[0].scatter(position + rng.uniform(-0.12, 0.12, len(values)), values,
                                s=12, alpha=0.25, color="black", linewidths=0)
            axes[0].set_xlabel("HLL precision p")
            axes[0].set_ylabel("eight-class ARI")
            axes[0].set_title("Muscle-merged tissue recovery")
            frequency_8 = historical_precision_8["exact_partition_frequency"].to_numpy()
            lower_8 = historical_precision_8["exact_partition_wilson95_low"].to_numpy()
            upper_8 = historical_precision_8["exact_partition_wilson95_high"].to_numpy()
            axes[1].errorbar(precisions, frequency_8,
                             yerr=[frequency_8 - lower_8, upper_8 - frequency_8],
                             fmt="o-", capsize=4)
            axes[1].axhline(0.5, color="0.6", linestyle="--", linewidth=1)
            axes[1].set_xticks(precisions)
            axes[1].set_ylim(0, 1)
            axes[1].set_xlabel("HLL precision p")
            axes[1].set_ylabel("exact eight-cluster recovery")
            axes[1].set_title("101 paired seeds; Wilson 95% intervals")
            figure.savefig(output / "historical_precision_extension_8class.png", dpi=180)
            plt.close(figure)
        transitions = pd.read_csv(summary / "rehash_precision_transitions.csv")
        if (len(transitions) == len(plotted_precisions) - 1 and
                (transitions["paired_seeds"] == expected_seeds).all()):
            right_precision = transitions["precision_right"].astype(int)
            paired = transitions["paired_seeds"]
            figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
            axes[0].plot(right_precision, transitions["partition_changed"] / paired,
                         "o-", label="ten-class partition changed")
            axes[0].plot(right_precision, transitions["hierarchy_changed"] / paired,
                         "o-", label="ranked hierarchy changed")
            axes[0].plot(right_precision, transitions["unranked_topology_changed"] / paired,
                         "o-", label="unranked topology changed")
            axes[0].plot(right_precision,
                         transitions["ari_unchanged_hierarchy_changed"] / paired,
                         "o-", label="ranked hierarchy changed, ARI unchanged")
            axes[0].set_ylim(0, 1.03)
            axes[0].set_xticks(right_precision)
            axes[0].set_xlabel("right precision in adjacent pair")
            axes[0].set_ylabel("fraction of paired seeds")
            axes[0].set_title("Same-seed state transitions")
            axes[0].legend(fontsize=8)
            clade_line = axes[1].plot(
                right_precision, transitions["median_adjacent_clade_distance"],
                "o-", label="normalized clade distance")
            cophenetic_axis = axes[1].twinx()
            cophenetic_line = cophenetic_axis.semilogy(
                right_precision, 1 - transitions["median_adjacent_cophenetic_pearson"],
                "o-", color="tab:orange", label="1 − cophenetic Pearson r")
            axes[1].set_xticks(right_precision)
            axes[1].set_xlabel("right precision in adjacent pair")
            axes[1].set_ylabel("median normalized clade distance")
            cophenetic_axis.set_ylabel("median 1 − cophenetic Pearson r (log scale)")
            axes[1].set_title("Adjacent-precision hierarchy dispersion")
            axes[1].legend(clade_line + cophenetic_line,
                           [line.get_label() for line in clade_line + cophenetic_line], fontsize=8)
            figure.savefig(output / "adjacent_precision_hierarchy_dispersion.png", dpi=180)
            plt.close(figure)
    print(f"wrote figures to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
