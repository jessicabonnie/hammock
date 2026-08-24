#!/usr/bin/env python3
"""Score exact and rehashed matrices under the frozen ten-class rules."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, cut_tree, dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import binomtest, chi2, kendalltau, pearsonr, spearmanr
from common import EXPERIMENT, atomic_json, fasta_paths, grid, load_config, resolve, sha256, strip_fasta


def choose2(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sum(values * (values - 1) / 2))


def adjusted_rand_score(labels_a, labels_b) -> float:
    table = pd.crosstab(pd.Series(labels_a, dtype="object"),
                        pd.Series(labels_b, dtype="object")).to_numpy()
    n = int(table.sum())
    index = choose2(table)
    rows, columns = choose2(table.sum(axis=1)), choose2(table.sum(axis=0))
    expected = rows * columns / math.comb(n, 2)
    maximum = (rows + columns) / 2
    return 0.0 if math.isclose(maximum, expected) else (index - expected) / (maximum - expected)


def normalized_mutual_info_score(labels_a, labels_b) -> float:
    table = pd.crosstab(pd.Series(labels_a, dtype="object"),
                        pd.Series(labels_b, dtype="object")).to_numpy(dtype=float)
    total = table.sum()
    row_p, column_p, joint_p = table.sum(axis=1) / total, table.sum(axis=0) / total, table / total
    product = np.outer(row_p, column_p)
    mask = joint_p > 0
    mutual_information = float(np.sum(joint_p[mask] * np.log(joint_p[mask] / product[mask])))
    row_h = float(-np.sum(row_p[row_p > 0] * np.log(row_p[row_p > 0])))
    column_h = float(-np.sum(column_p[column_p > 0] * np.log(column_p[column_p > 0])))
    return 0.0 if row_h == 0 or column_h == 0 else 2 * mutual_information / (row_h + column_h)


def tissue_labels(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t")
    if frame["File"].map(strip_fasta).duplicated().any():
        raise ValueError("filename key is not one-to-one")
    return dict(zip(frame["File"].map(strip_fasta), frame["Biosample_term_name"]))


def frame_matrix(frame: pd.DataFrame, metric: str, samples: list[str]) -> np.ndarray:
    index = {sample: i for i, sample in enumerate(samples)}
    matrix = np.full((len(samples), len(samples)), np.nan)
    for row in frame.itertuples(index=False):
        left, right = strip_fasta(row.file1), strip_fasta(row.file2)
        if left not in index or right not in index:
            raise ValueError(f"unknown sample in matrix: {left}, {right}")
        matrix[index[left], index[right]] = float(getattr(row, metric))
    if np.isnan(matrix).any():
        raise ValueError("matrix is incomplete")
    excursions = int(np.count_nonzero((matrix < 0) | (matrix > 1)))
    if np.any(matrix < -1e-12) or np.any(matrix > 1 + 1e-12):
        raise ValueError("matrix contains undocumented excursions")
    matrix = np.clip(matrix, 0, 1)
    if not np.allclose(matrix, matrix.T, rtol=1e-9, atol=1e-9):
        raise ValueError("matrix is asymmetric")
    if not np.allclose(np.diag(matrix), 1, rtol=0, atol=1e-12):
        raise ValueError("matrix diagonal is not one")
    return matrix


def cluster(matrix: np.ndarray, samples: list[str], truth: list[str], artifact_dir: Path,
            stem: str) -> dict[str, object]:
    hierarchy = linkage(squareform(1 - matrix, checks=False), method="average")
    partition = cut_tree(hierarchy, n_clusters=10).ravel()
    leaf_order = dendrogram(hierarchy, no_plot=True)["leaves"]
    merges_before_cut = len(samples) - 10
    lower = float(hierarchy[merges_before_cut - 1, 2]) if merges_before_cut else 0.0
    upper = float(hierarchy[merges_before_cut, 2])
    artifact_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"sample": samples, "tissue": truth, "cluster": partition}).to_csv(
        artifact_dir / f"{stem}_partition.csv", index=False)
    pd.DataFrame(hierarchy, columns=["left", "right", "height", "size"]).to_csv(
        artifact_dir / f"{stem}_linkage.csv", index=False)
    atomic_json(artifact_dir / f"{stem}_cluster.json", {
        "leaf_order": [samples[i] for i in leaf_order], "cut_clusters": 10,
        "height_below_cut": lower, "height_above_cut": upper, "cut_gap": upper - lower})
    return {"partition": partition, "hierarchy": hierarchy,
            "ari_10class": adjusted_rand_score(truth, partition),
            "nmi_10class": normalized_mutual_info_score(truth, partition),
            "cut_height_lower": lower, "cut_height_upper": upper, "cut_gap": upper - lower}


def unique_upper(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices_from(matrix, k=1)]


def errors(estimate: np.ndarray, target: np.ndarray) -> dict[str, float]:
    observed, expected = unique_upper(estimate), unique_upper(target)
    delta = observed - expected
    return {"signed_error": float(delta.mean()), "mae": float(np.abs(delta).mean()),
            "rmse": float(np.sqrt(np.mean(delta ** 2))), "max_abs_error": float(np.abs(delta).max()),
            "pearson": float(pearsonr(observed, expected).statistic),
            "spearman": float(spearmanr(observed, expected).statistic),
            "kendall": float(kendalltau(observed, expected).statistic)}


def hierarchy_clade_sequence(hierarchy: np.ndarray, sample_count: int) -> list[frozenset[int]]:
    nodes = {index: frozenset([index]) for index in range(sample_count)}
    clades = []
    for offset, row in enumerate(hierarchy):
        clade = nodes[int(row[0])] | nodes[int(row[1])]
        nodes[sample_count + offset] = clade
        if len(clade) < sample_count:
            clades.append(clade)
    return clades


def hierarchy_clades(hierarchy: np.ndarray, sample_count: int) -> set[frozenset[int]]:
    return set(hierarchy_clade_sequence(hierarchy, sample_count))


def hierarchy_agreement(estimate: np.ndarray, target: np.ndarray,
                        sample_count: int) -> dict[str, float | str]:
    estimate_sequence = hierarchy_clade_sequence(estimate, sample_count)
    estimate_clades = set(estimate_sequence)
    target_clades = hierarchy_clades(target, sample_count)
    denominator = max(sample_count - 2, 1)
    unranked_signature = ";".join(
        ",".join(map(str, sorted(clade)))
        for clade in sorted(estimate_clades, key=lambda value: (len(value), tuple(sorted(value)))))
    ranked_signature = ";".join(",".join(map(str, sorted(clade))) for clade in estimate_sequence)
    estimate_cophenetic, target_cophenetic = cophenet(estimate), cophenet(target)
    return {
        "hierarchy_signature": hashlib.sha256(ranked_signature.encode()).hexdigest()[:16],
        "unranked_topology_signature": hashlib.sha256(unranked_signature.encode()).hexdigest()[:16],
        "clade_distance_vs_exact": 1 - len(estimate_clades & target_clades) / denominator,
        "cophenetic_pearson_vs_exact": float(pearsonr(estimate_cophenetic, target_cophenetic).statistic),
        "cophenetic_spearman_vs_exact": float(spearmanr(estimate_cophenetic, target_cophenetic).statistic),
        "cophenetic_mae_vs_exact": float(np.abs(estimate_cophenetic - target_cophenetic).mean()),
    }


def completion_phase(metadata: dict, config: dict) -> str:
    """Recover phase for old manifests and use explicit phase for new ones."""
    if metadata.get("phase"):
        return str(metadata["phase"])
    if int(metadata["seed"]) in set(map(int, config["seeds"])):
        return ("primary" if int(metadata["precision"]) == int(config["primary_precision"])
                else "followup")
    return "extension"


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return math.nan, math.nan
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    half_width = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return centre - half_width, centre + half_width


def extension_run_count(config: dict) -> int:
    requested = set(range(int(config["extension"]["seed_start"]),
                          int(config["extension"]["seed_stop"]) + 1))
    requested.update(map(int, config["extension"].get("additional_seeds", [])))
    new_seeds = requested - set(map(int, config["seeds"]))
    return len(config["extension"]["precisions"]) * len(config["extension"]["cells"]) * len(new_seeds)


def interpolation_run_count(config: dict) -> int:
    requested = set(range(int(config["interpolation"]["seed_start"]),
                          int(config["interpolation"]["seed_stop"]) + 1))
    requested.update(map(int, config["interpolation"].get("additional_seeds", [])))
    return (len(config["interpolation"]["precisions"]) *
            len(config["interpolation"]["cells"]) * len(requested))


def expected_phase_keys(config: dict, phase: str) -> set[tuple[int, int, int, int]]:
    spec = config[phase]
    seeds = set(range(int(spec["seed_start"]), int(spec["seed_stop"]) + 1))
    seeds.update(map(int, spec.get("additional_seeds", [])))
    if phase == "extension":
        seeds -= set(map(int, config["seeds"]))
    return {(int(precision), int(seed), int(cell["k"]), int(cell["w"]))
            for precision in spec["precisions"] for seed in seeds for cell in spec["cells"]}


def cochrans_q(values: np.ndarray) -> tuple[float, float]:
    matrix = np.asarray(values, dtype=int)
    treatments = matrix.shape[1]
    columns = matrix.sum(axis=0)
    rows = matrix.sum(axis=1)
    denominator = treatments * rows.sum() - np.sum(rows ** 2)
    if denominator == 0:
        return 0.0, 1.0
    statistic = ((treatments - 1) *
                 (treatments * np.sum(columns ** 2) - columns.sum() ** 2) / denominator)
    return float(statistic), float(chi2.sf(statistic, treatments - 1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    args = parser.parse_args()
    config, base = load_config(args.config)
    samples = [strip_fasta(str(path)) for path in fasta_paths(config, base)]
    tissue_for = tissue_labels(resolve(base, config["inputs"]["filename_key"]))
    if set(samples) != set(tissue_for):
        raise SystemExit(f"sample/key mismatch: fasta-only={sorted(set(samples)-set(tissue_for))}, key-only={sorted(set(tissue_for)-set(samples))}")
    truth = [tissue_for[sample] for sample in samples]
    summary_dir = EXPERIMENT / "results" / "summaries"
    cluster_dir = summary_dir / "clusters"
    summary_dir.mkdir(parents=True, exist_ok=True)

    bed_frame = pd.read_csv(resolve(base, config["inputs"]["bedtools_reference"]), sep="\t")
    bed_matrix = frame_matrix(bed_frame.rename(columns={"jaccard": "bedtools_jaccard"}),
                              "bedtools_jaccard", samples)
    exact_rows, exact_matrices, exact_partitions, exact_hierarchies = [], {}, {}, {}
    for k, w in grid(config):
        path = EXPERIMENT / "results" / "exact" / f"k{k}_w{w}" / "pairs.csv"
        if not path.is_file():
            continue
        matrix = frame_matrix(pd.read_csv(path), "jaccard_exact", samples)
        scores = cluster(matrix, samples, truth, cluster_dir, f"exact_k{k}_w{w}")
        exact_matrices[(k, w)] = matrix
        exact_partitions[(k, w)] = scores.pop("partition")
        exact_hierarchies[(k, w)] = scores.pop("hierarchy")
        exact_rows.append({"k": k, "w": w, **scores,
                           **{f"bedtools_{key}": value for key, value in errors(matrix, bed_matrix).items()}})
    exact_scores = pd.DataFrame(exact_rows)
    if not exact_scores.empty:
        exact_scores.sort_values(["ari_10class", "nmi_10class", "k", "w"],
                                 ascending=[False, False, True, True]).to_csv(
                                     summary_dir / "exact_scores.csv", index=False)

    run_rows, run_hierarchies, seen_run_keys = [], {}, set()
    observed_primary_outputs = 0
    for completion in sorted((EXPERIMENT / "results" / "rehash").glob("p*/seed*/*.complete.json")):
        metadata = json.loads(completion.read_text())
        if metadata.get("status") != "complete":
            continue
        if metadata.get("hash_scheme") != config["hash_scheme"]:
            raise ValueError(f"unexpected hash scheme in {completion}")
        phase = completion_phase(metadata, config)
        if phase == "primary":
            observed_primary_outputs += 1
        k, w = int(metadata["k"]), int(metadata["w"])
        if (k, w) not in exact_matrices:
            continue
        recorded_output = Path(metadata["output"])
        output = completion.with_name(recorded_output.name)
        if not output.is_file():
            output = recorded_output
        if metadata.get("sha256") != sha256(output):
            raise ValueError(f"completion checksum mismatch: {output}")
        run_key = (int(metadata["precision"]), int(metadata["seed"]), k, w)
        if run_key in seen_run_keys:
            raise ValueError(f"duplicate precision/seed/cell output: {run_key}")
        seen_run_keys.add(run_key)
        matrix = frame_matrix(pd.read_csv(output), "jaccard_similarity_ie", samples)
        stem = f"rehash_p{metadata['precision']}_seed{int(metadata['seed']):05d}_k{k}_w{w}"
        scores = cluster(matrix, samples, truth, cluster_dir, stem)
        partition = scores.pop("partition")
        hierarchy = scores.pop("hierarchy")
        run_hierarchies[run_key] = hierarchy
        co_clustering = np.equal.outer(partition, partition)
        exact_error = errors(matrix, exact_matrices[(k, w)])
        bed_error = errors(matrix, bed_matrix)
        hierarchy_scores = hierarchy_agreement(hierarchy, exact_hierarchies[(k, w)], len(samples))
        run_rows.append({"phase": phase, "hash_scheme": metadata["hash_scheme"], "source_commit": metadata["source_commit"],
                         "precision": int(metadata["precision"]), "seed": int(metadata["seed"]),
                         "k": k, "w": w, **scores,
                         "partition_signature": hashlib.sha256(co_clustering.tobytes()).hexdigest()[:16],
                         "exact_partition_equal": bool(np.array_equal(
                             co_clustering,
                             np.equal.outer(exact_partitions[(k, w)], exact_partitions[(k, w)]))),
                         "partition_ari_vs_exact": adjusted_rand_score(exact_partitions[(k, w)], partition),
                         **hierarchy_scores,
                         **{f"exact_{key}": value for key, value in exact_error.items()},
                         **{f"bedtools_{key}": value for key, value in bed_error.items()},
                         "elapsed_s": metadata.get("elapsed_s"), "output": str(output)})
    run_scores = pd.DataFrame(run_rows)
    if not run_scores.empty:
        run_scores.sort_values(["precision", "k", "w", "seed"]).to_csv(
            summary_dir / "rehash_run_scores.csv", index=False)

    precision_aggregates = pd.DataFrame()
    if not run_scores.empty:
        precision_aggregates = run_scores.groupby(["precision", "k", "w"], as_index=False).agg(
            seeds_observed=("seed", "nunique"), mean_ari=("ari_10class", "mean"),
            std_ari=("ari_10class", "std"), q25_ari=("ari_10class", lambda values: values.quantile(0.25)),
            median_ari=("ari_10class", "median"),
            q75_ari=("ari_10class", lambda values: values.quantile(0.75)),
            min_ari=("ari_10class", "min"), max_ari=("ari_10class", "max"),
            mean_nmi=("nmi_10class", "mean"), std_nmi=("nmi_10class", "std"),
            q25_nmi=("nmi_10class", lambda values: values.quantile(0.25)),
            median_nmi=("nmi_10class", "median"),
            q75_nmi=("nmi_10class", lambda values: values.quantile(0.75)),
            min_nmi=("nmi_10class", "min"),
            max_nmi=("nmi_10class", "max"),
            distinct_ari_values=("ari_10class", "nunique"),
            distinct_partitions=("partition_signature", "nunique"),
            distinct_hierarchies=("hierarchy_signature", "nunique"),
            distinct_unranked_topologies=("unranked_topology_signature", "nunique"),
            median_clade_distance_vs_exact=("clade_distance_vs_exact", "median"),
            max_clade_distance_vs_exact=("clade_distance_vs_exact", "max"),
            median_cophenetic_pearson_vs_exact=("cophenetic_pearson_vs_exact", "median"),
            min_cophenetic_pearson_vs_exact=("cophenetic_pearson_vs_exact", "min"),
            exact_partition_count=("exact_partition_equal", "sum"),
            exact_partition_frequency=("exact_partition_equal", "mean"),
            min_partition_ari_vs_exact=("partition_ari_vs_exact", "min"),
            median_exact_mae=("exact_mae", "median"), max_exact_mae=("exact_mae", "max"),
            max_exact_max_abs_error=("exact_max_abs_error", "max"))
        intervals = [wilson_interval(int(row.exact_partition_count), int(row.seeds_observed))
                     for row in precision_aggregates.itertuples()]
        precision_aggregates["exact_partition_wilson95_low"] = [bounds[0] for bounds in intervals]
        precision_aggregates["exact_partition_wilson95_high"] = [bounds[1] for bounds in intervals]
        precision_aggregates.to_csv(summary_dir / "rehash_precision_aggregates.csv", index=False)
        outcome_counts = run_scores.groupby(
            ["precision", "k", "w", "ari_10class", "nmi_10class", "partition_signature"],
            as_index=False).agg(seed_count=("seed", "nunique"),
                                exact_partition=("exact_partition_equal", "all"))
        outcome_counts.to_csv(summary_dir / "rehash_precision_outcomes.csv", index=False)

        historical = run_scores[(run_scores.k == 10) & (run_scores.w == 30)]
        transition_rows = []
        precision_values = sorted(map(int, historical["precision"].unique()))
        for left_precision, right_precision in zip(precision_values, precision_values[1:]):
            left = historical[historical.precision == left_precision].set_index("seed")
            right = historical[historical.precision == right_precision].set_index("seed")
            common_seeds = sorted(set(left.index) & set(right.index))
            if not common_seeds:
                continue
            left, right = left.loc[common_seeds], right.loc[common_seeds]
            ari_delta = right["ari_10class"] - left["ari_10class"]
            nmi_delta = right["nmi_10class"] - left["nmi_10class"]
            exact_gain = int((~left["exact_partition_equal"] & right["exact_partition_equal"]).sum())
            exact_loss = int((left["exact_partition_equal"] & ~right["exact_partition_equal"]).sum())
            discordant = exact_gain + exact_loss
            adjacent_hierarchy = [
                hierarchy_agreement(
                    run_hierarchies[(int(right_precision), int(seed), 10, 30)],
                    run_hierarchies[(int(left_precision), int(seed), 10, 30)], len(samples))
                for seed in common_seeds]
            hierarchy_changed = left["hierarchy_signature"] != right["hierarchy_signature"]
            topology_changed = (left["unranked_topology_signature"] !=
                                right["unranked_topology_signature"])
            ari_unchanged = np.isclose(left["ari_10class"], right["ari_10class"], rtol=0, atol=1e-15)
            transition_rows.append({
                "precision_left": int(left_precision), "precision_right": int(right_precision),
                "paired_seeds": len(common_seeds), "ari_improved": int((ari_delta > 1e-15).sum()),
                "ari_declined": int((ari_delta < -1e-15).sum()),
                "ari_unchanged": int(np.isclose(ari_delta, 0, rtol=0, atol=1e-15).sum()),
                "mean_ari_delta": float(ari_delta.mean()), "median_ari_delta": float(ari_delta.median()),
                "mean_nmi_delta": float(nmi_delta.mean()), "median_nmi_delta": float(nmi_delta.median()),
                "exact_partition_gain": exact_gain, "exact_partition_loss": exact_loss,
                "mcnemar_exact_p": float(binomtest(exact_gain, discordant, 0.5).pvalue)
                if discordant else 1.0,
                "partition_changed": int((left["partition_signature"] != right["partition_signature"]).sum()),
                "hierarchy_changed": int(hierarchy_changed.sum()),
                "unranked_topology_changed": int(topology_changed.sum()),
                "ari_unchanged_hierarchy_changed": int((ari_unchanged & hierarchy_changed).sum()),
                "ari_unchanged_unranked_topology_changed": int(
                    (ari_unchanged & topology_changed).sum()),
                "median_adjacent_clade_distance": float(np.median(
                    [item["clade_distance_vs_exact"] for item in adjacent_hierarchy])),
                "median_adjacent_cophenetic_pearson": float(np.median(
                    [item["cophenetic_pearson_vs_exact"] for item in adjacent_hierarchy])),
                "median_adjacent_cophenetic_mae": float(np.median(
                    [item["cophenetic_mae_vs_exact"] for item in adjacent_hierarchy])),
            })
        pd.DataFrame(transition_rows).to_csv(summary_dir / "rehash_precision_transitions.csv", index=False)
        exact_pivot = historical.pivot(index="seed", columns="precision", values="exact_partition_equal")
        repeated_measures = {"status": "incomplete", "precisions": precision_values,
                             "complete_paired_seeds": int(exact_pivot.dropna().shape[0])}
        if len(precision_values) > 1 and not exact_pivot.isna().any().any():
            statistic, pvalue = cochrans_q(exact_pivot.to_numpy())
            repeated_measures.update({"status": "complete", "cochrans_q": statistic,
                                      "degrees_of_freedom": len(precision_values) - 1,
                                      "pvalue": pvalue})
        atomic_json(summary_dir / "precision_repeated_measures.json", repeated_measures)

    aggregates = pd.DataFrame()
    primary = run_scores[run_scores["phase"] == "primary"] if not run_scores.empty else run_scores
    if not primary.empty:
        aggregates = primary.groupby(["k", "w"], as_index=False).agg(
            seeds_observed=("seed", "nunique"), median_ari=("ari_10class", "median"),
            min_ari=("ari_10class", "min"), max_ari=("ari_10class", "max"),
            exact_partition_frequency=("exact_partition_equal", "mean"),
            median_nmi=("nmi_10class", "median"), median_exact_mae=("exact_mae", "median"),
            median_bedtools_spearman=("bedtools_spearman", "median"),
            median_bedtools_mae=("bedtools_mae", "median"), median_cut_gap=("cut_gap", "median"))
        aggregates.to_csv(summary_dir / "rehash_cell_aggregates.csv", index=False)

    leaders = {"status": "incomplete", "exact_leaders": [], "robustness_leaders": [],
               "required_exact_cells": 37, "observed_exact_cells": len(exact_scores),
               "required_primary_runs": 296, "observed_primary_outputs": observed_primary_outputs,
               "observed_primary_runs_scored": len(primary)}
    decision = "# Decision report\n\n**Classification: unresolved.** The complete exact and eight-seed evidence gates have not passed.\n"
    if len(exact_scores) == 37 and len(primary) == 296 and (aggregates["seeds_observed"] == 8).all():
        maximum_ari = exact_scores["ari_10class"].max()
        ari_ties = exact_scores[np.isclose(exact_scores["ari_10class"], maximum_ari, rtol=0, atol=1e-15)]
        maximum_nmi = ari_ties["nmi_10class"].max()
        exact_leaders = ari_ties[np.isclose(ari_ties["nmi_10class"], maximum_nmi, rtol=0, atol=1e-15)]
        ranked = aggregates.sort_values(
            ["median_ari", "min_ari", "exact_partition_frequency", "median_nmi", "median_exact_mae", "k", "w"],
            ascending=[False, False, False, False, True, True, True])
        leader_row = ranked.iloc[0]
        robust = ranked[(ranked["median_ari"] == leader_row["median_ari"]) &
                        (ranked["min_ari"] == leader_row["min_ari"]) &
                        (ranked["exact_partition_frequency"] == leader_row["exact_partition_frequency"]) &
                        (ranked["median_nmi"] == leader_row["median_nmi"]) &
                        (ranked["median_exact_mae"] == leader_row["median_exact_mae"])]
        leaders = {"status": "frozen", "exact_leaders": exact_leaders[["k", "w"]].to_dict("records"),
                   "robustness_leaders": robust[["k", "w"]].to_dict("records"),
                   "required_exact_cells": 37, "observed_exact_cells": 37,
                   "required_primary_runs": 296, "observed_primary_outputs": 296,
                   "observed_primary_runs_scored": 296}
        historical_exact = any((int(row.k), int(row.w)) == (10, 30) for row in exact_leaders.itertuples())
        followup_cells = {(10, 30),
                          *((int(row.k), int(row.w)) for row in exact_leaders.itertuples()),
                          *((int(row.k), int(row.w)) for row in robust.itertuples())}
        followup_precisions = {int(value) for value in config["precisions"]
                               if int(value) != int(config["primary_precision"])}
        followup = run_scores[
            (run_scores["phase"] == "followup") &
            run_scores["precision"].isin(followup_precisions) &
            pd.Series(list(zip(run_scores["k"], run_scores["w"])), index=run_scores.index).isin(followup_cells)]
        expected_followup_runs = len(followup_cells) * len(followup_precisions) * len(config["seeds"])
        followup_groups = followup.groupby(["precision", "k", "w"])["seed"].nunique()
        followup_complete = (len(followup) == expected_followup_runs and
                             len(followup_groups) == len(followup_cells) * len(followup_precisions) and
                             (followup_groups == len(config["seeds"])).all())

        historical_precision = precision_aggregates[
            (precision_aggregates.k == 10) & (precision_aggregates.w == 30)].sort_values("precision")
        table_lines = ["| p | seeds | ARI median [IQR] | ARI range | NMI median [IQR] | partitions / ranked hierarchies / topologies | median clade distance | exact partition (95% CI) | median exact MAE |",
                       "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for row in historical_precision.itertuples():
            table_lines.append(
                f"| {int(row.precision)} | {int(row.seeds_observed)} | {row.median_ari:.6f} "
                f"[{row.q25_ari:.6f}, {row.q75_ari:.6f}] | "
                f"{row.min_ari:.6f}–{row.max_ari:.6f} | {row.median_nmi:.6f} "
                f"[{row.q25_nmi:.6f}, {row.q75_nmi:.6f}] | {int(row.distinct_partitions)} / "
                f"{int(row.distinct_hierarchies)} / {int(row.distinct_unranked_topologies)} | "
                f"{row.median_clade_distance_vs_exact:.3f} | "
                f"{int(row.exact_partition_count)}/{int(row.seeds_observed)} "
                f"({row.exact_partition_wilson95_low:.3f}–{row.exact_partition_wilson95_high:.3f}) | "
                f"{row.median_exact_mae:.8g} |")

        extension_expected_keys = expected_phase_keys(config, "extension")
        extension_observed_keys = set(map(tuple, run_scores[run_scores["phase"] == "extension"][
            ["precision", "seed", "k", "w"]].astype(int).to_numpy()))
        extension_expected, extension_observed = len(extension_expected_keys), len(extension_observed_keys)
        extension_complete = extension_observed_keys == extension_expected_keys
        extension_text = (f"The exploratory seed extension is complete ({extension_observed}/{extension_expected} "
                          "new runs; 101 total seeds per precision)." if extension_complete else
                          f"The exploratory seed extension is incomplete ({extension_observed}/{extension_expected} new runs).")
        interpolation_expected_keys = expected_phase_keys(config, "interpolation")
        interpolation_observed_keys = set(map(tuple, run_scores[run_scores["phase"] == "interpolation"][
            ["precision", "seed", "k", "w"]].astype(int).to_numpy()))
        interpolation_expected = len(interpolation_expected_keys)
        interpolation_observed = len(interpolation_observed_keys)
        interpolation_complete = interpolation_observed_keys == interpolation_expected_keys
        interpolation_text = (
            f"The missing-precision interpolation is complete ({interpolation_observed}/{interpolation_expected} runs)."
            if interpolation_complete else
            f"The missing-precision interpolation is incomplete ({interpolation_observed}/{interpolation_expected} runs).")

        if not followup_complete:
            classification = "unresolved"
            gate_text = (f"Precision follow-up is incomplete: observed {len(followup)} of "
                         f"{expected_followup_runs} required non-primary runs.")
        elif len(exact_leaders) > 1 and historical_exact:
            classification = "part of an exact plateau"
            gate_text = "All prespecified exact, p=18 seed, and precision-follow-up gates passed."
        elif not historical_exact:
            classification = "not the exact-feature optimum"
            gate_text = "All prespecified exact, p=18 seed, and precision-follow-up gates passed."
        else:
            historical_runs = run_scores[(run_scores.phase.isin(["primary", "followup"])) &
                                         (run_scores.k == 10) & (run_scores.w == 30)]
            estimator_stable = (len(historical_runs) == len(config["precisions"]) * len(config["seeds"]) and
                                historical_runs["exact_partition_equal"].all())
            classification = "exact and robust" if estimator_stable else "exact but estimator-sensitive"
            gate_text = "All prespecified exact, p=18 seed, and precision-follow-up gates passed."
        decision = ("# Decision report\n\n"
                    f"**Classification: {classification}.**\n\n"
                    f"{gate_text} The historical `k=10,w=30` cell is "
                    f"{'an' if historical_exact else 'not an'} exact-feature leader. "
                    f"{extension_text} {interpolation_text}\n\n"
                    "## Historical-cell precision evidence\n\n" + "\n".join(table_lines) +
                    "\n\nTissue rankings are exploratory because selection and recovery use the same "
                    "20 labelled samples. No compatibility default or manuscript text was changed.\n")
    atomic_json(summary_dir / "leaders.json", leaders)
    (summary_dir / "DECISION.md").write_text(decision)
    print(json.dumps(leaders, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
