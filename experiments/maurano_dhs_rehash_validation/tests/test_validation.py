from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from validation import METRICS, validate_hll_csv
from run_rehash_sweep import (extension_seeds, normalize_csv_lf, parse_time_report,
                              plan_fingerprint, requested_seeds)
from analyze_rehash_sweep import (cochrans_q, completion_phase, expected_phase_keys, holm_adjust,
                                  extension_run_count, hierarchy_agreement, hierarchy_clades,
                                  interpolation_run_count, muscle_merged_labels, wilson_interval)
import numpy as np
from common import strip_fasta


def test_named_column_validator_accepts_complete_directional_matrix(tmp_path):
    path = tmp_path / "matrix.csv"
    fields = ["file1", "file2", "sketch_type", "mode", "precision", "num_hashes",
              "kmer_size", "window_size", *METRICS, "ref1", "ref2"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for left in ("a", "b"):
            for right in ("a", "b"):
                same = left == right
                row = {"file1": left + ".fa", "file2": right + ".fa", "sketch_type": "minimizer",
                       "mode": "D", "precision": 18, "num_hashes": "NA", "kmer_size": 10,
                       "window_size": 30, "ref1": "NA", "ref2": "NA"}
                row.update({metric: 1.0 if same else 0.25 for metric in METRICS})
                writer.writerow(row)
    result = validate_hll_csv(path, ["a", "b"], k=10, w=30, precision=18)
    assert result["rows"] == 4
    assert result["samples"] == 2


def test_gnu_time_report_parser(tmp_path):
    path = tmp_path / "time.txt"
    path.write_text("\n".join([
        "User time (seconds): 4.25", "System time (seconds): 0.75",
        "Maximum resident set size (kbytes): 12345",
    ]))
    assert parse_time_report(path) == {
        "user_cpu_s": 4.25, "system_cpu_s": 0.75,
        "maximum_rss_kb": 12345.0, "cpu_s": 5.0,
    }


def test_csv_line_ending_normalization(tmp_path):
    path = tmp_path / "matrix.csv"
    path.write_bytes(b"a,b\r\n1,2\r\n")
    assert normalize_csv_lf(path)
    assert path.read_bytes() == b"a,b\n1,2\n"
    assert not normalize_csv_lf(path)


def test_sample_stem_preserves_biological_merge_suffix():
    stem = "fBrain-DS14718.hotspot.twopass.fdr0.05.merge"
    assert strip_fasta(stem) == stem
    assert strip_fasta(stem + ".fa") == stem
    assert strip_fasta(stem + ".bed") == stem


def test_eight_class_endpoint_only_merges_fetal_muscle_labels():
    labels = ["fMuscle_arm", "fMuscle_back", "fMuscle_leg", "fHeart", "Muscle"]
    assert muscle_merged_labels(labels) == ["fMuscle", "fMuscle", "fMuscle", "fHeart", "Muscle"]


def test_extension_adds_93_seeds_without_changing_frozen_eight():
    config = {"seeds": [1, 2, 3, 17, 42, 43, 99, 31337],
              "extension": {"seed_start": 0, "seed_stop": 99,
                            "additional_seeds": [31337]}}
    seeds = extension_seeds(config)
    assert len(seeds) == 93
    assert 0 in seeds and 98 in seeds
    assert not set(seeds) & set(config["seeds"])


def test_completion_phase_preserves_old_manifests():
    config = {"seeds": [1, 2], "primary_precision": 18}
    assert completion_phase({"seed": 1, "precision": 18}, config) == "primary"
    assert completion_phase({"seed": 1, "precision": 24}, config) == "followup"
    assert completion_phase({"seed": 0, "precision": 18}, config) == "extension"
    assert completion_phase({"phase": "extension", "seed": 1, "precision": 18}, config) == "extension"


def test_wilson_interval_contains_observed_proportion():
    low, high = wilson_interval(50, 101)
    assert low < 50 / 101 < high


def test_extension_run_count_excludes_the_frozen_seeds():
    config = {"seeds": [1, 2, 3, 17, 42, 43, 99, 31337],
              "extension": {"seed_start": 0, "seed_stop": 99,
                            "additional_seeds": [31337],
                            "precisions": [12, 18, 22, 24], "cells": [{"k": 10, "w": 30}]}}
    assert extension_run_count(config) == 372


def test_interpolation_uses_all_101_seeds_at_nine_precisions():
    spec = {"seed_start": 0, "seed_stop": 99, "additional_seeds": [31337],
            "precisions": [13, 14, 15, 16, 17, 19, 20, 21, 23],
            "cells": [{"k": 10, "w": 30}]}
    config = {"interpolation": spec}
    assert len(requested_seeds(spec)) == 101
    assert interpolation_run_count(config) == 909


def test_hierarchy_agreement_detects_clade_reorganization():
    left_first = np.array([[0, 1, 0.1, 2], [2, 3, 0.2, 2], [4, 5, 0.5, 4]], dtype=float)
    crossed = np.array([[0, 2, 0.1, 2], [1, 3, 0.2, 2], [4, 5, 0.5, 4]], dtype=float)
    assert len(hierarchy_clades(left_first, 4)) == 2
    same = hierarchy_agreement(left_first, left_first, 4)
    changed = hierarchy_agreement(crossed, left_first, 4)
    assert same["clade_distance_vs_exact"] == 0
    assert np.isclose(same["cophenetic_pearson_vs_exact"], 1)
    assert changed["clade_distance_vs_exact"] == 1
    reordered = left_first[[1, 0, 2]].copy()
    reordered[0, 0:2] = [2, 3]
    reordered[1, 0:2] = [0, 1]
    reordered[:, 2] = [0.1, 0.2, 0.5]
    ranked = hierarchy_agreement(reordered, left_first, 4)
    assert ranked["clade_distance_vs_exact"] == 0
    assert ranked["hierarchy_signature"] != same["hierarchy_signature"]
    assert ranked["unranked_topology_signature"] == same["unranked_topology_signature"]


def test_plan_fingerprint_changes_when_an_index_identity_changes():
    job = {"phase": "interpolation", "precision": 13, "seed": 0,
           "output": "/tmp/result.csv", "command": ["hammock", "--precision", "13"]}
    changed = {**job, "precision": 14}
    assert plan_fingerprint([job]) != plan_fingerprint([changed])


def test_expected_interpolation_keys_are_exact_and_cochrans_q_is_paired():
    config = {"seeds": [1], "interpolation": {
        "seed_start": 0, "seed_stop": 1, "additional_seeds": [31337],
        "precisions": [13, 14], "cells": [{"k": 10, "w": 30}]}}
    keys = expected_phase_keys(config, "interpolation")
    assert len(keys) == 6
    assert (14, 31337, 10, 30) in keys
    statistic, pvalue = cochrans_q(np.array([[0, 1], [0, 1], [1, 1]]))
    assert np.isclose(statistic, 2)
    assert 0 < pvalue < 1
    json.dumps({"precisions": sorted(map(int, np.array([12, 13]))),
                "cochrans_q": statistic, "pvalue": pvalue})
    assert np.allclose(holm_adjust([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
