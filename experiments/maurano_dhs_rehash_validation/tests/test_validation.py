from __future__ import annotations

import csv
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from validation import METRICS, validate_hll_csv
from run_rehash_sweep import extension_seeds, normalize_csv_lf, parse_time_report
from analyze_rehash_sweep import completion_phase, extension_run_count, wilson_interval
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
