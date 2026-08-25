#!/usr/bin/env python3
"""Named-column validation for Mode D CSV products."""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from common import strip_fasta

METRICS = ["reg_eq_similarity", "jaccard_similarity_ie", "containment_AB",
           "containment_BA", "cosketch_geom", "cosketch_arith", "cosketch_max"]
REQUIRED = {"file1", "file2", "sketch_type", "mode", "precision",
            "kmer_size", "window_size", *METRICS}


def validate_hll_csv(path: Path, samples: list[str], *, k: int, w: int,
                     precision: int) -> dict[str, Any]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks named columns {sorted(missing)}")
        rows = list(reader)
    expected_pairs = len(samples) ** 2
    if len(rows) != expected_pairs:
        raise ValueError(f"{path} has {len(rows)} rows, expected {expected_pairs}")
    expected = set(samples)
    seen: dict[tuple[str, str], dict[str, float]] = {}
    excursion_count = 0
    for row in rows:
        left, right = strip_fasta(row["file1"]), strip_fasta(row["file2"])
        if left not in expected or right not in expected:
            raise ValueError(f"unexpected sample pair {left}, {right}")
        if (left, right) in seen:
            raise ValueError(f"duplicate ordered pair {left}, {right}")
        if row["mode"] != "D" or row["sketch_type"] != "minimizer":
            raise ValueError(f"wrong mode/sketch metadata in {path}")
        if int(row["precision"]) != precision or int(row["kmer_size"]) != k or int(row["window_size"]) != w:
            raise ValueError(f"parameter metadata mismatch in {path}")
        values = {}
        for metric in METRICS:
            value = float(row[metric])
            if not math.isfinite(value):
                raise ValueError(f"non-finite {metric} for {left}, {right}")
            if value < -1e-12 or value > 1 + 1e-12:
                raise ValueError(f"out-of-range {metric}={value} for {left}, {right}")
            if value < 0 or value > 1:
                excursion_count += 1
                value = min(1.0, max(0.0, value))
            values[metric] = value
        seen[(left, right)] = values
    if set(left for left, _ in seen) != expected or set(right for _, right in seen) != expected:
        raise ValueError(f"sample axes do not match input inventory in {path}")
    for left in samples:
        for right in samples:
            if (left, right) not in seen:
                raise ValueError(f"missing ordered pair {left}, {right}")
            for metric in METRICS:
                value = seen[(left, right)][metric]
                reverse = seen[(right, left)][metric]
                # Directional containments transpose into one another.
                reverse_metric = ({"containment_AB": "containment_BA",
                                   "containment_BA": "containment_AB"}.get(metric, metric))
                reverse = seen[(right, left)][reverse_metric]
                if not math.isclose(value, reverse, rel_tol=1e-9, abs_tol=1e-9):
                    raise ValueError(f"asymmetry in {metric} for {left}, {right}")
        if not math.isclose(seen[(left, left)]["jaccard_similarity_ie"], 1.0,
                            rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"non-unit IE diagonal for {left}")
    return {"rows": len(rows), "samples": len(samples), "clamps": excursion_count,
            "columns": reader.fieldnames}
