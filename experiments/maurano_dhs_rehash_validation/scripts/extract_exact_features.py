#!/usr/bin/env python3
"""Extract exact production selector identities and build Jaccard matrices."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path

import numpy as np
from Bio import SeqIO

from common import EXPERIMENT, atomic_json, fasta_paths, grid, load_config, strip_fasta, write_csv
from exact_features import exact_jaccard, read_features, write_features

try:
    from hammock.modes.sequence import window_minimizer
except ImportError as error:  # pragma: no cover - environment gate
    raise SystemExit(f"isolated hammock build is required: {error}")


def records(path: Path):
    if path.name.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            yield from SeqIO.parse(handle, "fasta")
    else:
        yield from SeqIO.parse(path, "fasta")


def extract(path: Path, k: int, w: int) -> tuple[set[int], set[bytes], dict[str, int]]:
    selectors: set[int] = set()
    fallbacks: set[bytes] = set()
    stats = {"records": 0, "selector_occurrences": 0, "fallback_occurrences": 0,
             "selector_failures": 0}
    for record in records(path):
        stats["records"] += 1
        sequence = str(record.seq)
        if not sequence:
            continue
        try:
            selected = window_minimizer(sequence, k=k, w=w, include_hash=True)
        except Exception:
            selected = []
            stats["selector_failures"] += 1
        if not selected:
            stats["fallback_occurrences"] += 1
            fallbacks.add(sequence.upper().encode())
            continue
        stats["selector_occurrences"] += len(selected)
        for _, raw_value in selected:
            value = int(raw_value)
            if not 0 <= value < 2**32:
                raise ValueError(f"digest selector is not uint32: {value} in {path}")
            selectors.add(value)
    return selectors, fallbacks, stats


def feature_path(cell_dir: Path, sample: str) -> Path:
    return cell_dir / "features" / f"{sample}.features.bin"


def build_matrix(cell_dir: Path, samples: list[str], k: int, w: int) -> None:
    features = {sample: read_features(feature_path(cell_dir, sample)) for sample in samples}
    rows = []
    matrix = np.empty((len(samples), len(samples)), dtype=float)
    for i, left in enumerate(samples):
        for j, right in enumerate(samples):
            similarity, intersection, union = exact_jaccard(features[left], features[right])
            matrix[i, j] = similarity
            rows.append({"file1": left, "file2": right, "k": k, "w": w,
                         "jaccard_exact": f"{similarity:.17g}",
                         "intersection_exact": intersection, "union_exact": union})
    if not (np.isfinite(matrix).all() and np.all((matrix >= 0) & (matrix <= 1))
            and np.array_equal(matrix, matrix.T) and np.array_equal(np.diag(matrix), np.ones(len(samples)))):
        raise ValueError(f"exact matrix invariants failed for k={k}, w={w}")
    write_csv(cell_dir / "pairs.csv",
              ["file1", "file2", "k", "w", "jaccard_exact", "intersection_exact", "union_exact"], rows)
    matrix_path = cell_dir / "matrix.csv"
    with matrix_path.with_name(matrix_path.name + f".tmp.{os.getpid()}").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sample", *samples])
        for sample, values in zip(samples, matrix):
            writer.writerow([sample, *(f"{value:.17g}" for value in values)])
        temporary = Path(handle.name)
    os.replace(temporary, matrix_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    parser.add_argument("--cell", nargs=2, type=int, metavar=("K", "W"))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config, base = load_config(args.config)
    cells = grid(config)
    if args.cell:
        requested = tuple(args.cell)
        if requested not in cells:
            raise SystemExit(f"cell {requested} is not in frozen grid")
        cells = [requested]
    fastas = fasta_paths(config, base)
    samples = [strip_fasta(str(path)) for path in fastas]
    plan = [{"cell_id": f"k{k}_w{w}", "k": k, "w": w, "sample": sample,
             "fasta": str(path), "output": str(feature_path(
                 EXPERIMENT / "results" / "exact" / f"k{k}_w{w}", sample))}
            for k, w in cells for sample, path in zip(samples, fastas)]
    write_csv(EXPERIMENT / "results" / "manifests" / "exact_jobs.csv",
              ["cell_id", "k", "w", "sample", "fasta", "output"], plan)
    if args.dry_run:
        print(f"planned {len(plan)} exact sample/cell extractions")
        return 0

    for k, w in cells:
        cell_dir = EXPERIMENT / "results" / "exact" / f"k{k}_w{w}"
        metadata_dir = cell_dir / "metadata"
        metadata_dir.mkdir(parents=True, exist_ok=True)
        for sample, fasta in zip(samples, fastas):
            output = feature_path(cell_dir, sample)
            metadata_path = metadata_dir / f"{sample}.json"
            if output.exists() and metadata_path.exists() and not args.force:
                continue
            selectors, fallbacks, stats = extract(fasta, k, w)
            artifact = write_features(output, selectors, fallbacks)
            metadata = {"sample": sample, "fasta": str(fasta), "k": k, "w": w,
                        **stats, **artifact,
                        "selector_duplicate_occurrences": stats["selector_occurrences"] - len(selectors),
                        "fallback_duplicate_occurrences": stats["fallback_occurrences"] - len(fallbacks)}
            atomic_json(metadata_path, metadata)
        build_matrix(cell_dir, samples, k, w)
        atomic_json(cell_dir / "complete.json", {"k": k, "w": w, "samples": samples,
                    "ordered_pairs": len(samples) ** 2, "status": "complete"})
        print(f"completed exact k={k}, w={w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
