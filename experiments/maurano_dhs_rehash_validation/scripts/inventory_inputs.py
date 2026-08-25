#!/usr/bin/env python3
"""Freeze the 37-cell grid and checksummed, read-only biological inputs."""
from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

from common import EXPERIMENT, atomic_json, fasta_paths, grid, load_config, resolve, sha256, strip_fasta, write_csv


def historical_cells(path: Path) -> set[tuple[int, int]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"k", "w"}.issubset(reader.fieldnames or []):
            raise ValueError("historical grid lacks named k/w columns")
        rows = [(int(row["k"]), int(row["w"])) for row in reader]
    if len(rows) != len(set(rows)):
        raise ValueError("historical grid has duplicate cells")
    return set(rows)


def git(args: list[str], cwd: Path) -> str:
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=EXPERIMENT / "config.yaml")
    parser.add_argument("--skip-fasta-hashes", action="store_true",
                        help="Development-only: record size but not SHA-256 for large FASTAs")
    args = parser.parse_args()
    config, base = load_config(args.config)
    cells = grid(config)
    historical = historical_cells(resolve(base, config["inputs"]["historical_grid"]))
    if set(cells) != historical:
        missing = sorted(historical - set(cells))
        extra = sorted(set(cells) - historical)
        raise SystemExit(f"37-cell gate failed: missing={missing}, extra={extra}")

    fastas = fasta_paths(config, base)
    expected = int(config["samples_expected"])
    if len(fastas) != expected or len({strip_fasta(str(p)) for p in fastas}) != expected:
        raise SystemExit(f"expected {expected} unique FASTAs, found {len(fastas)}")
    missing_paths = [str(path) for path in fastas if not path.is_file()]
    if missing_paths:
        raise SystemExit(f"missing FASTAs: {missing_paths}")

    key = resolve(base, config["inputs"]["filename_key"])
    bedtools = resolve(base, config["inputs"]["bedtools_reference"])
    list_path = resolve(base, config["inputs"]["fasta_list"])
    rows = []
    for role, paths in (("fasta", fastas), ("filename_key", [key]),
                        ("bedtools_reference", [bedtools]), ("fasta_list", [list_path])):
        for path in paths:
            rows.append({
                "role": role, "sample": strip_fasta(str(path)) if role == "fasta" else "",
                "resolved_path": str(path), "bytes": path.stat().st_size,
                "sha256": "SKIPPED" if args.skip_fasta_hashes and role == "fasta" else sha256(path),
            })
    manifest_dir = EXPERIMENT / "results" / "manifests"
    write_csv(manifest_dir / "inputs.csv",
              ["role", "sample", "resolved_path", "bytes", "sha256"], rows)
    write_csv(manifest_dir / "grid.csv", ["cell_id", "k", "w"],
              ({"cell_id": f"k{k}_w{w}", "k": k, "w": w} for k, w in cells))

    repo = EXPERIMENT.parents[1]
    packages = {}
    for name in ("numpy", "pandas", "scipy", "scikit-learn", "biopython", "xxhash", "pyyaml"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    provenance = {
        "argv": sys.argv, "cwd": os.getcwd(), "python": sys.version,
        "executable": sys.executable, "platform": platform.platform(),
        "hostname": platform.node(), "cpu_count": os.cpu_count(),
        "repository_commit": git(["rev-parse", "HEAD"], repo),
        "repository_status": git(["status", "--short"], repo),
        "config_sha256": sha256(args.config.resolve()), "packages": packages,
        "hash_scheme": config["hash_scheme"], "seeds": config["seeds"],
        "precision": config["primary_precision"], "threads": config["threads"],
    }
    atomic_json(manifest_dir / "provenance.json", provenance)
    print(json.dumps({"cells": len(cells), "fastas": len(fastas),
                      "manifest": str(manifest_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
