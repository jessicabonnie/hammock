#!/usr/bin/env python3
"""Shared contracts for the Maurano independent-rehash experiment."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable

import yaml

EXPERIMENT = Path(__file__).resolve().parents[1]


def load_config(path: Path) -> tuple[dict[str, Any], Path]:
    path = path.resolve()
    with path.open() as handle:
        config = yaml.safe_load(handle)
    return config, path.parent


def resolve(base: Path, value: str) -> Path:
    return (base / value).resolve()


def grid(config: dict[str, Any]) -> list[tuple[int, int]]:
    cells = [(int(k), int(w)) for k in config["k_values"]
             for w in config["w_values"] if int(w) >= int(k)]
    if len(cells) != 37 or len(cells) != len(set(cells)):
        raise ValueError(f"frozen grid must contain 37 distinct cells, got {len(cells)}")
    return cells


def strip_fasta(value: str) -> str:
    name = Path(value).name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".bed.gz",
                   ".fasta", ".fna", ".fa", ".bed"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return name


def fasta_paths(config: dict[str, Any], base: Path) -> list[Path]:
    list_path = resolve(base, config["inputs"]["fasta_list"])
    paths = []
    for raw in list_path.read_text().splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        candidate = Path(raw)
        if not candidate.is_absolute():
            candidate = (list_path.parent / candidate).resolve()
        paths.append(candidate)
    return paths


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)
