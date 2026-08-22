#!/usr/bin/env python3
"""Deterministic exact-feature encoding shared by extraction and tests."""
from __future__ import annotations

import hashlib
import json
import os
import struct
from pathlib import Path

import numpy as np

MAGIC = b"HMXSEL1\0"


def write_features(path: Path, selectors: set[int], fallbacks: set[bytes]) -> dict[str, object]:
    """Atomically write exact, domain-separated identities in stable order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    selector_array = np.asarray(sorted(selectors), dtype="<u4")
    fallback_values = sorted(fallbacks)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    digest = hashlib.sha256()
    with temporary.open("wb") as handle:
        def emit(value: bytes) -> None:
            handle.write(value)
            digest.update(value)
        emit(MAGIC)
        emit(struct.pack("<QQ", len(selector_array), len(fallback_values)))
        emit(selector_array.tobytes())
        for value in fallback_values:
            emit(struct.pack("<Q", len(value)))
            emit(value)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    return {
        "selector_distinct": len(selector_array), "fallback_distinct": len(fallback_values),
        "feature_distinct": len(selector_array) + len(fallback_values),
        "bytes": path.stat().st_size, "sha256": digest.hexdigest(),
    }


def read_features(path: Path) -> tuple[np.ndarray, set[bytes]]:
    with path.open("rb") as handle:
        if handle.read(len(MAGIC)) != MAGIC:
            raise ValueError(f"bad exact-feature magic: {path}")
        selector_count, fallback_count = struct.unpack("<QQ", handle.read(16))
        raw = handle.read(selector_count * 4)
        if len(raw) != selector_count * 4:
            raise ValueError(f"truncated selector payload: {path}")
        selectors = np.frombuffer(raw, dtype="<u4").copy()
        fallbacks = set()
        for _ in range(fallback_count):
            raw_length = handle.read(8)
            if len(raw_length) != 8:
                raise ValueError(f"truncated fallback length: {path}")
            length = struct.unpack("<Q", raw_length)[0]
            value = handle.read(length)
            if len(value) != length:
                raise ValueError(f"truncated fallback value: {path}")
            fallbacks.add(value)
        if handle.read(1):
            raise ValueError(f"trailing bytes in exact-feature file: {path}")
    if len(selectors) and (np.any(selectors[1:] <= selectors[:-1])):
        raise ValueError(f"selectors are not strictly sorted: {path}")
    return selectors, fallbacks


def exact_jaccard(left: tuple[np.ndarray, set[bytes]],
                  right: tuple[np.ndarray, set[bytes]]) -> tuple[float, int, int]:
    left_selectors, left_fallbacks = left
    right_selectors, right_fallbacks = right
    selector_intersection = int(np.intersect1d(
        left_selectors, right_selectors, assume_unique=True).size)
    intersection = selector_intersection + len(left_fallbacks & right_fallbacks)
    union = (len(left_selectors) + len(left_fallbacks) +
             len(right_selectors) + len(right_fallbacks) - intersection)
    similarity = 1.0 if union == 0 else intersection / union
    return similarity, intersection, union


def write_metadata(path: Path, value: dict[str, object]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
