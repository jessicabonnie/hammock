from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from exact_features import exact_jaccard, read_features, write_features


def test_exact_encoding_is_byte_deterministic(tmp_path):
    first = tmp_path / "first.bin"
    second = tmp_path / "second.bin"
    left = write_features(first, {7, 1, 2**32 - 1}, {b"ACGT", b"NN"})
    right = write_features(second, {2**32 - 1, 7, 1}, {b"NN", b"ACGT"})
    assert first.read_bytes() == second.read_bytes()
    assert left["sha256"] == right["sha256"]
    selectors, fallbacks = read_features(first)
    assert selectors.tolist() == [1, 7, 2**32 - 1]
    assert fallbacks == {b"ACGT", b"NN"}


def test_jaccard_keeps_selector_and_fallback_domains_separate():
    left = (np.array([1, 2, 5], dtype=np.uint32), {b"ACGT", b"1"})
    right = (np.array([2, 5, 9], dtype=np.uint32), {b"ACGT", b"2"})
    similarity, intersection, union = exact_jaccard(left, right)
    assert intersection == 3
    assert union == 7
    assert similarity == 3 / 7
