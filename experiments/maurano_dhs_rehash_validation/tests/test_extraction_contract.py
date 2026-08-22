from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from extract_exact_features import extract, window_minimizer


def test_fixture_matches_production_selector_sets(tmp_path):
    fasta = tmp_path / "fixture.fa"
    sequences = ["ACGTACGTACGTACGTACGT", "TTTTCCCCAAAAGGGGTTTT", "ACG"]
    fasta.write_text("".join(f">s{i}\n{sequence}\n" for i, sequence in enumerate(sequences)))
    selectors, fallbacks, stats = extract(fasta, k=4, w=5)
    expected_selectors = set()
    expected_fallbacks = set()
    for sequence in sequences:
        selected = window_minimizer(sequence, k=4, w=5, include_hash=True)
        if selected:
            expected_selectors.update(int(value) for _, value in selected)
        else:
            expected_fallbacks.add(sequence.upper().encode())
    assert selectors == expected_selectors
    assert fallbacks == expected_fallbacks
    assert stats["records"] == 3


def test_reverse_complement_fixture_selector_identity(tmp_path):
    sequence = "ACGTTGCAACGTTGCAACGT"
    complement = str.maketrans("ACGT", "TGCA")
    reverse_complement = sequence.translate(complement)[::-1]
    left = tmp_path / "left.fa"
    right = tmp_path / "right.fa"
    left.write_text(f">left\n{sequence}\n")
    right.write_text(f">right\n{reverse_complement}\n")
    left_selectors, _, _ = extract(left, 5, 6)
    right_selectors, _, _ = extract(right, 5, 6)
    assert left_selectors == right_selectors
