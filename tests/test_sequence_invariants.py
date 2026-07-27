"""Mode D sketch invariants that matter for BED→FASTA / cross-reference use.

These lock behaviours the bed2fasta feature relies on:
  * case-invariance — UCSC references are soft-masked (lowercase repeats), so
    a soft-masked reference must sketch identically to an unmasked one. The
    guarantee comes from `digest.window_minimizer` being case-insensitive; this
    test fails loudly if that ever changes.
  * N-run behaviour — assembly-gap 'N' runs are ingested (not skipped), which
    is why bed2fasta warns on high-N extractions.
"""
from __future__ import annotations

import pytest

from hammock.modes.sequence import MinimizerSketch, _DIGEST_AVAILABLE

digest_only = pytest.mark.skipif(not _DIGEST_AVAILABLE, reason="digest not installed")

# A sequence long enough (> window) to exercise the real minimizer path.
_SEQ = "ACGTTGCAACGTACGTTTGGCCAAACGTAGCTAGCTAGCTACGATCGATCGGGCCTTAA" * 4


@digest_only
def test_minimizer_sketch_is_case_insensitive() -> None:
    up = MinimizerSketch(kmer_size=8, window_size=40, precision=12)
    lo = MinimizerSketch(kmer_size=8, window_size=40, precision=12)
    up.add_string(_SEQ)
    lo.add_string(_SEQ.lower())
    # Soft-masked (lowercase) input must produce an identical sketch.
    assert up.minimizer_hll.estimate_jaccard(lo.minimizer_hll) == pytest.approx(1.0)
    assert up.merged().estimate_jaccard(lo.merged()) == pytest.approx(1.0)


@digest_only
def test_n_runs_are_ingested_not_skipped() -> None:
    """N is hashed, not dropped: an all-N run (longer than the window) still
    produces minimizers. That is why gap-heavy extractions carry a shared
    'N minimizer' and bed2fasta warns on high-N output. (The downstream
    similarity impact isn't asserted here — the HLL estimator is unreliable at
    the tiny cardinalities of these toy sketches.)"""
    alln = MinimizerSketch(kmer_size=8, window_size=40, precision=12)
    alln.add_string("N" * 200)
    assert alln.minimizer_hll.estimate_cardinality() > 0
