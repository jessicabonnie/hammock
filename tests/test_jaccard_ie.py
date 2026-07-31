"""Tests for the inclusion-exclusion Jaccard column (`jaccard_similarity_ie`).

The point of the column is that it estimates a *different quantity* from
`jaccard_similarity`: register-equality carries a chance-agreement floor whose
size depends on sketch load and on |A|/|B|, while inclusion-exclusion does not.
The regression that matters is that the two columns stay distinguishable in the
regime where the floor bites.

Numeric constants here are regime-dependent (they are functions of the load
factor n/m, not properties of the code), so assertions are stated as wide-margin
*relations* rather than golden values.
"""
from __future__ import annotations

import csv
import random
import subprocess
import sys
from pathlib import Path

import pytest

_core = pytest.importorskip("hammock._core")
from hammock.runner import _jaccard_ie_from_containments  # noqa: E402


def _sketch(precision: int, hashes) -> "_core.HLLSketch":
    s = _core.HLLSketch(precision)
    for h in hashes:
        s.add_hash64(h)
    return s


def _ie(a, b) -> float:
    _, c_ab, c_ba = _core.pairwise_metrics_hll([a], [b])
    return float(_jaccard_ie_from_containments(c_ab, c_ba)[0, 0])


def _regeq(a, b) -> float:
    jac, _, _ = _core.pairwise_metrics_hll([a], [b])
    return float(jac[0, 0])


# --- the helper in isolation ------------------------------------------------

def test_scalar_inputs():
    # Scalars only reach the helper from tests -- pairwise_metrics_hll always
    # returns 2-D arrays -- but the helper must not blow up on them. Assert on
    # value, not type: the scalar path returns a 0-d ndarray, not a float.
    assert float(_jaccard_ie_from_containments(0.5, 0.5)) == pytest.approx(1 / 3)
    assert float(_jaccard_ie_from_containments(1.0, 1.0)) == 1.0
    assert float(_jaccard_ie_from_containments(1, 1)) == 1.0  # integer dtype


def test_zero_containment_scores_zero_not_undefined():
    # A zero containment means the intersection estimate was zero -- genuinely
    # empty, or clamped from a negative. Scored 0.0, never NaN/inf.
    for c_ab, c_ba in ((0.0, 0.0), (0.0, 0.5), (0.5, 0.0)):
        v = float(_jaccard_ie_from_containments(c_ab, c_ba))
        assert v == 0.0


def test_containment_above_one_is_clamped():
    # Ertl noise can push a containment just past 1.0; unclamped that would
    # drive the denominator non-positive.
    v = float(_jaccard_ie_from_containments(1.0000000000050957, 1.0))
    assert v == 1.0


def test_no_warnings_emitted():
    import warnings

    import numpy as np

    # Include subnormals and non-finite values: 1/1e-320 overflows, and the
    # gate has to keep that off the live branch too.
    a = np.array([[0.0, 0.5], [1.0, 0.0], [1e-320, 5e-324], [np.nan, np.inf]])
    b = np.array([[0.5, 0.0], [1.0, 0.0], [1e-320, 1.0], [1.0, -1.0]])
    with warnings.catch_warnings(), np.errstate(all="raise"):
        warnings.simplefilter("error")
        out = _jaccard_ie_from_containments(a, b)
    assert np.all(np.isfinite(out)), out
    assert np.all((out >= 0.0) & (out <= 1.0)), out


# --- against real sketches --------------------------------------------------

@pytest.mark.parametrize("precision", [4, 8, 12, 16, 20])
def test_self_pair_is_exactly_one(precision: int):
    # For identical registers the union estimate equals both cardinalities
    # bitwise, so inter = 2c - c = c exactly and 1/(1 + 1 - 1) == 1.0.
    # Exact equality, not approx.
    s = _sketch(precision, [(i * 2654435761) % (1 << 64) for i in range(20000)])
    assert _ie(s, s) == 1.0
    assert _regeq(s, s) == 1.0


def test_empty_sketches_score_zero():
    # Note this is an *exception* to test_self_pair_is_exactly_one: two empty
    # sketches are identical, but with no active registers both estimators
    # report 0.0 rather than 1.0.
    empty, other = _core.HLLSketch(12), _sketch(12, range(1000))
    assert _ie(empty, empty) == 0.0
    assert _regeq(empty, empty) == 0.0
    assert _ie(empty, other) == 0.0
    assert _ie(other, empty) == 0.0


def _known_jaccard_pair(precision: int, n_a: int, n_b: int, n_shared: int):
    """Two sketches over integer key ranges with an exactly known set Jaccard."""
    assert n_shared <= min(n_a, n_b)
    # Distinct universes per element index, hashed so the low bits (which pick
    # the register) are unstructured.
    def h(i: int) -> int:
        return random.Random(i).getrandbits(64)
    shared = [h(i) for i in range(n_shared)]
    a = _sketch(precision, shared + [h(1_000_000 + i) for i in range(n_a - n_shared)])
    b = _sketch(precision, shared + [h(5_000_000 + i) for i in range(n_b - n_shared)])
    true_j = n_shared / (n_a + n_b - n_shared)
    return a, b, true_j


def test_ie_tracks_true_jaccard_on_asymmetric_sets():
    """Pins the *magnitude* on a pair whose cardinalities differ 10x.

    This is the assertion that actually constrains the formula: at |A| != |B|
    the two containments differ, so I-E is separated from the family of wrong
    expressions (geometric/harmonic/min of the containments) that coincide with
    it whenever c_ab == c_ba.
    """
    a, b, true_j = _known_jaccard_pair(18, n_a=200_000, n_b=20_000, n_shared=20_000)
    assert true_j == pytest.approx(0.1)          # B is a strict subset of A
    j_ie, j_re = _ie(a, b), _regeq(a, b)
    assert j_ie == pytest.approx(true_j, abs=0.01)
    # Register-equality is biased upward, but note the bias is SMALL here
    # (measured 0.111 vs 0.100): the chance floor shrinks as |A|/|B| grows --
    # ~0.058 at a 10:1 ratio versus 0.1699 at 1:1 -- which is exactly why it
    # cannot be treated as a constant offset. Assert the ordering that holds
    # in every regime rather than a margin that does not.
    assert j_re > true_j
    assert abs(j_ie - true_j) < abs(j_re - true_j)


def test_ie_tracks_true_jaccard_at_moderate_overlap():
    a, b, true_j = _known_jaccard_pair(18, n_a=100_000, n_b=100_000, n_shared=10_000)
    assert true_j == pytest.approx(10_000 / 190_000)
    assert _ie(a, b) == pytest.approx(true_j, abs=0.01)


def test_disjoint_inputs_separate_the_two_estimators():
    """On disjoint inputs at high load, register-equality sits near its chance
    floor while I-E collapses toward zero.

    Note the I-E side is a weak assertion by construction -- it is satisfied by
    a formula that always returns 0 -- so it is here to characterize the floor,
    not to pin the formula. The magnitude tests above do that.
    """
    n = 200_000
    rng = random.Random(20260731)
    # Must be genuinely random 64-bit values: HLL takes the register index from
    # the *low* bits, so a multiplicative sequence would place the two sets in
    # structured, non-overlapping register patterns and defeat the test.
    a = _sketch(16, [rng.getrandbits(64) for _ in range(n)])
    b = _sketch(16, [rng.getrandbits(64) for _ in range(n)])

    j_re, j_ie = _regeq(a, b), _ie(a, b)
    # Theory puts the high-load equal-cardinality floor at 0.1699; measured
    # 0.166 here. Wide margins -- it is a function of n/m, not a code constant.
    assert j_re > 0.10, j_re
    assert j_ie < 0.02, j_ie


def test_high_overlap_estimators_converge():
    # The floor vanishes as J -> 1, so the two columns must agree there.
    rng = random.Random(11)
    shared = [rng.getrandbits(64) for _ in range(50000)]
    a = _sketch(16, shared)
    b = _sketch(16, shared + [rng.getrandbits(64) for _ in range(500)])
    assert _ie(a, b) == pytest.approx(_regeq(a, b), abs=0.05)


# --- precision validation ---------------------------------------------------

@pytest.mark.parametrize("precision", [3, 25, 64, 65])
def test_out_of_range_precision_rejected(precision: int):
    # -p 64 used to return silently wrong estimates and -p 65 used to hang on a
    # ~1.8e19-iteration loop, because the member initializers ran before the
    # constructor body's check.
    with pytest.raises((ValueError, RuntimeError)):
        _core.HLLSketch(precision=precision)


@pytest.mark.parametrize("precision", [4, 18, 24])
def test_in_range_precision_accepted(precision: int):
    assert _core.HLLSketch(precision=precision).precision() == precision


def test_cli_rejects_out_of_range_precision(tmp_path: Path):
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t0\t100\n")
    lst = tmp_path / "files.txt"
    lst.write_text(f"{bed}\n")
    # Must fail fast (argparse, exit 2), not hang.
    r = subprocess.run(
        [sys.executable, "-m", "hammock.cli", str(lst), str(lst),
         "-p", "65", "-o", str(tmp_path / "out")],
        capture_output=True, text=True, timeout=60)
    assert r.returncode != 0
    assert "4..24" in r.stderr


# --- end to end -------------------------------------------------------------

def test_csv_has_both_columns_and_they_agree_at_high_overlap(tmp_path: Path):
    """High J and high precision -- I-E is noisiest at low J, so a tolerance
    assertion belongs in the regime where it is well behaved."""
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("".join(f"chr1\t{i * 1000}\t{i * 1000 + 900}\n"
                         for i in range(500)))
    # ~95% base-pair overlap with a.
    b.write_text("".join(f"chr1\t{i * 1000 + 45}\t{i * 1000 + 900}\n"
                         for i in range(500)))
    lst_a, lst_b = tmp_path / "la.txt", tmp_path / "lb.txt"
    lst_a.write_text(f"{a}\n")
    lst_b.write_text(f"{b}\n")

    subprocess.run(
        [sys.executable, "-m", "hammock.cli", str(lst_a), str(lst_b),
         "--mode", "B", "-p", "20", "-o", str(tmp_path / "out")],
        check=True, capture_output=True, text=True, timeout=300)

    rows = list(csv.DictReader(next(tmp_path.glob("out*.csv")).read_text().splitlines()))
    assert len(rows) == 1
    row = rows[0]
    assert "jaccard_similarity_ie" in row
    j_re = float(row["jaccard_similarity"])
    j_ie = float(row["jaccard_similarity_ie"])
    assert 0.0 <= j_ie <= 1.0
    # True bp Jaccard is ~0.95; I-E should land on it, register-equality high.
    assert j_ie == pytest.approx(0.95, abs=0.03), (j_re, j_ie)
    assert j_re >= j_ie - 0.03
