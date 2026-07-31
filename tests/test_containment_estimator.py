"""Numeric tests for the containment / co-sketch block.

The pre-existing coverage only asserted column presence, `0 <= v <= 1`, and
self-pair `== 1.0`. The self-pair case is exact by construction — comparing a
sketch with itself makes the union registers bit-identical to the operand's, so
the two Ertl estimates cancel — which means it validates column wiring and not
the estimator behind it. These tests pin the values against known ground-truth
intersections instead.

Background on the two estimators in play (they are NOT on the same scale, and
mixing them is a real analysis error): docs/jaccard-definitional-gap.md.
"""
from __future__ import annotations

import math

import pytest

_core = pytest.importorskip("hammock._core")


def range_sketch(precision: int, lo: int, hi: int) -> "_core.HLLSketch":
    """Sketch over the half-open key range [lo, hi)."""
    s = _core.HLLSketch(precision)
    for i in range(lo, hi):
        s.add_string(f"key_{i}")
    return s


def containments(a, b):
    """(containment_AB, containment_BA) the same way pairwise_metrics_hll does."""
    inter = a.estimate_intersection(b)
    return inter / a.estimate_cardinality(), inter / b.estimate_cardinality()


def test_known_half_overlap():
    # |A| = |B| = 20000, |A n B| = 10000 -> both containments 0.5,
    # true set Jaccard 10000/30000.
    a = range_sketch(16, 0, 20000)
    b = range_sketch(16, 10000, 30000)

    assert a.estimate_intersection(b) == pytest.approx(10000, rel=0.10)
    c_ab, c_ba = containments(a, b)
    assert c_ab == pytest.approx(0.5, rel=0.10)
    assert c_ba == pytest.approx(0.5, rel=0.10)

    j_ie = 1.0 / (1.0 / c_ab + 1.0 / c_ba - 1.0)
    assert j_ie == pytest.approx(1 / 3, rel=0.15)


def test_containment_asymmetry_under_strict_containment():
    # A strictly contains B: containment_BA is exactly 1.0 (register dominance
    # makes the union identical to A), containment_AB ~ |B|/|A| = 0.25.
    a = range_sketch(16, 0, 40000)
    b = range_sketch(16, 0, 10000)

    c_ab, c_ba = containments(a, b)
    assert c_ba == pytest.approx(1.0, rel=0.02)
    assert c_ab == pytest.approx(0.25, rel=0.15)
    assert c_ab < c_ba


def test_disjoint_containment_near_zero_but_jaccard_column_is_not():
    """The scale mismatch that makes the two column families non-comparable.

    Inclusion-exclusion has no chance-agreement term; register-equality does,
    and its floor is set by the load factor n/m. On disjoint inputs the
    containments read ~0 while `jaccard_similarity` sits well above it.
    """
    a = range_sketch(14, 0, 200_000)
    b = range_sketch(14, 200_000, 400_000)

    c_ab, _ = containments(a, b)
    assert c_ab < 0.02
    assert a.estimate_jaccard(b) > 0.10


def test_intersection_never_negative():
    # Inclusion-exclusion is a difference of large estimates and goes negative
    # for disjoint inputs roughly half the time; the >= 0 clamp must hold.
    for trial in range(8):
        base = trial * 100_000
        a = range_sketch(12, base, base + 30_000)
        b = range_sketch(12, base + 50_000, base + 80_000)
        assert a.estimate_intersection(b) >= 0.0


def test_cosketch_identities():
    from hammock.runner import _cosketch_from_containments

    c_ab, c_ba = 0.4, 0.9
    geom, arith, mx = _cosketch_from_containments(c_ab, c_ba)
    assert geom == pytest.approx(math.sqrt(c_ab * c_ba))
    assert arith == pytest.approx((c_ab + c_ba) / 2)
    assert mx == pytest.approx(max(c_ab, c_ba))


def test_mismatched_precision_raises_not_aborts():
    """Regression: this used to be a heap overread, or a SIGABRT via OpenMP.

    `union_with` indexed the operand's registers over *this* sketch's register
    count, so p=18 vs p=8 overread (and, on the reversed call, corrupted) the
    heap. Through `pairwise_metrics_hll` the resulting throw escaped an OpenMP
    structured block with the GIL released, terminating the interpreter instead
    of surfacing as a Python exception.
    """
    big = range_sketch(18, 0, 100)
    small = range_sketch(8, 0, 100)

    with pytest.raises(RuntimeError):
        big.estimate_intersection(small)
    with pytest.raises(RuntimeError):
        small.estimate_intersection(big)
    with pytest.raises(RuntimeError):
        big.estimate_jaccard(small)


def test_pairwise_metrics_mismatched_precision_raises_not_aborts():
    # Must raise a normal Python exception; a regression here shows up as the
    # test process dying with SIGABRT (exit 134) rather than as a failure.
    a = [range_sketch(14, 0, 100), range_sketch(8, 0, 100)]
    b = [range_sketch(14, 0, 100)]

    with pytest.raises((ValueError, RuntimeError)):
        _core.pairwise_metrics_hll(a, b)
    with pytest.raises((ValueError, RuntimeError)):
        _core.pairwise_jaccard_hll(a, b)


def test_pairwise_metrics_matches_scalar_path():
    """The parallel binding must agree with the per-pair methods it replaces."""
    a = [range_sketch(14, 0, 5000), range_sketch(14, 2500, 7500)]
    b = [range_sketch(14, 0, 5000), range_sketch(14, 5000, 10000)]

    jac, c_ab, c_ba = _core.pairwise_metrics_hll(a, b)
    for i, sa in enumerate(a):
        for j, sb in enumerate(b):
            inter = sa.estimate_intersection(sb)
            assert jac[i][j] == pytest.approx(sa.estimate_jaccard(sb))
            assert c_ab[i][j] == pytest.approx(inter / sa.estimate_cardinality())
            assert c_ba[i][j] == pytest.approx(inter / sb.estimate_cardinality())


def test_jaccard_ie_matches_direct_inclusion_exclusion():
    """The shipped `jaccard_similarity_ie` column must equal I/(|A|+|B|-I).

    This is the test that pins the definitional relationship between the new
    column and the containment block next to it. A "close to zero on disjoint
    inputs" assertion would pass for a range of wrong formulas; this doesn't.
    """
    from hammock.runner import _jaccard_ie_from_containments

    a = [range_sketch(14, 0, 5000), range_sketch(14, 2500, 7500)]
    b = [range_sketch(14, 0, 5000), range_sketch(14, 5000, 10000)]

    _, c_ab, c_ba = _core.pairwise_metrics_hll(a, b)
    jac_ie = _jaccard_ie_from_containments(c_ab, c_ba)
    for i, sa in enumerate(a):
        for j, sb in enumerate(b):
            inter = sa.estimate_intersection(sb)
            union = sa.estimate_cardinality() + sb.estimate_cardinality() - inter
            expected = inter / union if union > 0 else 0.0
            assert jac_ie[i][j] == pytest.approx(expected, rel=1e-12)
