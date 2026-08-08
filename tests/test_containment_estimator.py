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

    # Cardinalities must be ASYMMETRIC. With |A| == |B| the two containments
    # coincide and a whole family of wrong formulas (geometric mean of the
    # containments, min, harmonic mean, ...) agrees with I-E exactly, so a
    # symmetric fixture pins almost nothing.
    a = [range_sketch(14, 0, 40000), range_sketch(14, 0, 4000)]
    b = [range_sketch(14, 0, 4000), range_sketch(14, 2000, 42000)]

    _, c_ab, c_ba = _core.pairwise_metrics_hll(a, b)
    jac_ie = _jaccard_ie_from_containments(c_ab, c_ba)
    for i, sa in enumerate(a):
        for j, sb in enumerate(b):
            inter = sa.estimate_intersection(sb)
            union = sa.estimate_cardinality() + sb.estimate_cardinality() - inter
            expected = inter / union if union > 0 else 0.0
            assert jac_ie[i][j] == pytest.approx(expected, rel=1e-12)


def _hash_sketch(precision: int, n: int, seed: int) -> "_core.HLLSketch":
    """Sketch of `n` pseudo-random 64-bit hashes. Cheaper than range_sketch's
    per-key add_string, which makes the 1e7-element fixtures below tractable."""
    import random

    rnd = random.Random(seed)
    s = _core.HLLSketch(precision)
    for _ in range(n):
        s.add_hash64(rnd.getrandbits(64))
    return s


# The fused jaccard+union pass (HLLSketch::jaccard_and_union_cardinality) is
# claimed bit-identical to the route it replaced -- jaccard_similarity() plus
# union_with()->cardinality() reached via intersection_size(). That claim is an
# integer-multiset argument, so it deserves an exact test, not an approximate
# one: `estimate_intersection` and `estimate_cardinality` still expose the old
# scalar path untouched, which makes them an independent oracle.
#
# `==`, deliberately. `pytest.approx` here would pass through any drift up to
# its default 1e-6 relative tolerance and defeat the point. If this fails, the
# fused histogram has diverged from the materialized union -- fix that rather
# than loosening the comparison.
@pytest.mark.parametrize("precision", [4, 12, 18, 24])
def test_pairwise_metrics_exactly_matches_scalar_path(precision):
    a = [_hash_sketch(precision, 5000, 1), _hash_sketch(precision, 20000, 2)]
    b = [_hash_sketch(precision, 5000, 1), _hash_sketch(precision, 1000, 3)]

    jac, c_ab, c_ba = _core.pairwise_metrics_hll(a, b)
    for i, sa in enumerate(a):
        for j, sb in enumerate(b):
            inter = sa.estimate_intersection(sb)
            ca, cb = sa.estimate_cardinality(), sb.estimate_cardinality()
            assert jac[i][j] == sa.estimate_jaccard(sb)
            assert c_ab[i][j] == ((inter / ca) if ca > 0 else 0.0)
            assert c_ba[i][j] == ((inter / cb) if cb > 0 else 0.0)


def test_pairwise_metrics_exact_on_degenerate_sketches():
    """Empty (all-zero-register) sketches, self-pairs, and nested subsets.

    The all-zero case exercises cardinality()'s short-circuit, which the fused
    path has to reproduce off the union histogram rather than off the registers.
    """
    p = 14
    import random

    empty = _core.HLLSketch(p)
    small = _hash_sketch(p, 1000, 11)
    big = _hash_sketch(p, 50000, 12)
    # `nested` is a strict superset of `small`: it replays seed 11's first 1000
    # draws, then adds 9000 more from a different stream.
    nested = _hash_sketch(p, 1000, 11)
    rnd = random.Random(13)
    for _ in range(9000):
        nested.add_hash64(rnd.getrandbits(64))

    a = [empty, small, big, nested]
    b = [empty, small, nested]
    jac, c_ab, c_ba = _core.pairwise_metrics_hll(a, b)
    for i, sa in enumerate(a):
        for j, sb in enumerate(b):
            inter = sa.estimate_intersection(sb)
            ca, cb = sa.estimate_cardinality(), sb.estimate_cardinality()
            assert jac[i][j] == sa.estimate_jaccard(sb)
            assert c_ab[i][j] == ((inter / ca) if ca > 0 else 0.0)
            assert c_ba[i][j] == ((inter / cb) if cb > 0 else 0.0)


def test_pairwise_metrics_exact_at_extreme_size_ratio():
    """|A| ~ 1e3 against |B| ~ 1e7, A a subset of B.

    The regime where `inter` inherits ulp(cB) and the containment can exceed 1
    by ~3.4e4 ulp (see runner._jaccard_ie_from_containments). Exact agreement
    with the scalar path must hold there too.
    """
    p = 18
    tiny = _hash_sketch(p, 1000, 21)
    huge = _hash_sketch(p, 10_000_000, 22)
    a, b = [tiny], [huge]
    jac, c_ab, c_ba = _core.pairwise_metrics_hll(a, b)
    inter = tiny.estimate_intersection(huge)
    ca, cb = tiny.estimate_cardinality(), huge.estimate_cardinality()
    assert jac[0][0] == tiny.estimate_jaccard(huge)
    assert c_ab[0][0] == ((inter / ca) if ca > 0 else 0.0)
    assert c_ba[0][0] == ((inter / cb) if cb > 0 else 0.0)


# The pairwise loops used to ignore --threads entirely: omp_set_num_threads was
# called only by the standalone C++ binary, so `hammock --threads 4` inside a
# 4-CPU cgroup still ran this loop with a team per core on the whole node. The
# team size must not change any emitted value.
@pytest.mark.parametrize("threads", [0, 1, 2, 7])
def test_pairwise_metrics_values_are_thread_count_invariant(threads):
    a = [_hash_sketch(14, 20000, 31), _hash_sketch(14, 5000, 32)]
    b = [_hash_sketch(14, 20000, 31), _hash_sketch(14, 40000, 33)]

    ref = _core.pairwise_metrics_hll(a, b, threads=0)
    got = _core.pairwise_metrics_hll(a, b, threads=threads)
    for r, g in zip(ref, got):
        for i in range(len(a)):
            for j in range(len(b)):
                assert r[i][j] == g[i][j]


def test_pairwise_jaccard_accepts_threads_too():
    a = [_hash_sketch(12, 5000, 41)]
    b = [_hash_sketch(12, 5000, 42)]
    assert _core.pairwise_jaccard_hll(a, b, threads=3)[0][0] == \
        _core.pairwise_jaccard_hll(a, b)[0][0]
