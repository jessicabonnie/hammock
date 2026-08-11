"""Coverage for the `--threads` sketch-phase-bounding fix
(docs/seed-hammock-cpp-file-dispatch.md Part 1).

Before this fix, `_core.sketch_bed_file_hll` had no thread parameter at all,
so `process_bed_file_mode_b/c`'s OMP team always ran at the ambient default
regardless of `--threads`. Actually asserting the OMP team size stayed within
budget isn't portably unit-testable (it depends on the runtime's core count
and OMP configuration, not just this repo's code) -- see the seed doc's own
empirical isolation script for that check
(`experiments/bedtools_benchmark/measure_threads_isolation.py`). What *is*
testable here, cheaply and portably:

1. The pure arithmetic that splits a `--threads` budget between the outer
   file-dispatch pool and the inner OMP team (`runner._split_thread_budget`).
2. That the new `threads=` argument on `_core.sketch_bed_file_hll` is
   accepted, doesn't crash, and -- since thread count must only affect
   parallelism, never the result -- produces the identical cardinality
   regardless of what's passed.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hammock.runner import _split_thread_budget

_core = pytest.importorskip("hammock._core")

DATA = Path(__file__).parent / "data"
_TINY_A = str(DATA / "tiny_a.bed")


# --- 1. Pure arithmetic ------------------------------------------------

@pytest.mark.parametrize(
    "total_threads, n, expected",
    [
        # n<=1: no file-dispatch parallelism possible, so the whole budget
        # goes to the inner OMP team -- this is the case the seed doc's own
        # single-file isolation measurement exercises.
        (1, 1, (1, 1)),
        (8, 1, (1, 8)),
        (16, 0, (1, 16)),
        # total_threads<=1: always serial outer dispatch.
        (1, 5, (1, 1)),
        # Evenly-divisible splits.
        (8, 4, (4, 2)),
        (8, 8, (8, 1)),
        (3, 3, (3, 1)),
        # n < total_threads: outer pool is capped at n, leftover budget goes
        # into the inner OMP team per file.
        (16, 4, (4, 4)),
        # Not evenly divisible: floor division under-uses the budget rather
        # than exceeding it (8 threads / 3 workers -> 2 each, product 6, not
        # 8) -- documented trade-off, not a bug.
        (8, 3, (3, 2)),
        (5, 2, (2, 2)),
        # n >> total_threads: outer pool caps at total_threads, inner team
        # floors at 1 (never 0 -- num_threads(0) is an invalid OMP clause).
        (4, 100, (4, 1)),
    ],
)
def test_split_thread_budget(total_threads, n, expected):
    assert _split_thread_budget(total_threads, n) == expected


@pytest.mark.parametrize("total_threads", [1, 2, 3, 5, 8, 16, 64])
@pytest.mark.parametrize("n", [0, 1, 2, 3, 5, 32, 1000])
def test_split_thread_budget_never_exceeds_budget(total_threads, n):
    outer_workers, sketch_threads = _split_thread_budget(total_threads, n)
    assert outer_workers >= 1
    assert sketch_threads >= 1
    assert outer_workers * sketch_threads <= total_threads


# --- 2. Binding accepts `threads=` and it doesn't change the result ----

@pytest.mark.parametrize("mode", ["A", "B", "C"])
@pytest.mark.parametrize("threads", [0, 1, 2, 4, 1024])
def test_sketch_bed_file_hll_accepts_threads_arg(mode, threads):
    """`threads=1024` exercises the 4x-num-procs cap branch in
    `omp_team_size` (cpp/include/hammock/omp_util.hpp) without needing to
    know the actual core count."""
    sketch = _core.sketch_bed_file_hll(
        path=_TINY_A,
        mode=mode,
        precision=12,
        threads=threads,
    )
    assert isinstance(sketch, _core.HLLSketch)
    assert sketch.estimate_cardinality() > 0


@pytest.mark.parametrize("mode", ["B", "C"])  # the modes with a parallel region
def test_threads_arg_does_not_change_cardinality(mode):
    """Thread count controls parallelism only -- verifies the fix didn't
    accidentally make `threads` a semantic knob (e.g. by threading it into
    something other than the `num_threads()` clause)."""
    cards = [
        _core.sketch_bed_file_hll(
            path=_TINY_A, mode=mode, precision=14, threads=t,
        ).estimate_cardinality()
        for t in (0, 1, 2, 8)
    ]
    assert len(set(cards)) == 1, f"cardinality varied with threads: {cards}"
