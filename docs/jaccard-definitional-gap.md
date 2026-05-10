# The hammock vs bedtools Jaccard definitional gap

## Summary

`hammock-cpp` Mode B and `bedtools jaccard` report numerically different
Jaccard values for the same input, even at HLL precision high enough that
HLL noise is negligible. The two tools are estimating **different
quantities by design** — bedtools computes set Jaccard on base-pair
positions; hammock reports a *register-equality* statistic on the HLL
sketch, which is biased relative to the true set Jaccard.

This document exists because at first glance the disagreement looks like
a bug. It isn't, but it changes how we should set up and read precision
sweeps and other parity-style comparisons.

## The two estimators

**bedtools `jaccard`** — exact set-Jaccard on base-pair coverage:

```
J_bp = |A ∩ B in bp| / |A ∪ B in bp|
```

`A` and `B` are the sets of bp positions covered by intervals in each
file. Overlapping intervals within a file are merged before counting.

**hammock-cpp `--mode B`** — register-equality Jaccard on HyperLogLog
sketches built from those same bp positions
(`cpp/src/hll_sketch.cpp:40-60`):

```
active   = { i : R_A[i] != 0  ∨  R_B[i] != 0 }
matching = { i ∈ active : R_A[i] == R_B[i] }
J_re     = |matching| / |active|
```

This is the same formula the original Python `hammock` uses
(`estimate_jaccard_registers`). It does *not* go through inclusion-exclusion
on cardinality estimates — it compares register values directly.

## Why they disagree

Each HLL register stores `max ρ(h(x))` over elements `x` hashed into that
register, where `ρ` is the leading-zero count of the hash. Two registers
having the same value does not require them to have come from the same
element — they only need to have seen elements with the same leading-zero
count. With ρ geometrically distributed, "tie by chance" has a non-trivial
floor (~⅓ per active register pair).

For sets with low set-Jaccard but substantial union (the realistic case
for randomly-distributed BED intervals across the genome), most registers
are populated by both A and B but with values determined by the
*different* elements that happened to be argmax in each sketch. The
fraction of registers that nonetheless tie creates a positive bias floor
that no amount of HLL precision will erase.

The extreme cases agree:

| Input             | Set Jaccard | Register-equality Jaccard |
| ----------------- | ----------- | ------------------------- |
| A == B (identical sketches) | 1.0 | 1.0 (every active register matches by construction) |
| A ∩ B = ∅, |A∪B| ≪ m | 0.0 | 0.0 (no active register pairs to tie) |
| A ∩ B = ∅, |A∪B| ≫ m | 0.0 | > 0 (chance ties dominate) |
| 0 < J < 1, large union | J | f(J) > J |

The mapping `f(J)` is monotonic (more shared elements → more shared
argmaxes → more matching registers), but it is not the identity.

## Empirical evidence

From the precision smoke test in `experiments/bedtools_benchmark/` on two
1000-interval random BED files (seed 7, default 24-chrom × 10 Mb genome):

| Source | Jaccard |
| --- | --- |
| Exact bp-set computation | 0.010603 |
| `bedtools jaccard` | 0.010623 |
| hammock p=10 | 0.156 |
| hammock p=14 | 0.175 |
| hammock p=18 | 0.179 |

Hammock is converging to **~0.18**, not 0.01. The convergence pattern
(monotonic in p, plateauing as 1/√(2^p) noise vanishes) is exactly what
you'd expect from an unbiased estimator of register-equality Jaccard —
which is **not** set Jaccard.

## Where this matters in experiments

- **`experiments/bedtools_benchmark/sweep.py --axis precision`** —
  the precision sweep cannot use bedtools as ground truth for an HLL
  accuracy curve, because the systematic gap is much larger than HLL
  noise. We instead use hammock at maximum precision (e.g. p=18 or
  p=24) as the HLL ground truth, and *also* record bedtools Jaccards
  alongside so the gap is visible in the data and in plots.

- **Any "is hammock close to bedtools?" parity check** — only valid
  on inputs where set Jaccard is high. The byte-equal CSV parity tests
  (`tests/test_parity_against_original.py`) compare hammock-claude
  against the original hammock, which uses the **same** register-equality
  estimator, so those aren't affected.

- **Mode D FASTA work** — the same caveat applies; Mode D uses HLL on
  k-mers and reports the same register-equality Jaccard. Comparisons
  to `dashing` / `mash` / `sourmash` need to specify which estimator
  they're showing.

## How to read hammock vs bedtools plots

- Do **not** apply a constant offset (e.g. `J_hammock − 0.17`). The gap
  is not constant: it goes to 0 at J=1 and depends on |A∪B|/m at low J.
- The honest visualization is a **scatter of (J_bedtools, J_hammock)**
  across many pairs, with the y=x diagonal as a reference and a
  smoothed curve showing the empirical mapping. The shape of that curve
  is itself a finding — it tells you under what regimes hammock-Mode-B
  Jaccard does or doesn't track set Jaccard usefully.
- For HLL precision-accuracy claims (e.g. "p=14 is close enough"),
  compare hammock@p to hammock@p_max, not hammock@p to bedtools.

## See also

- `docs/bedtools-parallelism-caveat.md` — companion doc on the
  bedtools-via-GNU-parallel framing for hammock-vs-bedtools timing
  comparisons.

## References

- Implementation: `cpp/src/hll_sketch.cpp:40-60` (Jaccard formula)
- Algorithm parity note: `CLAUDE.md` ("register-equality Jaccard")
- Original Python source: `hammock` (orig) `intervals.py`,
  `estimate_jaccard_registers`
- Ertl, O. (2017). *New cardinality estimation algorithms for
  HyperLogLog sketches.* arXiv:1702.01284 — used for the cardinality
  estimator (`ertl_improved_estimate`), but the Jaccard formula
  hammock uses is the simpler register-equality variant inherited
  from the original Python tool, not Ertl's joint MLE.
