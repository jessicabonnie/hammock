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
fraction of registers that nonetheless tie creates a positive bias floor.

That floor is governed by the **load factor** `λ = n/m` (distinct elements
per register), not by precision as such. Raising `p` at fixed input does
eventually shrink it — but only once `m` overtakes `n`, at which point most
registers are empty or singly-active and there are few two-sided pairs left
to tie. Measured on ~5 Mbp inputs (see below), the fitted intercept `c` is
0.1800 / 0.1803 / 0.1803 at p = 12 / 16 / 20 — dead flat while `m ≪ n` —
then 0.0454 at p = 24, where `m` = 16.7 M exceeds `n` ≈ 5 M. So "more
precision won't fix it" is right in the regime the tool normally runs in,
and wrong as an absolute statement.

Analytically `c = Σ_k p_k²` where `p_k = e^(−λ2^(−k)) − e^(−λ2^(−(k−1)))` is
the probability a register lands on value `k`. At λ = 92 that sum is 0.169,
matching the measured 0.17. Note this is well below the `Σ_k (2^−k)² = 1/3`
figure quoted as the collision bound elsewhere — 1/3 is an upper bound that
is not approached at any realistic λ; `c` is quasi-periodic in `log₂ λ`,
sits around 0.15–0.18 for λ ≳ 20, and collapses toward 0 as λ → 0.

**`c` also depends on the cardinality ratio `|A|/|B|`, not just on λ.** For
same-size files at λ ≫ 1 it is ~0.17; at a 100:1 size ratio the same
simulation gives `J_re` ≈ 0.0142 against a true J of 0.010 — nearly
unbiased. Both effects mean `c` is **not** a constant you can calibrate away
with a single offset, and any quoted value is specific to a corpus, a
precision, *and* a size regime. The experiment below pins the ratio at ≈1 and
so does not probe this axis at all; asymmetric-size pairs are common in real
BED comparisons and remain untested.

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

## The containment columns are *not* affected — and give a way out

Only `jaccard_similarity` uses register equality. The `containment_AB` /
`containment_BA` columns (and the three `cosketch_*` derived from them) go
through `HLLSketch::intersection_size` (`cpp/src/hll_sketch.cpp:146`), which
uses **inclusion–exclusion** — `|A| + |B| − |A ∪ B|`, Ertl estimator on each,
union by register-wise max, clamped to `>= 0`. That path has no
chance-agreement term, so those columns estimate the true set quantities.
On disjoint inputs at p=16, n=2×10⁵: `jaccard_similarity` = 0.168 while
`containment_AB` = 0.000.

Consequently a set-Jaccard estimate is recoverable from output you already
have, with no rerun and no code change. From `C_AB = I/|A|` and
`C_BA = I/|B|`:

```
J_set = I/(|A| + |B| − I) = 1 / (1/C_AB + 1/C_BA − 1)
```

Measured against `bedtools jaccard` on 90 synthetic pairs spanning
J = 0.003 … 0.99 (matched row sets; undefined-I-E scored as 0):

| p | MAE `jaccard_similarity` | MAE reconstructed `J_set` | ratio | I-E undefined |
| --- | --- | --- | --- | --- |
| 12 | 0.15173 | 0.00822 | 18× | 25/90 |
| 16 | 0.15173 | 0.00207 | 73× | 0 |
| 20 | 0.15164 | 0.00054 | 283× | 0 |
| 24 | 0.03904 | 0.00011 | 361× | 0 |

In the low-J regime real BED comparisons occupy (J < 0.05, mean 0.0114),
register-equality reports **0.181** and the reconstruction reports
**0.0115**.

### Neither estimator dominates — MAE is not the whole story

The table above scores **calibration** (is the number right?). It does not
score **resolution** (can the number tell two similar pairs apart?), and on
resolution the ranking flips.

Because `J_re` is a near-perfectly affine function of true J, its error must
be rescaled by the fitted slope `1 − c` before its noise can be compared with
`J_ie`'s. Doing that, in true-J units:

| p | fitted c | J < 0.01 | 0.01–0.05 | 0.05–0.15 | winner at low J |
| --- | --- | --- | --- | --- | --- |
| 12 | 0.1800 | 0.00802 / 0.01023 | 0.00577 / 0.00924 | 0.00647 / 0.01335 | **register-equality** |
| 16 | 0.1803 | 0.00147 / 0.00242 | 0.00142 / 0.00235 | 0.00240 / 0.00491 | **register-equality** |
| 20 | 0.1803 | 0.00037 / 0.00061 | 0.00060 / 0.00059 | 0.00208 / 0.00060 | tie → I-E |
| 24 | 0.0454 | 0.00024 / 0.00012 | 0.00052 / 0.00012 | 0.00118 / 0.00009 | **inclusion–exclusion** |

(cells are `sd(err_re)/(1−c)` / `sd(err_ie)`)

So in the **saturated** regime — `m ≪ n`, which is where genome-scale BED
data actually sits — register-equality is the *more precise* discriminator at
low J even though it is far more biased, and inclusion–exclusion only takes
over once the sketch is over-provisioned relative to the data. Inclusion–
exclusion is also **censored**: it is clamped at 0 (`hll_sketch.cpp:157`),
which fires on 25/90 pairs at p=12, all at low J. Nothing distinguishes
"clamped from a large negative" from "genuinely zero", so conditional on
being non-zero it is biased *upward* near J=0 — a small positive floor of its
own (~0.002–0.014 on truly disjoint pairs, falling with p).

Rule of thumb: its absolute error is ≈ `0.6/√m` roughly flat in J, so its
*relative* error is ≈ `0.6/(J√m)` and it becomes uninformative below
J ≈ a few/√m (p=14 → J ≲ 0.02; p=20 → J ≲ 0.002).

**The honest summary:** `jaccard_similarity` is a low-variance, monotone, but
biased transform of set Jaccard; the reconstruction is near-unbiased but
noisier and censored at 0. Use the reconstruction when you need *magnitude*
(comparison against bedtools or another tool's absolute values); use
`jaccard_similarity` when you need *discrimination* (ranking, clustering,
nearest-neighbour). Report both; neither is strictly better.

Reproduction — generates the data, collects bedtools ground truth, sweeps
hammock, and prints every table above:

```bash
ml bedtools
python experiments/bedtools_benchmark/estimator_compare.py \
    --tmp-dir /tmp/ec --out results/estimator_compare_full.csv \
    --reps 3 --intervals 1000 --precisions 12,16,20,24
```

Per-pair rows land in the `--out` CSV (`experiments/bedtools_benchmark/results/`
is a symlink to scratch, so that file is not in git — rerun to regenerate).

Verified while building this: `bedtools jaccard` **merges each file
internally** before counting (three identical 100 bp intervals vs a 50 bp B
gives 0.5, not 0.33), and Mode B's `add_points_impl` (`cpp/src/stride.cpp:133`)
walks `pos = start; pos < end` half-open and hashes per position, so
within-file duplicates dedup in the HLL. **The two tools sketch the same
universe** — the comparison is not confounded by merge semantics.

**Known limits of this experiment, to fix before publishing any of it:**

- Synthetic data; each B file reuses a controlled fraction of its paired A
  file's intervals *verbatim*. (This turns out not to bias the estimator
  comparison — both estimators are functions of the bp *sets* only, so
  spatial arrangement is irrelevant — but it is not how real BED files share
  signal.) Repeat on the Maurano DHS corpus.
- The A×B cross product means most pairs sit at J ≈ 0.004; 66/90 rows fall in
  J < 0.05. The summary MAE is therefore close to just reporting `c`.
- No uncertainty quantification. Each file is sketched once at the default
  seed, so all pairs sharing a file share its sketch noise — the effective n
  is the number of *files*, not the number of pairs. Needs a seed sweep with
  CIs.
- Size ratio is pinned at ≈1 and λ is not recorded per row, though both
  govern `c`. Results should be indexed by λ, not by precision.

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
