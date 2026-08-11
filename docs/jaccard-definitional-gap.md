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
(`cpp/src/hll_sketch.cpp:63-75`):

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
floor. Its exact value at high load and equal cardinalities is
`c = Σ_k p_k² / (1 − p_0²) = 0.1699` (measured 0.16992 at p=18, n=2×10⁷).
The `Σ_k (2^−k)² = ⅓` figure sometimes quoted for this is a *different*
conditional — the probability two registers tie **given both hold exactly one
element**. That is not the λ→0 limit of `c`: as λ→0 the floor goes to *zero*
(0.0162 at λ=0.1, 0.0017 at λ=0.01), because almost every active register is
one-sided and cannot tie at all. ⅓ is not approached at any λ; it overstates
the floor ~2× at high load and by orders of magnitude at low load.

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

Analytically `c = Σ_{k≥1} p_k² / (1 − p_0²)` where
`p_k = e^(−λ2^(−k)) − e^(−λ2^(−(k−1)))` is the probability a register lands on
value `k`. Both restrictions matter: `k = 0` is excluded because
`matching/active` counts only registers where at least one side is non-zero,
and the `1 − p_0²` denominator conditions on that same event. (Summing from
`k = 0` gives 4.53 at λ = 0.1 — not a probability at all.) At λ = 92 the sum is
0.169, matching the measured 0.17. Note this is well below the `Σ_k (2^−k)² = 1/3`
figure quoted as the collision bound elsewhere — 1/3 is a different
conditional (see above), not approached at any realistic λ. `c` is
effectively a **step**: 0.1699 from λ ≈ 5 upward, flat to four decimals (any
log-periodic oscillation has peak-to-peak amplitude ~1.2×10⁻⁵), and collapsing
below that — 0.1681 at λ=3, 0.1588 at λ=2, 0.1178 at λ=1, 0.0456 at λ=0.3,
0.0162 at λ=0.1. These are the *conditional* values `Σp_k²/(1−p_0²)`, matching
what `matching/active` computes; the unconditional `Σp_k²` runs up to 5×
smaller at low λ and is not the right quantity. The practical consequence is
that heterogeneity only matters when an analysis *spans the λ ≈ 1–5 knee*; a
corpus entirely on the flat is internally consistent.

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
| A ∩ B = ∅, \|A∪B\| ≪ m | 0.0 | 0.0 (no active register pairs to tie) |
| A ∩ B = ∅, \|A∪B\| ≫ m | 0.0 | > 0 (chance ties dominate) |
| 0 < J < 1, large union | J | f(J) > J |

The mapping `f(J)` is monotonic **at fixed `(|A|, |B|)`** (more shared
elements → more shared argmaxes → more matching registers), but it is not the
identity — and it is a *different* mapping for each cardinality ratio. So
ordering is preserved when comparing pairs of the same geometry and **not**
in general: at p=16 a symmetric pair at J=0.100 and an asymmetric pair
(|A|=2×10⁵ ⊂ |B|=10⁶) at J=0.200 score `J_re` = 0.2645 and 0.2630 — a 0.0015
gap against a 2× difference in true similarity, and in the *wrong direction*.
(The symmetric value depends weakly on n: 0.2620 at |A|=|B|=2×10⁵, 0.2645 from
5×10⁵ up. Measured per-pair sd at p=16 is ≈0.0015, so the separation is about
1σ — recoverable in principle from many replicates, not from one comparison.) See "Rank fidelity" below.

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
through `pairwise_metrics_hll` (`bindings/_core.cpp`), which
uses **inclusion–exclusion** — `|A| + |B| − |A ∪ B|`, Ertl estimator on each,
union by register-wise max, clamped to `>= 0`. That path has no
chance-agreement term, so those columns estimate the true set quantities.
On disjoint inputs at p=16, n=2×10⁵: `jaccard_similarity` = 0.168 while
`containment_AB` = 0.000.

Consequently a set-Jaccard estimate is recoverable, with no rerun and no code
change, from any output that carries the containment block — but **not** from
output written before 2026-05-14, which has the orig `containment` placeholder
instead (see "Old CSVs cannot be back-filled" below). From `C_AB = I/|A|` and
`C_BA = I/|B|`:

```
J_set = I/(|A| + |B| − I) = 1 / (1/C_AB + 1/C_BA − 1)
```

Measured against `bedtools jaccard` on 90 synthetic pairs spanning
J = 0.003 … 0.99 (matched row sets; I-E clamped-to-zero rows scored as 0,
which is exactly what the shipped `jaccard_similarity_ie` column reports):

| p | MAE `jaccard_similarity` | MAE reconstructed `J_set` | ratio | I-E clamped |
| --- | --- | --- | --- | --- |
| 12 | 0.15173 | 0.00822 | 18× | 25/90 |
| 16 | 0.15173 | 0.00207 | 73× | 0 |
| 20 | 0.15164 | 0.00054 | 283× | 0 |
| 24 | 0.03904 | 0.00011 | 361× | 0 |

In the low-J regime real BED comparisons occupy (J < 0.05, mean 0.0114),
register-equality reports **0.181** and the reconstruction reports
**0.0115**.

### Where the two estimators change places

The table above scores **calibration** (is the number right?). It does not
score **resolution** (can the number tell two similar pairs apart?), and on
resolution the ranking does flip — but only in one corner, and precision
closes it.

Resolution is measured here by **Kendall τ against exact `bedtools jaccard`**,
not by a rescaled error sd. The reason is specific: an sd-based comparison has
to divide `sd(err_re)` by something to put the two columns on one axis, and any
binned form of it mixes within-bin true-J variance into the register-equality
column alone (see "the superseded table" below). Kendall τ is invariant under
*any* strictly monotone transform of the estimator, so `f` drops out of it
identically and there is nothing to contaminate.

| J < 0.05 | p=12 | p=16 | p=20 | p=24 |
| --- | --- | --- | --- | --- |
| τ, `jaccard_similarity` | 0.335 | **0.658** | **0.905** | 0.907 |
| τ, `jaccard_similarity_ie` | 0.289 | 0.562 | 0.794 | **0.967** |
| MAE, `jaccard_similarity` | 0.1696 | 0.1696 | 0.1694 | 0.0447 |
| MAE, `jaccard_similarity_ie` | 0.0081 | 0.0019 | 0.0005 | 0.0001 |

Above J = 0.05 both estimators are at or within one discordant comparison of
τ = 1.0000 from p = 16 up, so those strata carry no resolving power and no
winner should be read out of them. (Precisely: in 0.05 ≤ J < 0.2 both reach
1.0000 at p = 16; in J ≥ 0.2 inclusion–exclusion reaches it at p = 16 and
register-equality sits at 0.9804 — 1 discordant comparison out of 102 — until
p = 20. Verified against the generator 2026-08-06.) The whole disagreement lives
below J ≈ 0.05.

Generators: `experiments/bedtools_benchmark/estimator_rank_by_precision.py`
(table, plus a leave-one-replicate-out stability check — the ordering holds in
all three subsamples at p = 16, 20, 24 and flips at p = 12, where both columns
are near-useless at τ ≈ 0.3) and
`paper/estimator_crossover/plot_estimator_crossover.R` (figure).

Inclusion–exclusion is **censored**: it is clamped at 0 (`hll_sketch.cpp:183`),
which fires on 25/90 pairs at p=12, all at low J. Nothing distinguishes
"clamped from a large negative" from "genuinely zero", so conditional on being
non-zero it is biased *upward* near J=0 — a small positive floor of its own
(~0.002–0.014 on truly disjoint pairs, falling with p). Those clamped rows also
tie with one another, which is part of why τ_ie is lowest at p=12. On the
Maurano corpus nothing is clamped at any precision, so the concern belongs to
low-J data specifically.

Rule of thumb: its absolute error is ≈ `0.6/√m` roughly flat in J, so its
*relative* error is ≈ `0.6/(J√m)` and it becomes uninformative below
J ≈ a few/√m (p=14 → J ≲ 0.02; p=20 → J ≲ 0.002).

**The honest summary.** `jaccard_similarity` is a low-variance but biased
transform of set Jaccard, monotone only at fixed cardinality ratio;
`jaccard_similarity_ie` is near-unbiased but noisier and censored at 0. The
earlier framing "neither dominates, report both" overstated the symmetry:
inclusion–exclusion wins on magnitude everywhere by 12–388× (recomputed
directly from `experiments/bedtools_benchmark/results/estimator_compare_full.csv`:
per-(stratum, precision)-cell MAE ratios range from 11.9× at p=12, J≥0.05 to
387.7× at p=24, J<0.05; corrects an earlier "20–1700×" figure that did not
reproduce from this file), and the
register-equality advantage is confined to *ranking*, below J ≈ 0.05, at
p ≤ 20, among pairs of comparable size. Since precision is a knob the user sets
and true J is not knowable in advance, state it as a reading rule:

> **Read `jaccard_similarity_ie`. If your corpus is low-J and you need to rank
> rather than measure, raise `-p` to 24 rather than switching columns.**

Downstream checks on the repo's own low-J corpora are in
`docs/estimator-analysis-findings.md` §9.6–9.7: no published Mode D or
cross-species headline changes under either column, though at k=20 on the
primate H3K4me3 corpus inclusion–exclusion recovers one additional true clade.

<details>
<summary>The superseded resolution table, and why</summary>

The section previously reported `sd(err_re)/(1−c)` against `sd(err_ie)` in
binned true-J units, concluding that register-equality was the more precise
discriminator throughout the saturated regime with inclusion–exclusion taking
over only at p=24. Two defects:

1. **The estimand was contaminated.** Writing `J_re = f(J) + ε_re` and
   `J_ie = J + ε_ie`, the binned statistic is
   `sd_bin(err_re) = sqrt(Var_bin(f(J) − J) + σ_re²)` against
   `sd_bin(err_ie) = σ_ie`. The first term is real *signal* variation inside
   the bin, is seed-invariant, and inflates the register-equality column only.
   Because that floor is λ-independent while `σ_ie ∝ 1/√m`, it grows with
   precision — which is the shape of the crossover the table reported.
   Dividing by `1 − c` also assumes exact affinity, which the concavity
   documented below rules out; the correct divisor is the local slope `f′(J)`,
   running 0.962 → 0.740 across the range.
2. **The size ratio was pinned near 1**, and
   `experiments/bedtools_benchmark/estimator_compare.py` only generates
   equal-size files, so the ratio axis was untested.

Re-deriving that table against an OLS-residualized statistic flipped none of
its twelve winner cells, so the defects were in the method and its
interpretation rather than necessarily in the ordering. The τ table above
reaches a compatible conclusion by a route that has neither defect, which is
why the section is rewritten rather than merely retracted.

</details>

## Rank fidelity is not free — measure it per corpus

The "use `jaccard_similarity` for ranking" advice above holds only within a
fixed pair geometry. Because `c` varies with `|A|/|B|`, the estimator applies a
*different* transform to each pair, and ordering can invert. The resolution
limit for a corpus whose worst cardinality ratio is `r` is the true Jaccard at
which a ratio-`r` pair first outscores a *disjoint symmetric* pair — i.e. the
root of `J_re(r, J) = c₁`, `c₁ = 0.1699`:

| ratio `r` | 1.0 | 1.5 | 2.0 | 2.2 | 3.0 | 5.0 | 10 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `c` | 0.1699 | 0.1635 | 0.1520 | 0.1472 | 0.1293 | 0.0969 | 0.0584 |
| ΔJ below which rank is unreliable | 0 | 0.0067 | 0.0192 | 0.0245 | 0.0452 | 0.0861 | n/a |
| affine approx `(c₁−c_r)/(1−c₁)` | 0 | 0.0078 | 0.0216 | 0.0274 | 0.0490 | 0.0880 | 0.1343 |

Both rows are stable in λ (identical to 4 dp from λ=10 to λ=300). The affine
approximation — which earlier revisions of this file printed as *the* answer —
runs 12–15% high, because `J_re` is concave in `J` rather than a straight line
from `c_r` to 1 — it lies *above* the chord everywhere, peaking at +0.0297 at
J ≈ 0.5, with the local slope falling monotonically 0.962 → 0.740. At `r = 10` there is no root: a ratio-10 pair cannot exceed
`J = 1/r = 0.1`, and even at full containment it scores `J_re = 0.1375 < c₁`,
so **every** 10:1 pair ranks below **every** disjoint symmetric pair regardless
of true similarity. Quoting a ΔJ there would name a gap wider than the entire
feasible range.

Measured on the 20-sample Maurano fetal DHS corpus (bp coverage 63.4–139.4
Mbp, so `r` up to 2.20), 190 off-diagonal pairs at p=21:

| statistic | value |
| --- | --- |
| Pearson r | 0.9972 |
| Spearman ρ | 0.9939 |
| Kendall τ | 0.9511 |
| discordant comparisons | 439 / 17,955 (2.45%) |
| largest inverting gap | ΔJ_bedtools = 0.0250 |

The discordance is **systematic, not sketch noise**: the residual from the
affine model correlates at **−0.885** with `\|log(\|A\|/\|B\|)\|` — the further
apart two files are in size, the further the pair sits below the single-`c`
line. (Against the *signed* log ratio it is only −0.31; the effect is symmetric
in which file is larger, so the absolute ratio is the right regressor.) The
magnitude is fit-dependent and the fit must be named: −0.885 is against
`c = 0.180`, the value `estimator_compare.py` fits on its *synthetic* corpus.
Refitting on Maurano itself gives c = 0.1783 (OLS intercept), 0.1968
(constrained LS) or 0.2004 (mean pointwise), and correlations of −0.90, −0.83,
−0.82. The sign, the ordering and the significance are robust across all four;
only the second digit moves. The same statistics recur at p=18 and p=23 (τ = 0.9505, 0.9497)
on independent sketches, which HLL noise would not do. Note also that a
Pearson `r` over a full square matrix includes the self-pairs at (1,1), which
are pure leverage; quote the off-diagonal value.

## The `jaccard_similarity_ie` column (shipped in v0.4.0)

As of v0.4.0 the inclusion-exclusion estimate is emitted directly as
`jaccard_similarity_ie`, so the reconstruction below is no longer needed for
new runs. It is computed from the same `C_AB`/`C_BA` shown here and agrees with
the direct form to within ~1e-12. (Mode D briefly emitted a
`jaccard_similarity_ie_with_ends` counterpart; that whole `_with_ends` family
was removed in v0.6.0 — see CLAUDE.md divergence #8.)

Note the agreement is *not* a cross-implementation check: both sides derive
from the same two containment arrays in the same process, so it tests the CSV
float round-trip. `runner` also clamps containments to 1.0 before deriving
`jac_ie` but emits them unclamped, so compare with a tolerance rather than
expecting bit-equality.

Two things to know:

- **An exact `0.0` means "clamped or empty", never "measured zero."** The
  clamp is in `pairwise_metrics_hll` (`bindings/_core.cpp`),
  one level below the column.
- **Old CSVs cannot be back-filled.** Output written before 2026-05-14 carries
  the orig `containment` placeholder (constant 1.0) rather than the
  `containment_AB`/`containment_BA` block, so neither the column nor the
  reconstruction is recoverable from it. That includes the archived
  interval-mode A/B/C results.

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

- Implementation: `cpp/src/hll_sketch.cpp:63-75` (Jaccard formula)
- Algorithm parity note: `CLAUDE.md` ("register-equality Jaccard")
- Original Python source: `hammock` (orig) `lib/hyperloglog.py:640`,
  `estimate_jaccard_registers`
- Ertl, O. (2017). *New cardinality estimation algorithms for
  HyperLogLog sketches.* arXiv:1702.01284 — used for the cardinality
  estimator (`ertl_improved_estimate`), but the Jaccard formula
  hammock uses is the simpler register-equality variant inherited
  from the original Python tool, not Ertl's joint MLE.
