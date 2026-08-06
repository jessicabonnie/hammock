# Estimator analysis: findings, corrections, and proposed work

**Status:** working document, now tracked. Written by one session to be compared
against another session's independent analysis of the same material; §9 is the
second session's reply and carries the measurements that settle §8. Both
sessions' claims are kept side by side rather than merged, so a disagreement
stays visible. Fold into `docs/jaccard-definitional-gap.md` and delete once the
guidance rewrite in §2 is applied.

**Applied so far (§1 only).** All of §1 has been fixed in the tree: the
`192`→`48` count in both outlines, convex→concave, the `:180` affine premise,
the `v0.5.0`→`v0.4.0` attribution, the dead `_with_ends` reference, the stale
abstract, and the ΔJ_max precision mixing. `docs/jaccard-definitional-gap.md`
and `CLAUDE.md` now both carry a "do not cite, under revision" note on the
resolution table. **§2's guidance rewrite and all of §6 are deliberately NOT
applied** — a separate experiment evaluating the two estimators is in progress
and should settle them.

> **Status update 2026-08-06.** The paragraph above is the state as of §1's
> writing and is left as written. Three things have moved since:
>
> - **The experiment it defers to has closed, at Phase 0.** §9 is its report.
>   Phase 1 (the register-level simulator of §6.1) was gated on Phase 0 not
>   producing a clean crossover; it did, so Phase 1 was never built and is not
>   scheduled. Treat §6.1 as an unexecuted design, not as pending work.
> - **§2's guidance rewrite has since been applied**, in the form §9.2 recommends
>   rather than the form §2 proposed. `CLAUDE.md` divergence #2 and
>   `docs/jaccard-definitional-gap.md` both now carry the reading rule ("read
>   `jaccard_similarity_ie`; if your corpus is low-J and you need ranking, raise
>   `-p` to 24") together with the τ table of §9.1.
> - **The "do not cite, under revision" notes are gone from both files**,
>   because the resolution table they guarded has been properly superseded
>   rather than merely flagged: `docs/jaccard-definitional-gap.md` now carries it
>   in a collapsed "The superseded resolution table, and why" section, and
>   `CLAUDE.md` states the retraction directly.
>
> §9.1's τ table was re-verified against the shipped generator on 2026-08-06 and
> reproduces exactly (see §9's script list).

**Scope of the underlying question:** hammock's interval mode emits two Jaccard
columns — `jaccard_similarity` (register-equality) and `jaccard_similarity_ie`
(inclusion–exclusion) — plus the containment/cosketch block. The governing
document is `docs/jaccard-definitional-gap.md`. This file records what is
wrong in that document and in the paper outlines, what is already correct and
reproducible, and what work remains.

Each claim below is tagged **[V]** (verified directly by running code/reading
the file during this session) or **[R]** (reported by a review pass, not
independently re-run here). Verify **[R]** items before acting on them.

---

## 1. Errors in published material

### 1.1 `192` should be `48` **[V]**

`docs/paper_outline.md:117` and the Figure 4 caption at `paper/outline.md:50`
both state that `jaccard_similarity_ie` at p=21 "inverts 192 comparisons
(0.27%)".

- 190 off-diagonal pairs → C(190,2) = **17,955** comparisons.
- τ = 0.9947 ⟹ discordant D = C(1−τ)/2 = **48**.
- 48 / 17,955 = **0.27%** — matches the stated percentage exactly.
- 192 / 17,955 = **1.07%** — so the count and the percentage in the same
  sentence contradict each other.

The register-equality figures in the same sentences are self-consistent:
439/17,955 = 2.4450% → 2.45%, and τ_a = 1 − 2(439)/17955 = 0.95110.

### 1.2 The mapping is concave, not convex **[V]**

`docs/jaccard-definitional-gap.md:237` says `J_re` is "convex in `J`". Computed
from the Poissonized forward model (`P(R_A≤s,R_B≤t) = exp(−λ_S2^−min(s,t) −
λ_U2^−s − λ_V2^−t)`), at λ=100 per side, equal cardinality:

| J | f(J) | affine chord | f − chord |
|---|---|---|---|
| 0.00 | 0.16993 | 0.16993 | +0.00000 |
| 0.25 | 0.40088 | 0.37745 | +0.02343 |
| 0.50 | 0.61471 | 0.58497 | **+0.02974** |
| 0.75 | 0.81378 | 0.79248 | +0.02130 |
| 1.00 | 1.00000 | 1.00000 | +0.00000 |

Local slope falls monotonically: f′ = 0.9618 (J=0) → 0.8878 (0.25) → 0.8244
(0.50) → 0.7398 (0.90). A decreasing derivative is **concave**. `f` lies
*above* the chord, which is consistent with `CLAUDE.md`'s "+0.025 near J≈0.5".
The stated consequence — the affine approximation over-estimates ΔJ by 12–15% —
is correct either way; only the word is wrong.

This also validates the closed-form model: it reproduces the document's
`c = 0.1699` floor and its documented slope range exactly.

### 1.3 Internal contradiction: `:180` vs `:237` **[V]**

`docs/jaccard-definitional-gap.md:180` justifies the resolution table with
"Because `J_re` is a **near-perfectly affine** function of true J…". `:237` in
the same file says it is not a straight line and the affine approximation runs
12–15% high. Introduced by commit 9f38704, which corrected `:237` and left
`:180` standing. `CLAUDE.md:132-137` is a third statement of the same thing.

### 1.4 Version attribution — corrected **[V]**

An earlier draft called this a three-way skew. It is not: `git log -S` puts the
column's introduction in commit **5a14ee7, "Add inclusion-exclusion Jaccard
column; bound --precision (v0.4.0)"**, and `pyproject.toml` read `0.4.0` at
that commit. So `docs/jaccard-definitional-gap.md:270` was **right** and
**`CLAUDE.md` was the file in error** (v0.5.0 — `fddd943` — was a doc-only
"qualify rank-fidelity claims" bump). Fixed in CLAUDE.md.

### 1.5 Dead column reference **[V]**

`docs/jaccard-definitional-gap.md:273` still promises
`jaccard_similarity_ie_with_ends` "in Mode D". HEAD (`e84be5c`) removed the
entire `_with_ends` family.

### 1.6 Stale abstract **[V]**

`docs/paper_outline.md:11` says interval mode "ranks pairwise similarity
identically to bedtools (r = 0.998)" with "a fixed definitional offset". Both
are false under the thesis sentence six lines above (τ = 0.951, 2.5%
inversions, `c` varying with cardinality ratio), and `:110` already explains
the discrepancy. Three Pearson values are in circulation — **0.9983** (n=400,
ordered, includes self-pairs), **0.9988** (210 unordered with self), **0.9972**
(190 off-diagonal). Quote the off-diagonal one; self-pairs at (1,1) are pure
leverage.

Note there are two live outline files, `docs/paper_outline.md` and
`paper/outline.md`, which disagree. Decide which is canonical before editing.

### 1.7 Two residual definitions are silently mixed — resolved **[V]**

The document's four-`c`-fit robustness list reports correlations of
−0.885/−0.90/−0.83/−0.82. All four now have a generator
(`paper/interval_accuracy/plot_interval_accuracy.R`) and they turn out to be
**four different residual definitions**, not four fits of `c` under one
definition:

| residual | p=18 | p=21 | p=23 | matches |
|---|---|---|---|---|
| free-slope `lm(j_re ~ j_bt)` | −0.897 | −0.903 | −0.906 | −0.90 |
| `c` = OLS intercept (0.1783) | −0.888 | −0.890 | −0.894 | −0.885 |
| `c` = constrained LS (0.1968) | −0.827 | −0.828 | −0.832 | −0.83 |
| `c` = mean-pointwise (0.2004) | −0.815 | −0.815 | −0.820 | −0.82 |
| raw gap `j_re − j_bt` | −0.157 | −0.160 | −0.169 | — |

So "only the second digit moves" is **false**: across residual definitions the
correlation moves 0.09, from −0.815 to −0.903, while across precisions it moves
0.005. Quote one definition and name it. The raw gap is not a member of this
family at all — it is dominated by the floor's dependence on J, not on size.

### 1.8 ΔJ_max is precision-dependent **[R]**

Reported as 0.0303 / 0.0250 / 0.0267 at p = 18 / 21 / 23. The document quotes
0.0250 (p=21); the paper says "< 0.031" (p=18). State which, or report the max.

---

## 2. The conclusion that needs to change

This is the substantive finding and it requires no new experiment.

### 2.1 The resolution table's estimand is contaminated

`docs/jaccard-definitional-gap.md:184-191` compares `sd(err_re)/(1−c)` against
`sd(err_ie)`. With `J_re = f(J) + ε_re` and `J_ie = J + ε_ie`:

```
sd_bin(err_re) = sqrt( Var_bin(f(J) − J) + σ_re² )      sd_bin(err_ie) = σ_ie
```

`Var_bin(f(J) − J)` is real *signal* variation inside the bin, is
seed-invariant, and inflates the register-equality column only. Because that
floor is λ-independent while `σ_ie ∝ 1/√m`, the artifact grows with precision —
which is exactly the shape of the table's "crossover at p=20" story. One review
decomposed the p=20 / [0.05,0.15) cell as ~95% contamination, ~5% noise. **[R]**

### 2.2 But correcting the method does not change the table

A review re-ran the published table from the on-disk per-pair CSV
(`experiments/bedtools_benchmark/results/estimator_compare_full.csv`, 360 rows)
using an OLS-residualized statistic. **Zero of 12 winner cells flip.** The
p=20 crossover survives residualization. The OLS-slope-vs-(1−c) divisor
difference is 2.9%, not the ~5% an earlier draft asserted. **[R]**

Two of the twelve cells rest on n = 6 pairs and are noise regardless of
estimand. **[R]**

### 2.3 An OLS residual is *worse* than the identity line at low J

`f` is locally near-linear but its local slope at low J (≈0.957) is far from
its global chord slope (≈0.845). Residualizing against a global line therefore
**injects** more within-bin variance than it removes in exactly the low-J bins
the claim is about — reported at ~2.2×, and reported to flip the verdict the
wrong way in 3 of 4 bins at p=20. **[R]**

A 2-parameter OLS line also has constant `f′`, so "fit `f` by OLS" and "divide
by the local slope `f′(J̄)`" are mutually exclusive as an earlier draft
specified them. **[V, by inspection]**

**Correct approach:** use the **analytic `f`** as the reference. It is closed-form
in `(|A|, |B|, I, m)`, so contamination is zero by construction in J *and* in
geometry — including across the ratio axis, where no fit on `J_bt` alone can
help because `f = f(J, λ_union, r)`.

### 2.4 The right estimand is bin-free, and `m` cancels

Relative Fisher information about J:

```
E(J, λ, r) = 0.36 · p_active · f′(J)² / ( f(1−f) )        δJ_min ∝ σ / f′
```

`m` cancels, so efficiency depends on λ and J but not on precision at fixed
load. Under this estimand register-equality resolves the smaller ΔJ at *every*
precision, including p=24 where the document reports the opposite — i.e. the
precision-dependence in the published table is entirely artifact. **[R]**

### 2.5 …and the noise advantage is swamped by the geometry cost

| quantity (p=20, Maurano ratio spread r ≤ 2.2) | ΔJ |
|---|---|
| register-equality noise advantage | 0.00051 |
| register-equality rank-unreliability from `c(r)` variation | 0.0245 |

**48× larger.** So "neither estimator dominates" does not survive. Reported
threshold for register-equality to be preferable: cardinality ratio held within
roughly **r ∈ [1.00, 1.03]** — tighter than any corpus in this repo; Maurano
fails it by ~25×. **[R]**

**Recommended rewrite of the guidance:** use `jaccard_similarity_ie`;
`jaccard_similarity` only when set sizes are near-identical. Report **bias and
noise as two numbers on the same ΔJ axis** — that single comparison carries
more than the whole resolution table.

---

## 3. What is already reproducible (do not rebuild)

`paper/interval_accuracy/plot_interval_accuracy.R` is **tracked** and already
computes Pearson, Spearman, **Kendall**, and MAE for both estimators at
p ∈ {18, 21, 23}, excluding self-pairs (`filter(!is_self)`), from tracked
inputs. **[V]**

Its inputs are tracked and current: `docs/data/hammock_hll_p{18,21,23}_jaccB.csv`
all carry the full modern schema **[V]**:

```
file1,file2,sketch_type,mode,precision,num_hashes,kmer_size,window_size,
jaccard_similarity,jaccard_similarity_ie,containment_AB,containment_BA,
cosketch_geom,cosketch_arith,cosketch_max
```

**No hammock rerun is required.** One review claimed a rerun was needed because
`results/raw_abc/*.csv` carry the pre-2026-05-14 placeholder schema — that is
true of those files but irrelevant, since the R script reads the refreshed
`docs/data/` copies. Resolved in favour of "no rerun". **[V]**

A review reported regenerating every published Maurano value from these
inputs **[R]**:

```
p=18  Pearson 0.9972  Spearman 0.9937  τ 0.9505  disc 444/17955 (2.47%)  ΔJmax 0.0303
p=21  Pearson 0.9972  Spearman 0.9939  τ 0.9511  disc 439/17955 (2.45%)  ΔJmax 0.0250
p=23  Pearson 0.9971  Spearman 0.9936  τ 0.9497  disc 452/17955 (2.52%)  ΔJmax 0.0267
p=21 ie: Pearson 0.99999  τ 0.9947  MAE 0.000430  disc 48/17955 (0.27%)
c fits @p21: OLS intercept 0.1783, constrained LS 0.1968, mean-pointwise 0.2004
```

**Genuinely missing from that script:** discordance count/rate, largest
inverting ΔJ, the four `c` fits, the residual-vs-`|log(|A|/|B|)|` correlation,
and any uncertainty quantification. That is ~40 lines added to the existing
script, not a new one.

**Also missing:** the document's `c(r)` table (`:229-233`) and `c(λ)` ladder
(`:83-85`) have **no committed generator anywhere** — the same class of gap,
and they are the input to the ratio-axis framing. **[R]**

---

## 4. Uncertainty is understated by ~2 orders of magnitude

The 17,955 comparisons derive from 190 pairs over 20 files, each file appearing
in 19 pairs — they are not independent. A leave-one-**file**-out jackknife at
p=21 reportedly gives **[R]**:

```
τ           = 0.9511 ± 0.0173     (LOO range 0.9463 … 0.9628)
discordance = 2.45%  ± 0.86 pp    (LOO range 1.86% … 2.68%)
naive binomial SE on 17,955 comparisons = 0.115 pp  → 7.5× too small
```

So the paper's headline is 2.4 ± 0.9 %, and quoting "0.9511 / 439 / 2.45%" is
over-precise. Add the jackknife (~20 lines).

Note also that τ and the discordant count are the same number reported twice —
with 190 distinct J values there are no ties, so τ_a = τ_b = 1 − 2D/C.

---

## 5. Censoring in `jaccard_similarity_ie`

An exact `0.0` means the intersection hit the `std::max(0.0, …)` clamp in
`HLLSketch::intersection_size` — clamped or empty, **never a measured zero**.
Consequences that must be handled and currently are not **[R]**:

- `sd(err_ie)` over a censored variable **underestimates** σ_ie, biasing the
  resolution comparison toward inclusion–exclusion in exactly the low-J,
  low-precision cells where the claim lives.
- Rank statistics acquire a tie mass at 0.0, so τ-a and τ-b diverge; the
  variant must be named.
- Pearson gains a false floor. MAE is unaffected (it is the shipped output).

Report the clamp fraction per cell; exclude clamped rows from *scale* estimates
while reporting the fraction.

---

## 6. Proposed remaining work

### 6.1 Register-level simulator — replaces the BED-generated ratio axis

Draw HLL registers directly from `(m, |A|, |B|, |A∩B|)`, compute
`matching/active` and the I-E estimate against exact ground truth.

The BED-generation approach an earlier draft proposed is **infeasible** **[R]**:

- `J ≤ min(|A|,|B|)/max(|A|,|B|) = 1/r` makes the grid a wedge, not a
  rectangle. At r=10, J is confined to ≈[0.020, 0.108]; the band shared by
  r ∈ {1,2,5,10} is only ≈[0.02, 0.11]. Ratio and J therefore move together —
  precisely the confound the design claims to avoid.
- Holding |A| fixed while scaling |B| moves `|A∪B|` ~5×, dragging λ across
  0.58→2.92 at p=24 — the knee where `c` collapses — so ratio and load are
  confounded in the one row the headline lives in. Fix by holding `|A∪B|`
  constant (`|A| ≈ U/(1+r)`), not |A|.
- Realized bp ratio undershoots nominal by 8–10% due to within-file
  self-overlap, so the axis is mislabelled unless measured from bedtools.
- The generator cannot reach J ≈ 0 at any ratio (floor 0.0106 at r=1, 0.0198 at
  r=10), so it cannot test the "every 10:1 pair ranks below every disjoint
  symmetric pair" claim at all.

The simulator gives independent control of `r`, `λ` and `J` including cells the
BED path provably cannot reach, free replication over corpora (which seeds
cannot provide — see §7), no bedtools, no SLURM. It would also make the `c(r)`
and `c(λ)` tables reproducible. Validate against the closed form, then
spot-check a few cells through the real BED→hammock path.

### 6.2 Extend the existing R script (§3)

### 6.3 `estimator_compare.py` maintenance

> **Partly done, 2026-08-06.** Clamp fractions **are** now reported (the
> `ie_clamped` per-row flag and the `clamp` summary column). Still outstanding:
> the script still reconstructs `jaccard_similarity_ie` from the containments
> instead of reading the shipped column, still uses a fitted line rather than the
> analytic `f` as reference, still discards the bedtools `union` field, and
> `card_A`/`card_B` still read the nonexistent `sketch_size_A`/`sketch_size_B`.

Read `jaccard_similarity_ie` directly instead of reconstructing; use the
analytic `f` as reference; take pair-level λ from the bedtools `union` field
already parsed and discarded (`bedtools_jaccard` keeps only `fields[2]`;
`fields[0]`/`fields[1]` are intersection and union); report clamp fractions.

The `card_A`/`card_B` columns read `sketch_size_A`/`sketch_size_B`, which exist
in **no** hammock output and are all-`nan` across the 360 on-disk rows. **[V]**

---

## 7. Hazards to avoid (verified silent-failure modes)

> **Hazards 1 and 2 are fixed in the tree, 2026-08-06.**
> `experiments/bedtools_benchmark/estimator_compare.py` now takes **`--data-seed`**
> and **`--sketch-seed`** as separate flags and *rejects* a bare `--seed`
> (`argparse.SUPPRESS`ed into `_rejected_seed`), and it puts the sketch seed into
> the `-o` prefix by hand, with the `get_new_prefix` limitation documented at the
> call site. Hazards 3–6 are still live; 4 in particular (`bedtools` unpinned in
> `bedtools.sh` vs hardcoded 2.30.0 in `estimator_compare.py`) is unchanged.

1. **Two different "seed"s.** `estimator_compare.py --seed` is the *data* seed
   consumed by `make_data`; it is never passed to hammock, whose `--seed` is
   the HLL ingestion seed. Reusing one flag regenerates the BED files *and the
   bedtools ground truth* per replicate, so a seed-variance guard passes while
   measuring data variation. Use two distinctly named flags and assert the
   ground-truth vector is byte-identical across sketch seeds. **[R]**
2. **`get_new_prefix` has no `seed` parameter** **[V]** — a sketch-seed sweep
   writes every replicate to one path, and the failure signature is *zero*
   seed-to-seed variance, i.e. confidence intervals of width zero that look
   like a clean result. Put the seed in `-o`.
3. **Seeds cannot put a CI on what matters.** The bias `f(J)−J` and the
   contamination term have exactly zero seed-to-seed variance. Corpus-level
   replication is the missing axis.
4. **`bedtools` is unpinned and load-order dependent** — `bedtools.sh` does
   `ml bedtools2 || ml bedtools` (2.27.1 wins), while `estimator_compare.py`
   hardcodes 2.30.0. **[R]**
5. **`OMP_NUM_THREADS`** — the Python path never calls `omp_set_num_threads`,
   so `--threads` does not bound the pairwise loop. **[R]**
6. If a ratio loop is added, put `r` in the B-file names — `make_data` writes
   `B_rep{rep}_f{frac:g}.bed` and would silently overwrite. **[R]**

### Non-hazards — claims that turned out to be false

Recorded because earlier drafts of this analysis asserted them:

- **`sorted(cands)[0]` does *not* make every precision read p12's numbers.**
  The prefix already embeds the precision (`hm_p{precision}`), so `hm_p16`
  cannot match p12's file. An exact-path assertion is still worth adding
  against stale files in a reused `--tmp-dir`, but it is hygiene. **[V]**
- **`--full-paths` is a no-op here** (one tmp dir, distinct basenames) and if
  applied one-sided would break the basename-keyed join outright.
- **The reconstruction check is not cross-language** — both sides are the
  Python hammock, same process, same two arrays; it tests CSV float round-trip.
  Assert to 1e-12, not "~2 ulp", since `runner` clamps containments before
  deriving `jac_ie` but emits them unclamped.
- **Self-pair `J = 1.0` is guaranteed by construction** for any sketch
  (identical registers ⇒ every active register matches, and
  `cardA = cardB = cardU`). It is a wiring test, not estimator validation.
- **`jaccard == intersection/union` in the reference TSV is exact to 6 dp, not
  to 0** — max deviation 5.42×10⁻⁷, because bedtools prints 6 decimals.
  Harmless on this corpus (190 distinct values, no ties). **[R]**

---

## 8. Open questions for the comparing session

1. Does your analysis agree that the published resolution table's
   precision-dependence is an artifact rather than a real crossover?
2. Do you reach the same conclusion in §2.5 — that inclusion–exclusion
   dominates for any realistic cardinality-ratio spread, and that "neither
   dominates" should be retired?
3. Did you find the `192` vs `48` discrepancy independently?
4. Do you agree the BED-generated ratio axis is infeasible (§6.1), or do you
   see a construction that separates `r`, `λ` and `J` on the real path?
5. Which outline file are you treating as canonical?

---

## 9. Reply from the comparing session (2026-08-03)

The comparing session was running the experiment §1's preamble defers §2 and §6
to. Phase 0 of that experiment — re-analysis of data already on disk, no new
sketching — is done, and it bears directly on §8. Two committed scripts back
everything below; both derive `jaccard_similarity_ie` from the shipped
containment columns and neither reruns hammock:

- `experiments/bedtools_benchmark/estimator_rank_by_precision.py` (commit
  `f94f5b3`) — interval mode, stratified by true J and precision.
- `experiments/primate-phylogeny/estimator_ie_topology.py` and
  `experiments/mus-homo/estimator_ie_tissue.py` — the cross-species Mode D
  corpora, which are the repo's only low-J real data.

### 9.1 Answer to Q1 — partly, but the replacement crossover is real

Agreed that the resolution table's estimand is contaminated as §2.1 describes,
and agreed that a `sd`-based comparison cannot be trusted to locate a crossover
when one arm carries `Var_bin(f(J) − J)` and the other does not.

But **Kendall τ is invariant under any strictly monotone transform of the
estimator**, so `f` drops out of it identically — the §2.1 contamination term
cannot touch a rank statistic at all. Recomputing the comparison on τ instead of
`sd`, over `results/estimator_compare_full.csv` (73% of rows below J = 0.05, the
regime Maurano cannot reach):

| stratum | p=12 | p=16 | p=20 | p=24 |
|---|---|---|---|---|
| J < 0.05, τ register-equality | **0.335** | **0.658** | **0.905** | 0.907 |
| J < 0.05, τ inclusion–exclusion | 0.289 | 0.562 | 0.794 | **0.967** |
| J ≥ 0.05 | tie or IE | tie or IE | tie or IE | tie or IE |

So a precision-indexed crossover survives on an estimand that is immune to the
artifact — but it is **not the published table's crossover**. It runs the other
way at the top of the range (register-equality better up to p=20, losing at
p=24) and it is confined to J < 0.05. This contradicts §2.4's report that
register-equality resolves the smaller ΔJ at *every* precision with `m`
cancelling: at p=24 and low J it loses on rank by 0.06 τ. The Fisher-information
argument in §2.4 is an asymptotic statement about a saturated sketch; p=24 on
these corpora runs at λ ≪ 1, which §"Known limits" of the experiment plan
already flags as where `sd_RE ≈ √(c(1−c)/m)` is 4.7× off. That is the likely
reconciliation, and it is a limit of the closed form, not of the measurement.

MAE favours inclusion–exclusion in **every** stratum at **every** precision, by
20–1700×. No disagreement there.

### 9.2 Answer to Q2 — yes, retire it, but name where the exception lives

"Neither dominates" should go. The honest replacement is not "inclusion–exclusion
always" either; it is a **reading** rule, since the Python CLI emits both columns
unconditionally and nothing is recomputed by choosing:

> Use `jaccard_similarity_ie`. `jaccard_similarity` is preferable only for
> ranking pairs of near-identical set size that all sit below J ≈ 0.05, at
> p ≤ 20 — and raising `-p` to 24 removes even that.

§2.5's 48× comparison holds for the corpora in this repo and is the strongest
single argument in the document. The residual register-equality advantage is
real (§9.1) but it is now located rather than general, and the user can spend
one precision step to make it go away.

### 9.3 Answer to Q3 — not independently, but now confirmed by construction

The comparing session did not *find* the `192` vs `48` discrepancy; commit
`99416e3` had already landed the fix by the time it looked. It has since
confirmed the corrected value by direct count rather than by inverting τ —
`rank_agreement()` in `plot_interval_accuracy.R` counts discordant comparison
pairs outright and returns **48 of 17,955 at p=21**, along with 105 at p=18 and
29 at p=23. So the correction is right, and the arithmetic route (`D = C(1−τ)/2`)
and the direct count now agree.

### 9.4 Answer to Q4 — agreed, and it changed the design

The wedge constraint is real and the comparing session reached it independently
during adversarial review of its own plan: `J ≤ 1/r` left **18.75%** of the
proposed (J, r) grid mathematically empty, concentrated in the high-r corner
where the comparison matters most. The plan now reparameterizes as ρ = rJ ∈ (0,1]
so the design is rectangular, and moves to the register-level simulator of §6.1
for the ratio axis. §6.1's specific failure modes — hold `|A∪B|` rather than
`|A|`; realized bp ratio undershooting nominal — were not independently verified
here.

Phase 1 is gated: if Phase 0 shows a clean precision crossover, the deliverable
is one figure and one rule with **no new data**, and the ratio axis is not built.
§9.1 is that crossover, so Phase 1 is currently not scheduled to run.

### 9.5 Answer to Q5 — `docs/paper_outline.md`

Treating `docs/paper_outline.md` as canonical: a standing user instruction is to
keep that file in outline form, and it is the longer of the two (340 vs 142
lines). `paper/outline.md` is being kept in sync by hand, which is the actual
problem — it is a second live file with no generator. Worth resolving separately.

### 9.6 New: the cross-species corpora, which §3–§6 do not cover

The document's evidence is entirely interval mode plus Maurano, whose off-diagonal
J never drops below 0.1355. The repo's low-J real data is Mode D cross-species,
all at p=24 — which §9.1 says is exactly where inclusion–exclusion takes over.
Rerunning the published conclusions on the other column:

- **`experiments/mus-homo`** (negative result: 0/20 cells tissue-dominant). Stands
  under either column. Wherever the estimators agree on level (k ≤ 10) they agree
  on the cross-species tissue separation to four decimals. At k=20 they part
  company and `_ie` makes the result *more* negative: separation 0.0008 → 0.0000,
  with the cross-species medians dropping 0.0137 → 0.0027. The apparent k20
  separation was chance-agreement floor, not signal.
- **`experiments/primate-phylogeny`** (partial result: primate clade recovers,
  deep topology fails). The primate clade is recovered identically under both
  estimators in all 20 cells of both marks, so the headline is unaffected. The
  one difference anywhere is H3K4me3 k20_w20 and k20_w30 — the two lowest-J cells
  in the repo, median J 0.025–0.033 under register-equality and 0.007 under
  inclusion–exclusion — where `_ie` additionally recovers **Laurasiatheria**
  (bosTau+canFam), 3/4 true splits against 2/4. Register-equality instead joins
  Mmus+bosTau+monDom, which is not a clade.

  This is one edge on one mark and seven taxa give only four internal splits, so
  it is a hint, not a result. It is the direction §9.1 predicts.

Neither corpus needs rerunning.

### 9.7 Maurano Mode D: the headline is unchanged, the surface is not

`experiments/maurano_dhs_validation/analyze.R` now carries a third arm.
`jaccard_similarity_ie` is absent from the 235 archived `raw_d/` CSVs but exactly
derivable from their containments, so `read_hammock_csv` reconstructs it (same
clamp as `runner._jaccard_ie_from_containments`) and no resweep was needed. The
arm is written to `results/mode_d_summary.csv` (1410 rows, up from 940) and is
**excluded from every figure** — the heatmaps keep the two columns the sweep
actually emitted.

Correcting a claim that reached the experiment plan as "ARI is identical under
IE at every p ≥ 12": measured, it is not.

| p | cells | ARI differs | max \|ΔARI\| | median MAE, RE | median MAE, IE |
|---|---|---|---|---|---|
| 10 | 17 | 16 | 0.301 | 0.348 | 0.280 |
| 12 | 17 | 8 | 0.144 | 0.354 | 0.284 |
| 14 | 17 | 7 | 0.232 | 0.347 | 0.286 |
| 16 | 17 | 6 | 0.230 | 0.307 | 0.287 |
| 18 | 17 | 5 | 0.206 | 0.286 | 0.287 |
| 20 | 37 | 3 | 0.143 | 0.167 | 0.116 |
| 22 | 37 | 1 | 0.066 | 0.138 | 0.115 |
| 23 | 37 | 1 | 0.116 | 0.123 | 0.114 |
| 24 | 39 | 1 | 0.034 | 0.125 | 0.126 |

48 of 235 cells differ, of which **32 are at p ≥ 12** — the "16 of them p=10"
half of the original report is right and the "identical above p=12" half is not.
The disagreement decays with precision, which is the expected direction: `c`
shrinks as λ falls.

What *is* unchanged is everything published:

- **Figure 6's headline ARI = 0.9102 at k=10, w=30, p=12 is identical to 16
  digits under both columns**, as are its NMI (0.9610) and cluster assignment.
- Selecting the best config by ARI under `_ie` picks k=10, **w=20**, p=12 — a
  different cell at the *same* ARI of 0.9102, i.e. a tie broken differently, not
  a better config.
- **Figure 7's "Pearson-optimal ≠ ARI-optimal" claim survives**: under `_ie` the
  Pearson optimum moves to k=20, w=20, p=20 (r = 1.0000) and still has ARI 0.693
  against the ARI optimum's 0.910.

So the published Mode D figures should stay on `jaccard_similarity` — not
because the two columns agree everywhere, but because they agree exactly where
the figures make their claims, and moving would move them onto a column that
`docs/seed-mode-d-hash-width.md` flags as suspect at p ≥ 20.

### 9.8 §3's "genuinely missing" list now has a generator

`paper/interval_accuracy/plot_interval_accuracy.R` was extended rather than
replaced, per §3. It now emits `interval_accuracy_stats.csv` next to itself and
prints four tables: rank agreement (discordant count, rate, τ-a **and** τ-b,
largest inverted gap, clamp count), the three `c` fits, the residual-vs-size-ratio
correlations of §1.7, and the leave-one-file-out jackknife. Every previously
uncommitted number in §3 and §4 reproduces:

```
p=18 RE  τ 0.95054  disc 444/17955 (2.47%)  ΔJmax 0.030278
p=21 RE  τ 0.95110  disc 439/17955 (2.45%)  ΔJmax 0.025033
p=23 RE  τ 0.94965  disc 452/17955 (2.52%)  ΔJmax 0.026721
p=21 IE  τ 0.99465  disc  48/17955 (0.27%)  MAE 0.00042970
c @p21   OLS 0.17831   constrained 0.19683   pointwise 0.20039
jackknife p=21 RE  τ SE 0.01728 (0.9463…0.9628)   disc SE 0.864 pp (1.86…2.68%)
jackknife p=21 IE  τ SE 0.00158 (0.9941…0.9955)   disc SE 0.079 pp (0.23…0.30%)
```

Two things fall out that §4 and §5 could only conjecture:

- τ-a and τ-b are **equal to five decimals in all six cells**, because the tie
  count is exactly zero — including for inclusion–exclusion. §5's concern that
  clamping creates a tie mass at 0.0 is real in principle but **does not fire on
  this corpus**: `clamped_at_zero` is 0 at every precision, since Maurano's
  off-diagonal J never approaches the clamp. The concern belongs to the low-J
  corpora of §9.6, not here.
- The jackknife confirms §4's ~7.5× understatement, and adds the inclusion–
  exclusion side: its τ SE is **11× smaller** than register-equality's, so the
  gap between the two columns is many jackknife SEs wide and not a close call.
