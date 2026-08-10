> # ⛔ EVERY BEDTOOLS-RELATIVE NUMBER BELOW IS SUPERSEDED (2026-08-09)
>
> **Do not quote any speedup-over-bedtools figure from this file.** Four defects
> were found in the benchmark harness in one day, **all inflating the bedtools
> baseline, i.e. all in hammock's favour**:
>
> 1. `bedtools.sh` ran three processes per pair instead of one (~2.1×).
> 2. bedtools was pinned to 2.27.1, whose jaccard union is order-dependent.
> 3. `ml` module loads ran inside the timed region.
> 4. **The bedtools module exports `LD_LIBRARY_PATH` with 17 directories, of
>    which bedtools uses 4. The dynamic linker searches the other 13 — on GPFS —
>    on every `execve`, and a pairwise workflow is N² execs.** This alone
>    inflated bedtools 2.4–2.8×.
>
> Measured effect at N=512: bedtools 1978 s → 716 s. **The headline speedup falls
> from 27.61× to roughly 10.2×**, and that replacement is itself provisional —
> it rests on a single job and has not been through the checks below.
>
> Also retracted: the "process creation caps near **123 exec/s** and does not
> scale with cores" mechanism. Measured 16-way: **364 exec/s** clean, 196 slow.
> The `md5sum` control that appeared to confirm it ran in the same polluted
> environment, so it confirmed nothing.
>
> **What survives:** anything hammock-internal — mixed-stride vs hash-threshold,
> the precision sweeps, the fused-pass A/B, Mode D threading. Those divide two
> hammock numbers measured in one run and are unaffected.
>
> Status and worklist: `docs/bedtools-baseline-retraction.md`.

# subB / subB-method sweep — results notes

Findings from the 2026-05-12 (3-method × subB × size grid) and 2026-05-13
(bedtools reference) runs. Numbers below come from
`results/summary_synthetic_20260512_144455.csv`,
`results/summary_maurano_20260512_173720.csv`, and
`results/sweep_maurano_bedtools_20260513_132550.csv`. Figures live in
`figures/`.

Design and how-to-rerun live in `README.md`.

> **Which Jaccard column (noted 2026-08-06).** Every accuracy number on this
> page is `jaccard_similarity` — the register-equality column, which is not set
> Jaccard and carries a chance-agreement floor (`CLAUDE.md` divergence #2). The
> sweep predates `jaccard_similarity_ie`, so nothing here is a calibrated
> comparison against bedtools; the bedtools comparison in §4.1 is a *wall-time*
> comparison only. The one place the column choice changes the conclusion — the
> 10k MAE cliff — is called out in its own box below.
>
> **Answered 2026-08-08 — see §4.6.** The IE arm now exists, measured against
> exact bedtools truth on this corpus: `--subB` does **not** meaningfully degrade
> `jaccard_similarity_ie` anywhere in 0.01 ≤ subB ≤ 1, at p=18 or p=23. The
> guidance on this page stands, and its stated reason — that part of the
> register-equality error is the floor moving rather than the sample degrading —
> is now confirmed on real data rather than inferred.
>
> **Binary version.** These runs used a pre-0.7.0 `hammock-cpp`, when 3-column
> output was the default and `--no-metrics` did not exist. `run_sweep.py` has
> since been updated to require ≥ 0.7.0 and to pass `--no-metrics` explicitly,
> so a re-run produces the same 3-column shape and comparable timings — but it
> is a different binary, not the one that produced the tables below.

---

## TL;DR

> **SUPERSEDED 2026-08-09 — re-measured as job 29652709; see the correction
> block below before quoting anything in this file against bedtools.** The
> bedtools baseline every ratio here divides by was taken with a `bedtools.sh`
> that ran three processes per pair instead of one. Ratios *among hammock
> settings* are unaffected; ratios *against bedtools* are all ~1.5× too
> generous.

On real DHS data, `--subB-method mixed-stride --subB 0.1` is **1.64× faster
than `bedtools jaccard`** (was reported as 2.13×) with mean absolute Jaccard
error < 0.001 vs the no-subsample reference. At `subB=1.0` (no subsampling)
hammock-cpp is **1.12× *slower*** than bedtools on this corpus — the earlier
claim of "1.16× faster" is retracted. Twenty files is well below the N ≈ 64
crossover, so the sketch path has not yet amortised its construction cost over
enough pairs; subsampling is what buys the win at this corpus size, and scale
is what buys it without subsampling (27.6× at N=512, `bedtools_benchmark`).

![](figures/headline_maurano_pareto.png)

The three `--subB-method` modes give very similar accuracy, but only
**mixed-stride** actually trades subsampling for speed.

**Reading the zigzag.** The lines aren't sorted by error along the x-axis;
they connect points in subB sweep order (`0.5 → 0.25 → 0.1 → 0.05 → 0.01`,
i.e. each line moves *along the subB axis*, which is not a plot axis). So
whenever MAE isn't monotone in subB, the line doubles back. That's most
visible on mixed-stride: as subB goes from 0.5 → 0.25 the error jumps up,
then drops back at subB=0.1. The cause is that mixed-stride is a
*deterministic* chr-keyed sampler — at each subB it picks a specific stride
length, and the resulting sample sets are not nested across subB values.
Whether a given subB happens to align with the interval/chromosome
structure varies discretely. Per-replicate IQR across 5 replicates is
effectively zero in every cell (algorithm is deterministic given the
seed), so the zigzag is structural, not noise. See
`figures/pareto_variants.pdf` (generated by `pareto_variants.R`) for
path-resorted and no-path versions if you'd prefer to read the same data
without the doubling-back.

![](figures/headline_hammock_vs_bedtools.png)

---

## Speed: mixed-stride is the only method whose wall time tracks subB

![](figures/headline_synthetic_scaling.png)

At `subB=0.1` on synthetic data, mixed-stride gives a 3–4× speedup that
grows with corpus size:

| size | hash-threshold | mixed-stride | single-hash |
|------|----------------|--------------|-------------|
| 10k  | 1.19× | **3.09×** | 1.18× |
| 100k | 1.14× | **3.71×** | 1.12× |
| 1M   | 1.14× | **3.93×** | 1.11× |

The other two methods plateau around 1.1–1.2× regardless of corpus size or
how aggressive subB gets — their per-point gate-hash cost masks the savings
from fewer HLL ingests. Mixed-stride's chromosome-keyed deterministic
stride avoids that per-point hash, so its savings actually compound.

Pushing further down to `subB=0.01` on synthetic 1M, mixed-stride hits
**14×** (vs ~1.4× for the other two). The headline plot (Maurano) and the
diagnostic plots show the same pattern on real data, just at smaller
magnitudes (peak 2.62× vs no-subsample at subB=0.01; that same point is
3.05× vs bedtools — the baseline used for the §4.1 headline figure).

![](figures/synthetic_speedup_vs_nosub.png)

**The subB=0.5 hump:** every method on every corpus is *slower* than no
subsample at `subB=0.5` — speedup 0.64–0.96× across all 12 (method, corpus,
size) cells. Real wins only start at `subB ≤ 0.25` and grow from there.

**Two different mechanisms, not one** (corrected 2026-08-06 — this paragraph
previously blamed both on "the gate's per-point hash cost", which cannot
explain the mixed-stride cells, since mixed-stride pays no per-point hash):

- `hash-threshold`/`single-hash` (0.64–0.68× and 0.70–0.74×) do pay the gate.
  Each position costs an extra `xxh32` plus an unpredictable 50/50 branch,
  while only half the `xxh64`+insert work is skipped — ≈1.5 hash-units of
  work against 1.0. Predicted 1.5×, observed 14.668/9.555 = **1.535×**.
- `mixed-stride` (0.94–0.96× synthetic, 0.82× Maurano) visits only half the
  positions and should therefore be ~2× *faster*. It isn't, because the
  strided loop formats each position with `std::to_chars`
  (`cpp/src/stride.cpp:114`) while the no-subsample path uses the in-place
  ASCII digit-increment fast path, measured at `stride.cpp:127` as ~55% of
  inner-loop instructions. Half the iterations at roughly double the
  per-iteration cost nets ≈1.0×.

(The 0.64–0.96× range is itself a correction: this line previously read
"0.65–0.86×", which excludes the three mixed-stride synthetic cells.)

---

## Accuracy: methods are tied; real data is roughly 2× worse than synthetic

Per-pair MAE vs the subB=1.0 ground truth (averaged across methods):

![](figures/maurano_mae_vs_subB.png)

| corpus    | size  | hash-threshold     | mixed-stride       | single-hash        |
|-----------|-------|--------------------|--------------------|--------------------|
| synthetic | 10k   | 1.1 × 10⁻³ / 2.1 × 10⁻² | 6.0 × 10⁻⁴ / 2.0 × 10⁻² | 9.7 × 10⁻⁴ / 2.0 × 10⁻² |
| synthetic | 100k  | 1.4 × 10⁻³ / 1.2 × 10⁻³ | 6.1 × 10⁻⁴ / 8.4 × 10⁻⁴ | 1.5 × 10⁻³ / 1.8 × 10⁻³ |
| synthetic | 1M    | ~2 × 10⁻⁵ / ~2 × 10⁻⁵ | ~2 × 10⁻⁵ / ~2 × 10⁻⁵ | ~2 × 10⁻⁵ / ~2 × 10⁻⁵ |
| maurano   | 20 DHS | 9.3 × 10⁻⁴ / 1.9 × 10⁻³ | 9.4 × 10⁻⁴ / 1.5 × 10⁻³ | 1.1 × 10⁻³ / 2.3 × 10⁻³ |

(each cell: MAE @ subB=0.1 / MAE @ subB=0.01)

Three observations:

1. **Method spread is small.** Within each (corpus, size) the three methods
   are within ~2× of each other on MAE. Mixed-stride is the recommended
   default for performance reasons, not accuracy.
2. **Real data has 2× worse accuracy** than synthetic at comparable file
   size. Max-error tails are also fatter: max |error| ≈ 10⁻² on Maurano at
   subB=0.01 vs ~3 × 10⁻³ on synthetic 100k. Clumpier spatial structure
   means subsampling has higher per-pair variance.
3. **The 10k cliff doesn't repeat on Maurano.** Synthetic 10k jumps to
   MAE ≈ 0.02 at subB=0.01; the Maurano corpus (per-file size somewhere
   between synthetic 100k and 1M) stays at MAE ≈ 0.002.

> **The 10k cliff is an estimator artifact, not a sampling failure**
> (established 2026-07-31; these MAEs predate `jaccard_similarity_ie`).
>
> **Confirmed on real data 2026-08-08 (§4.6), with two corrections to how this
> box should be read.** The mechanism holds: on Maurano at p=23, register-equality
> drifts up to 0.124 as λ crosses the knee while IE stays at its unsubsampled
> floor — a 670× gap at subB=0.01. But (a) the underlying CSV for the 4×4 table
> below was never persisted, so it is not reproducible and is retained only as the
> original observation; and (b) its "IE is *lower* at subB=0.01 than at 0.1"
> reading is not supportable — that corpus is 15 pairs from 6 files whose
> between-pair J spread (sd 0.0020) equals the MAE being compared, so the
> 0.0022 → 0.0018 step is about one standard error. Cite §4.6, which has 190 pairs
> and jackknife CIs, for anything load-bearing.
> Every MAE above is measured on `jaccard_similarity`, which is
> register-equality and carries a chance-agreement floor `c` that is a step
> function of the load factor λ = n/m. Subsampling lowers `n`, so at low
> enough subB it drags λ through the knee at λ ≈ 5 and `c` — the estimator's
> own offset — moves. The MAE then reports the yardstick shifting, not the
> sample degrading.
>
> Measured directly on a 4×4 synthetic 10k corpus at p=18, reporting both
> estimators from the same runs:
>
> | subB | λ = n/m | `c` (theory) | MAE `jaccard_similarity` | MAE `jaccard_similarity_ie` |
> | ---- | ------- | ------------ | ------------------------ | --------------------------- |
> | 1.0  | 173.8   | 0.1699       | — (reference)            | — (reference)               |
> | 0.1  | 17.3    | 0.1699       | 0.0009                   | 0.0022                      |
> | 0.01 | 1.73    | 0.1525       | **0.0199**               | **0.0018**                  |
>
> The register-equality column reproduces the published ≈0.02 cliff. The
> inclusion-exclusion column does not move at all — it is *lower* at
> subB=0.01 than at subB=0.1. The set-Jaccard estimate is fine; only the
> offset moved. Mean `jaccard_similarity` falls 0.2690 → 0.2491 while mean
> `jaccard_similarity_ie` holds at 0.1042 → 0.1060.
>
> The λ mechanism accounts for the bulk of the shift but not all of it: the
> affine model predicts `(c(173.8) − c(1.73))·(1 − J̄) = 0.0156` against 0.0199
> measured, the gap being the model's known ±0.025 curvature. It also explains
> observation 3 without appealing to spatial structure — Maurano files cover
> 63–139 Mbp, so at p=18 their λ is 240–531 and subB=0.01 only brings them to
> 2.4–5.3, at or above the knee (`c` = 0.1644 at λ=2.4), a far smaller shift
> than the synthetic corpus's drop to λ=1.73.

The per-pair error box gives the distribution view, including the heavier
real-data tails:

![](figures/maurano_jaccard_error_vs_subB.png)

---

## Absolute comparison: hammock vs bedtools on Maurano

Bedtools wall (8-way GNU parallel, sort time not charged): **11.08 s
median**, very tight (10.99–11.10 s across 5 reps).

> The `bedtools / hammock` column below is **superseded**: it divides by an
> 11.08 s baseline measured with three processes per pair. The corrected
> baseline is **7.34 s**. Corrected values, mixed-stride, job 29652709:
> subB 1.00 → **0.90×** (1.12× slower), 0.10 → **1.64×**, 0.01 → **2.58×**.
> The hammock wall column is also ~1.16–1.28× stale (fused pairwise pass,
> reserved allocation). Left in place because the *within-hammock* method
> comparison it supports — mixed-stride beating hash-threshold and single-hash
> at equal subB — is a ratio among these rows and survives unchanged.

| method         | subB | hammock wall | bedtools / hammock |
|----------------|------|--------------|--------------------|
| (any)          | 1.00 | ~9.5 s       | **1.16×** (now 0.90×) |
| hash-threshold | 0.10 |  8.95 s      | 1.24×              |
| **mixed-stride** | **0.10** | **5.20 s** | **2.13×** (now 1.64×) |
| single-hash    | 0.10 |  9.21 s      | 1.20×              |
| mixed-stride   | 0.05 |  4.37 s      | 2.54×              |
| mixed-stride   | 0.01 |  3.64 s      | **3.04×**          |

The headline plot (top of this document) captures this in a single line
graph: bedtools is one fixed reference; only mixed-stride's wall time
drops meaningfully as subB shrinks.

The accuracy/speed tradeoff against bedtools is in the Pareto plot at the
top — at MAE < 10⁻³ mixed-stride sits at ~2× faster than bedtools; pushing
to MAE ≈ 2 × 10⁻³ buys 3×. Hash-threshold and single-hash never get more
than ~1.5× faster than bedtools even at the most aggressive subB tested.

---

## Recommendation

**`--subB-method mixed-stride --subB 0.1`** is the right default on this
evidence:

- Synthetic: 3–4× faster than no-subsample hammock, MAE ≤ 10⁻³ except the
  small-file (10k) cliff at very low subB — which the box above shows is a
  register-equality artifact; on `jaccard_similarity_ie` there is no cliff.
- Maurano (real): 1.83× faster than no-subsample hammock and **2.13×
  faster than bedtools**, MAE ~10⁻³, max error ~3 × 10⁻³.

Aggressive `subB=0.01` is worth considering only on large
(≥1M-interval), high-Jaccard corpora where MAE stays at ~10⁻⁵. On real,
moderate-sized data expect roughly a doubling of max-error per decade of
subB drop below 0.1. **Caveat, same root cause:** that guidance is calibrated
on `jaccard_similarity`, so part of what it is protecting against is the floor
moving rather than the sample degrading. If you read
`jaccard_similarity_ie`, the binding constraint at low subB is instead that λ
must stay high enough for the sketch itself to be informative. **Re-derived
2026-08-08 — see §4.6. On the IE column the guidance is not binding at all:
error stays within 1.66× its unsubsampled floor all the way down to subB=0.01,
at both p=18 and p=23.**

`hash-threshold` and `single-hash` are not competitive speed defaults on
either corpus — at best ~1.5× speedup vs no-subsample even at `subB=0.01`
(1.54× is the maximum over the `hash-threshold`/`single-hash` cells,
hash-threshold on synthetic 10k; Maurano tops out at 1.32×. For contrast,
mixed-stride reaches 14.0× on synthetic 1M at the same subB). Keep them for
parity reasons (`hash-threshold` matches the original hammock gate) or for
codepath simplicity (`single-hash`), not performance.

---

## Does subB affect `jaccard_similarity_ie`? (2026-08-08)

**No, not meaningfully.** Across 0.01 ≤ subB ≤ 1, all three methods, p ∈ {18, 23},
IE error against exact bedtools truth never exceeds **1.66×** its own
unsubsampled floor, and 26 of 30 subsampled cells sit within 1.3×.

Generators: `run_ie_subb.py` (32 cells, unmodified v0.7.0 `hammock-cpp`, ~7 min)
and `analyze_ie_subb.py` → `results/ie_maurano_summary.csv`,
`results/ie_maurano_pairs.csv`. This needed no changes to the sweep harness:
`hammock-cpp` emits `jaccard_similarity_ie` natively by default since v0.7.0, and
exact per-pair truth already existed in `docs/data/maurano_bedtools_ref.tsv`.
`--reps 1`, because the archived replicates are byte-identical (all 3420
(method, subB, pair) cells hold exactly one distinct Jaccard across the five).

Unlike everything above, the primary statistic here is **accuracy** — |IE −
bedtools| — not drift against hammock's own subB=1.0 run. IE is on bedtools'
scale; register-equality is not, so its numbers below remain *drift* and the two
must not be compared by level.

| p | unsubsampled IE floor | worst subsampled cell | λ at subB=1 → 0.01 |
|---|---|---|---|
| 18 | 1.152e-3 ± 1.0e-4 | 1.556e-3 (1.35×, mixed-stride @ 0.1) | 409 → 4.1 |
| 23 | 2.022e-4 ± 2.2e-5 | 3.358e-4 (1.66×, hash-threshold @ 0.01) | 12.8 → 0.13 |

Uncertainty is leave-one-file-out jackknife over the 20 files (the 190 pairs are
not independent — each file is in 19 of them).

**The effect is real but tiny.** 7 of 30 subsampled cells show a bias more than
3 SE from zero — largest is mixed-stride at p=18, subB=0.1 (−1.33e-3, 8.2 SE),
and hash-threshold at p=23, subB=0.01 (+2.63e-4, 4.7 SE, λ=0.13, i.e. the sketch
far sparser than its register array). So subsampling *is* detectable on IE; it is
simply small enough not to matter against a corpus whose true J spans
0.135–0.627. Note the default method, `mixed-stride`, is the worst of the three
for IE at p=18 subB=0.1 — it is chosen for speed, and that trade is unchanged.

**The contrast with register-equality is the real finding.** At p=23, as λ falls
through the knee, register-equality drifts from its own baseline by up to 0.124
while IE stays flat at its floor:

| subB | λ | IE error | RE drift | ratio |
|---|---|---|---|---|
| 0.5 | 6.39 | 1.95e-4 | 1.63e-4 | 1× |
| 0.25 | 3.20 | 1.89e-4 | 4.49e-3 | 24× |
| 0.1 | 1.28 | 1.83e-4 | 4.25e-2 | 232× |
| 0.05 | 0.64 | 1.72e-4 | 7.97e-2 | 462× |
| 0.01 | 0.13 | 1.86e-4 | 1.24e-1 | 670× |

(mixed-stride; the other two methods track it.) IE at every one of those cells is
at or **below** its unsubsampled floor of 2.022e-4. This is the λ mechanism the
10k-cliff box below infers from synthetic data, now measured on real data against
exact truth: what moves at low subB is the register-equality offset, not the
set-Jaccard estimate.

**Harness validation.** The new runs reproduce all 15 archived per-method
register-equality drift cells at p=18 to five significant figures (ratio 1.000
throughout) — so there is no pre-0.7.0 → v0.7.0 binary drift, and `--reps 1` loses
nothing.

**Scope.** Maurano only, Mode B, p ∈ {18, 23}, J ≥ 0.135. The low-J regime where
IE is censored at 0 and uninformative (J ≲ a few/√m) is **not** tested, and it is
where subsampling would be most likely to bite. A null here means "no effect above
a minimum detectable ≈ 4.1e-4 at p=18", not "no effect".

---

## Diagnostic figures

The headline figures above were generated by `headline_figures.R`. The
diagnostic figures from `analyze.R` cover the same data with per-method
breakdowns:

- `figures/<corpus>_wall_time_vs_subB.png` — raw wall times per (method, size).
- `figures/<corpus>_speedup_vs_nosub.png` — speedup vs subB=1.0 cross-method
  baseline (bedtools overlaid as a dashed reference on the maurano panel).
- `figures/<corpus>_speedup_within_method.png` — each method against its
  own subB=1.0 (useful for spotting the subB=0.5 hump).
- `figures/<corpus>_mae_vs_subB.png` and `<corpus>_jaccard_error_vs_subB.png`
  — accuracy curves.

Full per-cell numbers: `results/summary_<corpus>_<stamp>.csv`.
