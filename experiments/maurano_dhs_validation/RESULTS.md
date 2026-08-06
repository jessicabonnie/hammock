# maurano_dhs_validation — results notes

Generated from the 2026-05-14 sweep (209 Mode D configs × 2 columns ×
2 references + 24 A/B/C configs). Numbers below come from
`results/abc_summary.csv` and `results/mode_d_summary.csv`; figures live
in `figures/`.

> **2026-05-14 rerun after Mode D bug fix.** The earlier sweep's Mode D
> `jaccard_similarity` (no_ends) column was contaminated by a refactor
> bug in `MinimizerSketch.add_string` that silently dropped most
> minimizer hashes (see `CLAUDE.md` "Intentional divergences" §6 and
> `memory/project_modeD_no_ends_zero_bug.md`). Pre-fix numbers
> understated `no_ends` performance dramatically — best Pearson r jumped
> from 0.966 to **0.9996** at the post-fix optimum. All Mode D numbers
> below are post-fix. Pre-fix CSVs are kept at `results/raw_d_buggy_pre_fix/`.

> **Audit note, 2026-08-06 — the summary CSV has grown since this was
> written; the headline numbers have not moved.** `mode_d_summary.csv` now
> holds **235** Mode D configs (not 209) × **3** columns (`jaccard_similarity`,
> `jaccard_similarity_with_ends`, and `jaccard_similarity_ie`, added in
> v0.5.0) × 2 references, and `abc_summary.csv` holds **29** A/B/C configs
> (not 24 — the Mode C sweep was extended to expA = 4.0 and down to
> subB = 1e-4). Re-counting the two census claims below against the current
> file: `r > 0.99` is still **47** configs, and `MAE < 0.05` is now **29**
> (was 21), the increase coming from the high-k × large-w fill-in run
> (`sbatch_fill_highk_w.sh`), not from any change to the numbers. Every
> named optimum below (Pearson 0.9996 at k=20/25 w=100 p=24 — w=50 ties;
> MAE 0.0061 at k=25 w=100 p=24; ARI 0.9102 / NMI 0.9610 at k=10 w=30 for
> every p ≥ 12; with_ends topping out at r = 0.9740 at k=20 w=20 p=20)
> re-verifies exactly against the current CSV.
>
> **Not analysed below: `jaccard_similarity_ie` is the stronger column.** On
> the same corpus it reaches Pearson r = 1.0000 and MAE = 0.0019 (k=20,
> w=20–30, p=24) against bedtools, versus 0.9996 / 0.0061 for the
> register-equality column this document reports. That is the expected
> direction — `jaccard_similarity` is register-equality and is not on
> bedtools' scale (`CLAUDE.md` divergence #2) — but the prose below was
> written before the column existed and has not been re-derived on it.
>
> **`jaccard_similarity_with_ends` no longer exists.** It was removed in
> v0.6.0 (`CLAUDE.md` divergence #8). The `with_ends` results below are
> retained as the record of why; they cannot be regenerated from a current
> build, and the archived `results/raw_d/` CSVs are the only source.

---

## Modes A, B, C vs bedtools

![](figures/abc_pearson_vs_bedtools.png)

**Mode B tracks bedtools closely (r = 0.998 at every precision), and every
Mode C variant inherits that.** Mode A is the only mode that breaks down
(r ≈ 0.82) — the interval-only HLL throws away the length information that
bedtools cares about. So when we later ask "does Mode D recapitulate
bedtools," Mode B is a serviceable stand-in.

> **Two caveats on both numbers, added 2026-07-31.** (i) These `r` values are
> computed over the full square matrix, which includes 20 self-pairs pinned at
> (1,1) — pure leverage. Off-diagonal, Mode B is r = 0.9972 and Mode A is
> r = **0.069**, i.e. Mode A is not "degraded", it is uninformative here and
> its results should not be quoted. (ii) `jaccard_similarity` is
> register-equality and is not on bedtools' scale, so a high `r` means the two
> track each other affinely, not that the values agree; Kendall τ for Mode B is
> 0.951, so 2.45% of comparisons invert. "Essentially perfect" overstated both.
> See `docs/jaccard-definitional-gap.md`.

---

## Mode C interpolation between Mode A and Mode B

![](figures/mode_c_subB_interpolation_agg.png)

**subB is a sharp knob between Mode-A-like (low subB) and Mode-B-like
(high subB), with the crossover near subB ≈ 0.005.** Below that Mode C
closely tracks Mode A; by subB = 0.01 it has flipped to Mode-B-like
(|C − B| = 0.14 vs |C − A| = 0.29). Above subB ≈ 0.05 Mode C and Mode B
are indistinguishable at the per-pair level.

![](figures/mode_c_expA_interpolation_agg.png)

**expA is a weaker but still real knob — the bend only shows up past
expA ≈ 2, with the A/B crossover sitting between expA = 3.0 and 3.5.**
Below expA = 2, |C − Mode B| stays under 0.05 (Mode-B regime). By
expA = 4, |C − Mode A| has dropped to 0.17 while |C − Mode B| has climbed
to 0.26 (Mode-A regime). Compared to subB, the transition is much
shallower — you need to push expA all the way to ~4 to see the same
swing that subB delivers between 0.001 and 0.005.

---

## Mode D vs bedtools — sweep maps

![](figures/mode_d_pearson_heatmap.png)

**The `no_ends` column achieves near-perfect agreement at high k + high w
+ high precision.** 47 of 209 configs exceed Pearson r > 0.99 against
bedtools, peaking at **r = 0.9996** at k = 20/25, w = 100, p = 24. The
ridge in the heatmap runs along the right edge — large k and large w —
not the modest (k=10, w=20) ridge the pre-fix sweep had shown. The
`with_ends` column tops out around r = 0.97 at k = 20, w = 20.

![](figures/mode_d_mae_heatmap.png)

**Numerical agreement is excellent at the Pearson optimum.** The MAE-best
config is k = 25, w = 100, p = 24 with **MAE = 0.0061** — essentially
indistinguishable from bedtools. 21 of 209 cells have MAE < 0.05.
Pearson-best and MAE-best are now in the same region of the (k, w) plane
(both at high k + high w + p = 24), unlike the pre-fix story where they
disagreed.

---

## Pearson vs Spearman per config (does scale matter?)

![](figures/mode_d_pearson_vs_spearman_scatter.png)

**Most configs sit just below the y = x line: Spearman is consistently
0.05–0.10 below Pearson.** That gap means the bedtools→Mode-D relationship
isn't *purely* monotone — there's some rank shuffling, not just a scale
offset. A few low-precision configs that have decent Pearson collapse to
much worse Spearman (the bottom-left cluster), which is the more honest
signal that they aren't reliably ordering pairs.

---

## Mode D vs bedtools vs Mode B (per config)

![](figures/mode_d_bedtools_vs_modeB_scatter.png)

**Points lie tightly on y = x: r(D, bedtools) ≈ r(D, Mode B) for every
config.** That's expected — Mode B correlates with bedtools at r = 0.998
(0.9972 off-diagonal), so it is a good *correlation* reference. Note the
weaker claim: it is interchangeable for questions about covariation, not for
questions about magnitude, since Mode B's `jaccard_similarity` sits a floor
`c·(1 − J)` above bedtools. This is more useful as a sanity check than a
discovery.

---

## The headline tradeoff — numerical agreement vs clustering recovery

![](figures/mode_d_metric_tradeoff.png)

**Pearson-best and ARI-best are still in different corners.** The
Pearson-best configs (large k + large w + p = 24) cluster far right
(Pearson ≈ 0.9996) but mid-vertical (ARI ≈ 0.69). The ARI-best config
sits high-left (Pearson ≈ 0.946 but ARI = 0.910). Numerical perfection
doesn't buy you the best tissue clustering — different metrics still
pick different optima.

---

## Clustering recovery — the actual tissue clustering

The ARI-best Mode D config is **k = 10, w = 30, p ≥ 12** with
**ARI = 0.910, NMI = 0.961** against the 10 fetal-tissue labels. (Pre-fix
this peaked at k = 8, w = 10, p = 12 with the same ARI/NMI — the
clustering signal moved to a different cell of the grid but the achievable
quality is identical.)

![](figures/mode_d_best_dendrogram.png)

**The dendrogram is annotated with the predicted 10 clusters (blue
boxes).** Most tissue groups end up as their own clade: the 5
fIntestine_Sm samples, the 3 fHeart, the 3 fLung, the kidney/skin/stomach
singletons. The two errors visible are (a) the two fBrain samples land
in different clusters and (b) one fMuscle_back is grouped with the two
fMuscle_legs — a sub-tissue split that arguably shouldn't be penalised.
Compare to the bedtools reference below.

![](figures/bedtools_dendrogram.png)

**Bedtools' own dendrogram cuts the same way for most tissues** — the
intestine block, the heart block, the muscle clade — modulo the same
brain-pair split and the same muscle_back/leg lumping. Mode D isn't
ignoring biology to score 0.91 on ARI; it's reproducing the bedtools
dendrogram structure.

### Cluster contingency

Predicted cluster × true tissue label at the post-fix ARI-best config
(`results/best_cluster_contingency.csv`):

| predicted | fBrain | fHeart | fIntestine_Sm | fKidney | fLung | fMuscle_arm | fMuscle_back | fMuscle_leg | fSkin | fStomach |
|----------:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | 1 | . | . | . | . | . | . | . | . | . |
| 2 | 1 | . | . | . | . | . | . | . | . | . |
| 3 | . | 3 | . | . | . | . | . | . | . | . |
| 4 | . | . | 5 | . | . | . | . | . | . | . |
| 5 | . | . | . | 1 | . | . | . | . | . | . |
| 6 | . | . | . | . | 3 | . | . | . | . | . |
| 7 | . | . | . | . | . | 1 | . | . | . | . |
| 8 | . | . | . | . | . | . | 1 | 2 | . | . |
| 9 | . | . | . | . | . | . | . | . | 1 | . |
| 10 | . | . | . | . | . | . | . | . | . | 1 |

Same shape as pre-fix: brain-pair split + muscle_back/leg lumping;
everything else clean. Per-sample assignments are in
`results/best_cluster_assignment.csv`.

![](figures/mode_d_cluster_confusion.png)

**Same contingency as a heatmap.** Every column has its mass concentrated
on a single row (or two adjacent rows for the muscle subtypes / brain
pair). That's why ARI gets to 0.91.

---

## Sweep maps for the clustering metrics

![](figures/mode_d_clustering_ari.png)

**ARI peaks in a different region from Pearson** — clustering is
preserved over k ∈ {10}, w ∈ {30}, **all precisions p ≥ 12**, with smaller
plateaus at neighbouring (k, w). `with_ends` (right column) is worse at
clustering because it inflates pair similarities disproportionately for
files with many short sequences, blurring the between-tissue contrast.

![](figures/mode_d_clustering_nmi.png)

**NMI tells the same story** — peak 0.96 at the same ARI-best config.

---

## Best-config scatter vs both references

![](figures/mode_d_best_scatter.png)

**At the ARI-best config, predictions sit near the y = x diagonal** —
the post-fix `no_ends` Jaccard is essentially calibrated against both
bedtools (left) and Mode B (right). The pre-fix systematic upward bias
is gone: the bug had been *suppressing* the sketch cardinality and now
the recovered Jaccards agree numerically.

---

## Summary — which metric for which question?

| If you want to answer... | Use | At config | Value |
|---|---|---|---|
| Do predictions covary with bedtools?   | Pearson r | k=20/25, w=100, p=24 | **0.9996** |
| Do predictions rank-order the same?    | Spearman ρ | k=20/25, w=100, p=24 | **0.998** |
| (interval mode's rank fidelity, for contrast) | Kendall τ | Mode B, p=21 | 0.951 (2.45% inverted) |
| Are the absolute Jaccards close?       | MAE        | k=25, w=100, p=24 | **0.0061** |
| Is the biology preserved?              | ARI        | k=10, w=30, p≥12 | **0.910** |
| Is the biology preserved (info-theory)?| NMI        | k=10, w=30, p≥12 | **0.961** |

**Post-fix headline: at high k + high w + high precision, Mode D's
`no_ends` column is a near-perfect drop-in replacement for bedtools
pairwise Jaccard** — 47 of 209 configs exceed r > 0.99 against bedtools,
with the best four-decimal-place agreement at MAE = 0.0061. The
clustering optimum sits at a *different* (smaller-k, smaller-w) cell of
the grid; numerical perfection and clustering quality remain
non-coincident knobs.
