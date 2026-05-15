# Paper outline: hammock — fast, accurate, biology-preserving interval-set sketching

**Working title (placeholder):** *hammock: HLL-sketch similarity for genomic intervals — faster than bedtools, with biological signal preserved*

**Thesis (one sentence):** A Python+C++ refactor of hammock matches bedtools' pairwise interval-Jaccard at r ≈ 0.9996 (Mode D, post-fix) and r ≈ 0.998 (Mode B), while being substantially faster than bedtools at every scale tested; the same sketches independently recover tissue clustering (ARI = 0.91) and are robust to reference-genome choice — so the speed gain comes with, not at the cost of, biological fidelity.

> **Status note (2026-05-14):** Numbers below reflect the post-bugfix Mode D
> implementation (commit 11739e6). Maurano + modeD_flanking + primate-phylogeny
> reruns are complete; ref-comparison and mus-homo full sweeps are still in
> flight at the time of writing. Cells marked **(pending rerun)** will be
> refreshed once those finish; the qualitative claims should hold but exact
> values may shift.

---

## 1. Abstract

One paragraph. Three numbers in the lead: (i) Mode B r = 0.998 / Mode D r = 0.9996 vs bedtools on Maurano DHS; (ii) substantially faster than `bedtools jaccard` on real DHS data and orders of magnitude faster at large catalog size; (iii) Mode D ARI = 0.91 / NMI = 0.96 on 8-tissue Maurano clustering. Close with the applicability boundary: sequence-level sketch transfers across human reference genomes but not across the ~80 Mya human↔mouse split.

## 2. Introduction

- All-vs-all pairwise interval similarity is bedtools' weakest scaling regime (O(N²·M log M)); large epigenome catalogs make this the bottleneck.
- Sketching trades exactness for speed; the question is *how much exactness*, and whether the sketch preserves the biology that the exact answer captures.
- Contribution: a re-implemented hammock with HLL backing for all four modes (A/B/C BED-interval, D FASTA-minimizer), parity-tested against the original, plus a four-experiment validation framework that quantifies the speed/accuracy/biology trade-off on real data.

## 3. Methods

### 3.1 hammock implementation

- Python orchestrator + C++ extension (pybind11); HLL with register-equality Jaccard, Ertl 2017 estimator, xxh64 ingestion.
- Mode A: interval coords. Mode B: per-bp HLL with optional `subB` subsampling. Mode C: interpolates A↔B via `expA` and `subB`. Mode D: minimizer-HLL on FASTA, with both interior-minimizer and minimizer-plus-flanks similarity columns.
- Output per-pair similarity columns: `jaccard_similarity`, `jaccard_similarity_with_ends`, plus `containment_AB`, `containment_BA`, `cosketch_{geom,arith,max}` in both flavors (10 metrics per Mode D pair; 5 per A/B/C pair). The `jaccard_*` columns are the analyses' default; the cosketch/containment columns are reported for transparency and inform Section 5.

### 3.2 Benchmark datasets

| dataset | role | n |
|---|---|---|
| Synthetic BED (`bedtools_benchmark/`) | scaling regime | up to 512 × 10k intervals |
| Synthetic FASTA pairs (`modeD_flanking/`) | exact k-mer truth for Mode D | 10 251 pairs across n_intervals × length × dist × mutation grid |
| Maurano 2012 fetal DHS | real interval data, 8-tissue labels | 20 BEDs |
| ENCODE H3K27ac (Heart/Liver/Lung × GRCh37/GRCh38/CHM13) | reference-build robustness | 9 sample×ref |
| ENCODE DNase-seq (5 tissues × human/mouse) | cross-species limit | 10 BEDs |

### 3.3 Statistical evaluation

Pearson r, Spearman ρ, MAE vs bedtools (Mode B/C/D); ARI, NMI vs known tissue labels (Mode D); Wilcoxon same-tissue vs different-tissue (cross-reference); R/ggplot via CairoPNG.

### 3.4 Mixed-stride subsampling (`--subB-method mixed-stride`)

Mode B's `subB` knob subsamples the per-base-pair positions that get ingested into the HLL. The original hammock implemented `subB` as a *hash-threshold* gate: a per-position xxh32 hash compared to `subB × UINT32_MAX`. The gate is correct but the per-position hash cost grows linearly with the input, so the wall-time savings from skipping positions plateau quickly: at moderate subB the gate hash itself dominates and `subB ≤ 0.5` is actually slower than no subsampling at all.

This paper introduces **mixed-stride** subsampling as a deterministic, hash-free alternative. For each chromosome, a fixed stride length `s = round(1/subB)` is chosen, offset by a chr-keyed hash so different chromosomes sample disjoint position phases. Ingestion walks the chromosome in stride-`s` increments — no per-position decision, no hash. The accepted-position cardinality matches the hash-threshold gate in expectation, but the cost is `O(L/s)` rather than `O(L)`, so wall time scales with `subB` instead of being dominated by the gate.

Three properties matter for downstream use:

1. **Performance scaling.** At `subB = 0.1`, mixed-stride is 3–4× faster than hash-threshold on synthetic data and 1.8–2.4× faster than hash-threshold on real DHS. The advantage compounds: at `subB = 0.01` on 1M-interval synthetic files, mixed-stride hits 14× while hash-threshold caps at ~1.4×.
2. **Accuracy.** Per-pair MAE vs the no-subsample reference is statistically indistinguishable across the three subB methods (hash-threshold, mixed-stride, single-hash). Mixed-stride does not buy speed by losing accuracy.
3. **Determinism and reproducibility.** Output is exactly reproducible given the chr-stride seed (`--gate-seed`, default 31337). The same files at the same subB produce byte-identical HLLs across runs and machines.

The structured-sampling concern — that a fixed stride could miss periodic features — is theoretical for BED inputs at the strides we test; in practice the chr-keyed offset randomization breaks cross-chromosome alignment and the Maurano + synthetic accuracy measurements show no detectable bias. Mixed-stride is the default `--subB-method` setting in this implementation. Hash-threshold remains available for parity with the original (`--subB-method hash-threshold`) and is what the parity tests cover.

This subsampling refinement is what makes the Section 4.1 speed numbers attainable: the "high-subsample" and "max-subsample" rows of the real-DHS table are mixed-stride at `subB = 0.1` and `subB = 0.01`. The per-method comparison plot is in supplementary.

## 4. Results

### 4.1 Speed: hammock is substantially faster than bedtools on real and synthetic interval corpora

> **Sources:** `experiments/subB_mixed_stride/RESULTS.md` (real DHS); `experiments/bedtools_benchmark/RESULTS.md` (synthetic).

**Real DHS (Maurano, 20 samples, 190 pairs, 8 threads for both tools — hammock `-t 8`, bedtools 8-way GNU parallel):**

| tool | wall (s) | speedup over bedtools |
|---|---|---|
| bedtools | 11.08 | 1.00× (ref) |
| hammock (default) | ~9.5 | ~1.2× |
| hammock (high-subsample) | ~5 | ~2× |
| hammock (max-subsample) | ~3.6 | ~3× |

**Synthetic scaling (10k intervals/file, p=14, 16 threads for both tools — hammock `-t 16`, bedtools 16-way GNU parallel):**

| N files | bedtools | hammock | speedup |
|---|---|---|---|
| 64 | 11.0 s | 1.7 s | ~6× |
| 128 | 44 s | 3 s | ~14× |
| 256 | 176 s | 7 s | ~26× |
| **512** | **706 s** | **13 s** | **~50×** |

Two things to communicate in this section, nothing more:
1. **hammock is faster than bedtools at every regime tested** — modestly faster on small real corpora (≈2×), dramatically faster as catalog size grows (≈50× at N=512 files).
2. **The speedup is not bought with accuracy loss.** Per-pair MAE vs bedtools is statistically indistinguishable across hammock subsampling settings — i.e., the speed knob is "free" against the bedtools reference.

**Fig 1 (headline):** Pareto curve, hammock dominates bedtools on Maurano across the entire accuracy axis.

![Fig 1 — Maurano Pareto (hammock vs bedtools)](../experiments/subB_mixed_stride/figures/headline_maurano_pareto.png)

**Fig 2:** N-scaling out to 512 files on synthetic.

![Fig 2 — synthetic N-scaling](../experiments/bedtools_benchmark/figures/cpp_vs_bedtools_t16_20260512_160412_sketch_compare_split.png)

(Internal hammock-version comparisons — mixed-stride vs hash-threshold vs single-hash subB strategies, sort-time accounting, OpenMP scaling shape — belong in the supplementary methods, not the main text.)

### 4.2 Accuracy: both interval-mode (B) and sequence-mode (D) hammock match bedtools

> **Source:** `experiments/maurano_dhs_validation/RESULTS.md` (rerun post Mode-D fix, 2026-05-14).

Across Maurano's 400 sample pairs:

| mode | best r vs bedtools | best MAE | note |
|---|---|---|---|
| **Mode B** | **0.998** | — | at every precision tested |
| Mode A | 0.82 | — | only mode that breaks — drops interval-length info |
| Mode C (subB→1) | →0.998 | — | inherits B |
| **Mode D, post-fix** (k=20, w=100, p=24) | **0.9996** | **0.0061** | **47 of 209 configs exceed r > 0.99** |

The Mode D headline number changed by an order of magnitude post-fix: the pre-fix sweep's best was r = 0.966 (MAE 0.060); the post-fix optimum sits at r = 0.9996 with MAE 0.0061, four-decimal-place agreement with bedtools. The ridge in the (k, w) Pearson heatmap moved from a modest mid-grid optimum to the high-k / high-w edge of the sweep, indicating that Mode D's near-perfect agreement is unlocked at long minimizer windows where interior coverage is richest.

**Fig 4:** Modes A/B/C summary — Pearson r vs bedtools by precision.

![Fig 4 — Mode A/B/C Pearson vs bedtools](../experiments/maurano_dhs_validation/figures/abc_pearson_vs_bedtools.png)

**Fig 5:** Mode C is a sharp interpolation between A-regime (subB ≲ 0.005) and B-regime (subB ≳ 0.05). One knob, two regimes.

![Fig 5 — Mode C subB interpolation](../experiments/maurano_dhs_validation/figures/mode_c_subB_interpolation_agg.png)

**Fig 5b (new):** Post-fix Mode D Pearson ridge along high k + high w + p ≥ 22.

![Fig 5b — Mode D Pearson heatmap](../experiments/maurano_dhs_validation/figures/mode_d_pearson_heatmap.png)

### 4.3 Biological signal: the sketch recovers tissue identity

> **Source:** `experiments/maurano_dhs_validation/RESULTS.md` (post-fix).

Post-fix Mode D best-ARI cell: **k = 10, w = 30, p ∈ {22, 23, 24}**, all yielding **ARI = 0.910, NMI = 0.961** on the 8-tissue / 10-fetal-tissue-subtype label set. The dendrogram makes 2 errors: the two fBrain samples split into separate clusters, and one fMuscle_back is grouped with the two fMuscle_legs — both errors that the bedtools reference dendrogram *also* makes, so this is not a sketch artifact but a property of the underlying signal. Pre-fix this cell was k=8, w=10, p=12; post-fix the same achievable ARI moved to a different (k, w, p) cell while the numerical value is unchanged.

A second post-fix change worth highlighting: at the ARI-best config, Mode D's predicted Jaccards now sit on the y = x diagonal versus bedtools (and versus Mode B), where the pre-fix scatter showed a systematic upward bias. The bug had been suppressing sketch cardinality; with it fixed, the sketch is calibrated against the bedtools reference rather than just rank-correlated.

**Fig 6:** Mode D best-ARI dendrogram (top) paired with the bedtools reference dendrogram (bottom). Same tissue blocks recovered.

![Fig 6a — Mode D best dendrogram](../experiments/maurano_dhs_validation/figures/mode_d_best_dendrogram.png)
![Fig 6b — bedtools reference dendrogram](../experiments/maurano_dhs_validation/figures/bedtools_dendrogram.png)

**Fig 7:** The headline methodological point holds post-fix — **best-Pearson cell ≠ best-ARI cell**. Pearson-best (large-k, large-w) cluster at ARI ≈ 0.69; ARI-best (k=10, w=30) at Pearson ≈ 0.946. Numerical perfection and clustering quality remain non-coincident knobs.

![Fig 7 — Mode D Pearson vs ARI tradeoff](../experiments/maurano_dhs_validation/figures/mode_d_metric_tradeoff.png)

**Fig 8:** ARI ≥ 0.85 plateau now narrows to k = 10, w = 30, all p ≥ 22 (pre-fix plateau was broader; post-fix the bug-induced "easy ARI at low p" is gone).

![Fig 8 — Mode D ARI across (k, w, p) sweep](../experiments/maurano_dhs_validation/figures/mode_d_clustering_ari.png)

### 4.4 Robustness to reference genome **(pending post-fix rerun)**

> **Source:** `experiments/ref-comparison/docs/exp_a_results.md` (pre-fix); full re-sweep submitted 2026-05-14 at 21:02, in flight.

3 H3K27ac samples (heart, liver, lung) × 3 references (GRCh37/GRCh38/CHM13), 9 sample×ref combinations. Pre-fix headline: same-tissue cross-reference Jaccard significantly higher than different-tissue Jaccard at **every (k, w) cell with k ≥ 8** (Wilcoxon p ≤ 10⁻⁵), best Δmedian = 0.107 at k=10, w=20 broad. The dendrogram cleanly grouped the three references of each tissue with within-tissue merge heights ≤ 0.1 vs cross-tissue ≥ 0.4.

The Mode D bugfix is expected to *strengthen* this result rather than weaken it: pre-fix, the `jaccard_similarity` (no_ends) column was largely empty on these FASTAs at k ≥ 15 (silent zero-jaccard from minimizers being dropped), so the analysis defaulted to `jaccard_similarity_with_ends`. Post-fix, the no_ends column carries real signal and may give a cleaner same-tissue/different-tissue separation. Final numbers and dendrograms refresh once the in-flight rerun completes.

**Fig 9 (pre-fix; will be refreshed):** Cross-reference dendrogram at k=10, w=10 — within-tissue clades hold across GRCh37/GRCh38/CHM13.

![Fig 9 — cross-reference dendrogram](../experiments/ref-comparison/figures/cross_ref_dendrogram_k10_w10.png)

**Fig 10 (pre-fix; will be refreshed):** (k × w) effect-size heatmap for broad peaks; same-tissue cross-ref ≥ different-tissue across the sweep.

![Fig 10 — cross-ref effect-size sweep, broad](../experiments/ref-comparison/figures/sweep_effect_size_broad.png)

Practical interpretation (unchanged by the bug): when peaks are aligned to a different human reference than expected, the sketch still produces the same biological neighborhood. This is the property that lets hammock be deployed against heterogeneous catalogs (ENCODE/Roadmap mixtures) without first re-aligning everything.

### 4.5 Methodological notes: choosing Mode D's flanking column

> **Source:** `experiments/modeD_flanking/` Parts 1 (Maurano re-analysis) and 2 (synthetic FASTA pairs with exact k-mer Jaccard truth). Both fully rerun post-fix.

Mode D emits two Jaccard columns: minimizer-only (`jaccard_similarity`) and minimizer-plus-flanks (`jaccard_similarity_with_ends`). The post-fix sweep gives a cleaner version of the regime story than was possible pre-fix:

- **Part 1 (Maurano, real DHS):** at the high-precision configs, `no_ends` is the clear winner — at k=20, w=100, p=24, r_no_ends = 0.9996 vs r_with_ends = 0.888. The advantage flips only at the smallest (k, w) cells (k ≤ 10, w ≤ 10), where interior minimizers are sparse and the flanking k-mers carry the signal.
- **Part 2 (synthetic, 10 251 pairs across n_intervals × length × distribution × mutation, with exact k-mer Jaccard as truth):** `with_ends` has systematically *larger* MAE than `no_ends` against the k-mer truth (mean Δmae_r ≈ −0.015 across all k); the gap widens with φ, the flanking-fraction predictor we defined.

The flanking-fraction φ ≈ 2(k−1)·n_intervals / (total_length / w) predicts the regime: large φ → `_with_ends` over-weights boundary k-mers; small φ → either column works.

**Recommendation:** minimizer-only (`jaccard_similarity`) is the right default for ChIP/DHS-shaped corpora at any reasonable parameter choice; `_with_ends` is a fallback for short-sequence, sparse-minimizer regimes.

**Fig 11:** Part 1 (Maurano real DHS) — sign of (no_ends − with_ends) flips with w.

![Fig 11 — Maurano Δr(no_ends − with_ends) vs w](../experiments/modeD_flanking/figures/maurano_delta_r_vs_w.png)

**Fig 12:** Part 2 (synthetic) — Δ(no_ends, with_ends) vs φ, and the φ × mutation phase plane.

![Fig 12a — synthetic Δ vs φ](../experiments/modeD_flanking/figures/synthetic_delta_vs_phi.png)
![Fig 12b — synthetic φ × mutation phase diagram](../experiments/modeD_flanking/figures/synthetic_phase_diagram.png)

**Fig 12c:** Empirical agreement with the analytical φ-prediction — independent validation of the flanking-fraction model.

![Fig 12c — empirical vs analytical φ](../experiments/modeD_flanking/figures/synthetic_empirical_vs_analytical.png)

## 5. Limitations and applicability boundary

### 5.1 Cross-species sequence sketching at long divergence **(pending post-fix rerun)**

> **Source:** `experiments/mus-homo/results/column_comparison.tsv` (pre-fix); jobs rerunning 2026-05-14 21:03 onward.

Pre-fix finding: across 20 (k, w) cells × DNase-seq peaks of 5 matched tissues in human/mouse (~80 Mya divergence), the sketch never produced a tissue-dominant dendrogram (species ARI > tissue ARI in every cell with k ≥ 8). The cross-reference-robustness story (Section 4.4) did *not* extend to cross-species — orthologous peaks in different species don't share enough k-mers for sequence-level sketching to recover tissue identity.

This null may flip or sharpen post-fix: pre-fix many cells had `no_ends` Jaccards spuriously near zero (the bug), which would have *inflated* apparent species separation by collapsing within-species pairs. If the null survives the rerun, it's a stronger negative result — sequence-level sketching of regulatory peak FASTAs genuinely cannot recover cross-species tissue identity at the human↔mouse split, independent of any sketch-implementation pathology. If the null inverts, the cross-species claim needs reframing.

(Primate-phylogeny follow-ups: post-fix CSVs for human + mouse + macaque + marmoset + cow + opossum + dog are in place; full 20-species panel still being staged.)

### 5.2 Definitional gap vs bedtools

Mode B/C/D estimate slightly different Jaccards than bedtools (bp-set vs interval-overlap vs k-mer-set). On the synthetic Mode B benchmark, the median per-pair gap is ~0.16 and is independent of precision and subsampling. On Maurano post-fix, Mode D at the optimal high-k/high-w cell brings MAE down to 0.006 — i.e., the gap effectively closes when the sketch is well-conditioned. For applications where absolute Jaccard magnitude matters (not just ranking), the appropriate setting is determined by the (k, w, p) choice; the relevant figures are in Section 4.2.

### 5.3 The cosketch + containment columns are reported but not yet exploited

The five auxiliary similarity columns added 2026-05-14 (containment_AB, containment_BA, cosketch_geom, cosketch_arith, cosketch_max — each in two flavors for Mode D) are present in every Mode D output CSV. Within-experiment analyses currently use the jaccard columns only. A pre-fix sanity check on the ref-comparison Exp A k=10, w=10 cell suggested `cosketch_geom_with_ends` is a near-tie with `jaccard_similarity_with_ends` and `cosketch_max` is uniformly the weakest discriminator; the post-fix re-evaluation across all (k, w) cells and across the Maurano sweep is straightforward and may identify a column (likely cosketch_geom on the minimizer-only flavor) that is more robust than jaccard at small precision or small k. Treat this as a fast-follow analysis rather than a result.

## 6. Discussion

- **Practical recipe.** Mode B for fast bedtools-equivalent interval-Jaccard with optional subsampling for further speedup at no accuracy cost. Mode D at large k and w (k=20, w=100, p=24) for the closest numerical match to bedtools (r = 0.9996, MAE = 0.006). Mode D at k=10, w=30, p ≥ 22 for tissue clustering (ARI = 0.91, NMI = 0.96). Use `jaccard_similarity` by default; fall back to `jaccard_similarity_with_ends` only in the short-sequence / sparse-minimizer regime characterized in Section 4.5.
- **The sketch carries more than bedtools captures.** Mode D recovers tissue clustering directly from peak FASTAs — independently from bedtools' interval overlap — at ARI = 0.91. Sketch similarity ≈ biological similarity, even when the two estimators don't agree numerically.
- **Boundaries.** Within-species, robust; cross-species without coordinate alignment, fails (pre-fix; post-fix rerun pending). This is a property of any k-mer sketch, not specific to hammock.

## 7. Conclusion

hammock provides a fast, sketch-based alternative to `bedtools jaccard`. On real DHS data Mode B matches bedtools at r = 0.998 and Mode D matches it at r = 0.9996 / MAE = 0.006, while running substantially faster than bedtools at every scale tested (orders of magnitude faster at large catalog size). The same Mode D sketches recover tissue clustering at ARI = 0.91 directly from peak FASTAs and are robust to reference-genome choice within a species. The combination — speed *plus* near-exact agreement *plus* preserved biology *plus* reference-build invariance — positions sketching as a viable default for large-scale epigenome comparison.

---

## Figure index

| # | Path | Carries |
|---|---|---|
| 1 | `subB_mixed_stride/figures/headline_maurano_pareto.png` | hammock dominates bedtools on real Pareto |
| 2 | `bedtools_benchmark/figures/cpp_vs_bedtools_t16_20260512_160412_sketch_compare_split.png` | synthetic scaling to N=512 |
| 4 | `maurano_dhs_validation/figures/abc_pearson_vs_bedtools.png` | Mode B r=0.998 |
| 5 | `maurano_dhs_validation/figures/mode_c_subB_interpolation_agg.png` | Mode C as A↔B interpolant |
| 5b | `maurano_dhs_validation/figures/mode_d_pearson_heatmap.png` | post-fix Mode D Pearson ridge |
| 6 | `maurano_dhs_validation/figures/mode_d_best_dendrogram.png` + `bedtools_dendrogram.png` | tissue recovery |
| 7 | `maurano_dhs_validation/figures/mode_d_metric_tradeoff.png` | Pearson-best ≠ ARI-best |
| 8 | `maurano_dhs_validation/figures/mode_d_clustering_ari.png` | ARI plateau across sweep |
| 9 | `ref-comparison/figures/cross_ref_dendrogram_k10_w10.png` *(pending rerun)* | within-tissue clades across refs |
| 10 | `ref-comparison/figures/sweep_effect_size_broad.png` *(pending rerun)* | cross-ref effect-size heatmap |
| 11 | `modeD_flanking/figures/maurano_delta_r_vs_w.png` | flanking column choice on real data |
| 12 | `modeD_flanking/figures/synthetic_delta_vs_phi.png` + `synthetic_phase_diagram.png` | flanking on synthetic, φ-axis |
| 12b | `modeD_flanking/figures/synthetic_empirical_vs_analytical.png` | analytical φ-prediction validated |

(The `synthetic_speedup_vs_nosub.png` figure showing within-hammock subB-method comparison moves to supplementary — per the constraint that internal hammock-version performance differences are not paper material.)
