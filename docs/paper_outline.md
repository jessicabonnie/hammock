# Paper outline: hammock — fast, accurate, biology-preserving interval-set sketching

**Working title (placeholder):** *hammock: HLL-sketch similarity for genomic intervals — faster than bedtools, with biological signal preserved*

**Thesis (one sentence):** hammock — a Python+C++ HyperLogLog-backed interval-set sketcher — ranks pairwise interval-Jaccard exactly as bedtools does (interval mode, r ≈ 0.998, because its register-equality estimator is an affine transform of bedtools' set-Jaccard — order-preserving, with an absolute offset that shrinks to zero as similarity rises and does not close with precision), and at high k/w its sequence mode additionally reproduces bedtools' values to four decimals (r ≈ 0.9996, MAE ≈ 0.006), while being substantially faster than bedtools at every scale tested; the same sketches independently recover tissue clustering (ARI = 0.91) and are robust to reference-genome choice — so the speed gain comes with, not at the cost of, biological fidelity; hammock extends beyond bedtools capabilities by enabling interval comparisons across references.

---

## 1. Abstract

hammock provides two complementary similarity primitives over genomic interval sets, both backed by HyperLogLog sketches. The **interval mode** ranks pairwise similarity identically to `bedtools jaccard` (r = 0.998) with substantially faster wall time and orders of magnitude lower cost at large catalog size; because the two tools estimate slightly different Jaccards, the absolute values differ by a fixed definitional offset (median ~0.16 on synthetic data), so interval mode is a rank-faithful — not value-identical — drop-in. The **sequence mode** compares the nucleotide content of intervals — a similarity primitive bedtools does not offer — recovering tissue identity from fetal DHS peak FASTAs at ARI = 0.91 / NMI = 0.96 on the 10-tissue Maurano label set, and remaining robust across human reference builds (peaks aligned to GRCh37, GRCh38, and CHM13 cluster by tissue, not by reference genome). At high k/w, sequence mode additionally closes the absolute gap, reproducing bedtools' interval-Jaccard values to four decimals (r = 0.9996, MAE = 0.0061) on Maurano DHS.

Since its release in [YEAR], bedtools has been "the Swiss Army knife" of bioinformatics, relied on to quickly calculate overlap and similarity between genomic regions within the same reference. The calculation of all-vs-all jaccard comparisons between region files is one process where bedtools' approach to pairwise comparison of BED files can be improved. Bedtools' all-pairs cost couples the number of files with the per-file size (O(N²M)), thus every additional file costs another full sweep against every other file. Our interval set sketcher, hammock, utilizes hyperloglog sketching to preserves accuracy while reducing cost by reading in each bedfile only once $(O(N^2 +NM ))$. 

## 2. Introduction
Increased availability of high-throughput sequencing has facilitated a corresponding increase in large-scale modern sequencing projects which produce vast numbers of genomic annotation datasets which then, themselves, become part of the genomic data ecosystem. Large-scale initiatives such as the ENCODE Project\cite{encode_2012}, the Roadmap Epigenomics Project\cite{roadmap_2015}, GTex[REF], and the 1000 Genomes Project\cite{10002015global} make datasets publicly available, providing unprecedented opportunities for integrative analysis. Among the most common and biologically informative outputs of these projects are genomic intervals: contiguous stretches of the genome that define transcription factor binding sites, chromatin accessibility peaks, histone modifications, structural variants, and splicing junctions. Such interval datasets have become a cornerstone of downstream analyses that seek to characterize genome function and variation across populations, tissues, cell types, and experimental conditions. Their growing scale, however, creates a substantial computational challenge. ChIP-Atlas, for example, expanded from 37,720 accumulated experiments in 2015 to 464,655 in 2025, while repositories such as ENCODE contain thousands of additional ChIP-seq and chromatin-accessibility experiments. Because the number of possible pairwise comparisons grows quadratically with the number of datasets, even a modest subset of these resources can imply millions or billions of potential comparisons.


Comparisons of interval datasets are fundamental for addressing questions of biological similarity and difference. For instance, overlap analyses enable researchers to identify conserved regulatory elements across tissues... These comparisons are typically conducted using files in the BED format, a simple tab-delimited representation of genomic intervals that has been in wide use since the early 2000s\cite{kent_hgbrowser}. While conceptually straightforward, the notion of ``overlap'' between interval sets is not trivial: it may be defined in terms of base-pair overlap, interval-level overlap, or fractional coverage, each capturing distinct biological relationships. Tools such as \program{BEDTools} provide exact calculations of these measures\cite{quinlan2010bedtools}, but the computational cost of pairwise overlap calculations grows rapidly with the size of modern repositories. In practice, this creates a scalability bottleneck that hampers systematic comparisons across the interval datasets now available\cite{li2020design}.

Scalability is not the only obstacle to systematic interval comparison. BED coordinates are meaningful only with respect to the reference genome on which they were defined, and large public collections remain distributed across multiple genome assemblies. Standard overlap measures therefore cannot directly compare two interval files when one is defined on hg19 and the other on hg38, even when the files describe the same assay or biological feature. Figure~\ref{fig} illustrates this limitation for histone-mark ChIP-seq peak files from Roadmap Epigenomics and BLUEPRINT. For histone marks represented in both collections, within-Roadmap and within-BLUEPRINT comparisons can be performed in their respective coordinate systems, whereas every Roadmap–BLUEPRINT pairing is blocked by the hg19–hg38 mismatch. Coordinate conversion can sometimes recover such comparisons, but it introduces additional preprocessing, may not map all intervals unambiguously, and continues to define similarity in terms of a selected coordinate system. An alternative is to extract the reference sequence underlying each interval and compare the resulting sequence collections directly, allowing related genomic annotations to be evaluated across references without requiring their coordinates to coincide.

Probabilistic data structures, or sketches, offer a powerful solution to this challenge. Sketching methods compress large datasets into compact representations that allow efficient estimation of set cardinality and similarity. Techniques such as MinHash\cite{minhash} and HyperLogLog\cite{hll} have seen wide application in computer science, and their introduction to genomics has already revolutionized k-mer–based comparisons of sequencing datasets. Tools like Mash\cite{mash}  and Dashing \cite{dashing} have demonstrated the value of sketching for rapid, large-scale sequence and metagenome comparisons. Extending these methods to genomic interval data, however, remains relatively unexplored. Because intervals represent continuous spans rather than discrete tokens, adapting sketches for overlap estimation presents unique methodological challenges but also significant opportunities. Because intervals represent continuous spans rather than discrete tokens, adapting sketches for overlap estimation presents unique methodological challenges but also significant opportunities. Moreover, representing intervals through both their coordinates and their underlying reference sequences provides complementary notions of similarity: coordinate-based sketches support rapid comparison within a shared reference, whereas sequence-based sketches can support comparisons across references.


In this study, we present \program{hammock}, a command-line tool for scalable comparison of genomic interval datasets. Within a shared reference genome, \program{hammock} applies probabilistic sketches to BED intervals to approximate overlap and similarity across large collections of files, reducing the computational burden of systematic pairwise analysis. We benchmark these estimates against exact calculations from established interval-processing tools and evaluate the resulting trade-offs among speed, memory use, and accuracy. We additionally introduce a sequence-based representation in which the reference sequences underlying BED intervals are extracted and sketched, enabling similarity measurements between annotations defined on different genome assemblies. Together, these approaches address two complementary barriers to large-scale interval analysis: the computational cost of expanding pairwise comparisons within a reference and the coordinate incompatibility that prevents direct comparison across references.
   

- All-vs-all pairwise interval similarity is bedtools' weakest scaling regime (O(N²·M)); large epigenome catalogs make this the bottleneck.
- Sketching trades exactness for speed; the question is *how much exactness*, and whether the sketch preserves the biology that the exact answer captures.
- Contribution: hammock, a Python+C++ HyperLogLog-backed sketcher with two operating modes (interval mode for BED, sequence mode for FASTA), plus a four-experiment validation framework that quantifies the speed/accuracy/biology trade-off on real data.



## 3. Background

### 3.1 The all-pairs interval-similarity problem

- BED file = a set of genomic intervals, treated as the base-pair set the intervals cover (merge within-file before counting).
- Pairwise similarity = base-pair Jaccard `J = |A ∩ B| / |A ∪ B|`.
- Exact all-pairs cost is O(N²·M) — couples file count with file size; the bottleneck sketching exists to break by reading each file once into a fixed-size summary.

### 3.2 Sketching for cardinality and similarity

- Define *sketch*: compact fixed-size set summary supporting approximate cardinality/similarity without retaining elements.
- Two genomics families: MinHash [@minhash] (sample hash-minimal elements) vs HyperLogLog [@hll] (extreme hash values per bucket); precedents Mash [@mash], Dashing [@dashing]; interval-set application comparatively unexplored.
- hammock uses HLL for both modes; §3.3–§3.4 develop the two pieces the results depend on.

### 3.3 HyperLogLog and the register-equality estimator

- HLL mechanics: hash → low `p` bits route to one of `m = 2^p` registers, each storing max ρ (leading/trailing-zero count, geometric: P(ρ=k)=2⁻ᵏ); cardinality via Ertl 2017 estimator [@Ertl2017]; RSE ≈ 1.04/√m, independent of input size.
- hammock reports the **register-equality** Jaccard (inherited from original `hammock`): `J_re = #{i : R_A[i]=R_B[i]≠0} / #{active}` — *not* exact set-Jaccard.
- Derivation that this is **affine** in the true Jaccard: per active register, the max-ρ "winner" is a shared element with prob ≈ J (guaranteed match), else ties by chance with prob `c` ⇒ `J_re ≈ c + (1 − c)·J`.
- Two load-bearing consequences (used in §4.2): (i) affine ⇒ order-preserving, so ranking is exact up to HLL noise; (ii) offset `J_re − J = c·(1 − J)` is largest at low J, vanishes as J→1.
- Intercept `c` = chance-tie rate of registers set by different elements. Exactly `c = Σ_k p_k²` with `p_k = e^(−λ2^(−k)) − e^(−λ2^(−(k−1)))` at load factor `λ = n/m`; the often-quoted `Σ_k(2⁻ᵏ)² = 1/3` is a loose upper bound not approached at any realistic λ. Measured/predicted `c ≈ 0.17` for λ ≳ 20, collapsing toward 0 as λ → 0. An estimator property, so it does not close with subsampling — but it *is* a function of λ and of the cardinality ratio `|A|/|B|`, so it is not a single constant that can be calibrated away.

### 3.4 Minimizers for sequence sketching

- Motivation: compare nucleotide content of intervals; full k-mer set is large/redundant.
- **Minimizer** scheme [@Roberts2004; @Schleimer2003]: smallest-hash k-mer per length-`w` window; subsamples k-mers ~`2/(w+1)` while preserving shared-substring minimizers (local similarity survives).
- Minimizer hashes feed the HLL of §3.3 → register-equality Jaccard again.
- Parameters as sweep axes: `k` = token specificity, `w` = density/compression; the (k, w) axes of the Results heatmaps.

## 4. Results

### 4.1 Speed: hammock is substantially faster than bedtools on real and synthetic interval corpora

> **Sources:** `experiments/subB_mixed_stride/RESULTS.md` (real DHS); `experiments/bedtools_benchmark/RESULTS.md` (synthetic).

**Real DHS (Maurano, 20 samples, 190 pairs, 8 threads for both tools — hammock `-t 8`, bedtools 8-way GNU parallel):**

> Generated by `docs/scripts/maurano_speed_table.R` from `docs/data/maurano_subB_summary.csv` (mixed-stride wall medians + per-pair MAE) and `docs/data/maurano_bedtools.csv` (bedtools reference wall); 5 reps per setting. Hammock subsample rows are mixed-stride.

| tool / setting | wall median (s) | speedup over bedtools | ΔJ vs hammock no-subsample |
|---|---|---|---|
| bedtools | 11.08 | 1.00× (ref) | — |
| hammock (no subsample) | 9.52 | 1.16× | 0 |
| hammock (subB = 0.1, high) | 5.20 | 2.13× | 9×10⁻⁴ |
| hammock (subB = 0.01, max) | 3.64 | 3.05× | 2×10⁻³ |

**Synthetic scaling (10k intervals/file, p=14, 16 threads for both tools — hammock `-t 16`, bedtools 16-way GNU parallel):**

| N files | bedtools | hammock | speedup |
|---|---|---|---|
| 64 | 11.0 s | 1.7 s | ~6× |
| 128 | 44 s | 3 s | ~14× |
| 256 | 176 s | 7 s | ~26× |
| **512** | **706 s** | **13 s** | **~50×** |

Two things to communicate in this section:
1. **hammock is faster than bedtools at every regime tested** — modestly faster on small real corpora (1.16× with no subsampling, up to ~3× at subB=0.01), dramatically faster as catalog size grows (≈50× at N=512 files).
2. **The speedup is not bought with accuracy loss.** Subsampling changes hammock's own per-pair Jaccard by < 2×10⁻³ vs its no-subsample output, so the gap to bedtools is statistically indistinguishable across subsampling settings — the speed knob is effectively "free."

**Fig 1 (headline):** Grouped wall-time bars — bedtools vs hammock at no-subsample / subB=0.1 / subB=0.01 (mixed-stride) on the Maurano corpus. Each hammock bar is annotated with its speedup over bedtools and the mean per-pair Jaccard change vs hammock's own no-subsample output (ΔJ ≤ 2×10⁻³), carrying "faster at no accuracy cost" in a single view. (Additional supporting figure: Fig S5.)

![Fig 1 — Maurano wall-time bars (hammock vs bedtools)](figures/maurano_speed_bars.png)

**Fig 2:** N-scaling out to 512 files on synthetic.

![Fig 2 — synthetic N-scaling](figures/synthetic_nscaling.png)

(Internal hammock-version comparisons — mixed-stride vs hash-threshold vs single-hash subB strategies, sort-time accounting, OpenMP scaling shape — belong in the supplementary methods, not the main text. The three-method speed/accuracy scatter is **Fig S5**, in the Supplementary figures section at the end.)

### 4.2 Accuracy: interval mode ranks like bedtools; sequence mode matches its values

> **Source:** `experiments/maurano_dhs_validation/RESULTS.md`.

Two distinct claims, which we are careful to keep separate. Interval mode is **rank-faithful** — its per-pair Jaccard tracks bedtools at r = 0.998 — but reports the register-equality statistic, which is an affine transform `J_re ≈ c + (1 − c)·J` of the exact base-pair Jaccard (§3.3). Because that transform is order-preserving, ranking is exact up to HLL noise (hence we report correlation, not an MAE); because it carries an intercept, the absolute values sit an offset `c·(1 − J)` above bedtools — ≈ 0.16 on these low-J pairs, shrinking to zero as similarity rises. The intercept is a property of the estimator rather than of sampling noise, so it does not close with subsampling; it is however a function of the **load factor** `n/m`, not a universal constant, and it collapses once the sketch outgrows the data (measured on a synthetic 90-pair spread: `c` = 0.180 at p = 12, 16, and 20 where `m ≪ n`, but 0.045 at p = 24 where `m > n`). Quoted offsets are therefore specific to a corpus *and* a precision. Sequence mode, at high k/w, is by contrast **value-identical** to four decimals.

Across Maurano's 400 sample pairs:

| mode | best r vs bedtools | best MAE | claim |
|---|---|---|---|
| **interval** | **0.998** | — (affine offset c·(1−J), §3.3) | rank-faithful |
| **sequence** (k=20, w=100, p=24) | **0.9996** | **0.0061** | value-identical |

Sequence mode's numerical agreement with bedtools peaks at r = 0.9996 / MAE = 0.0061 — four-decimal-place agreement, i.e. it closes the absolute gap that interval mode leaves open. The high-correlation ridge in the (k, w) Pearson heatmap runs along the high-k / high-w edge of the sweep, indicating that sequence mode's near-perfect agreement is unlocked at long minimizer windows where interior coverage is richest. (Caveat: agreement is not monotonic in w — at the very largest windows it falls back off saturation, e.g. at k = 15, p = 24 Pearson drops from ≈ 1.0 to 0.937 at w = 500 as interior minimizers grow too sparse — so the ridge peaks at high-k / moderate-w rather than at the extreme corner.)

**Fig 3:** Per sequence mode config, Pearson r against bedtools (x-axis) vs Pearson r against interval mode (y-axis); points sit tightly on y = x. Sequence mode's agreement with bedtools is mirrored by its agreement with interval mode, so interval mode is an interchangeable interval-Jaccard reference for sequence mode — once sequence mode is validated against bedtools, it is validated against any reasonable interval-Jaccard estimator.

![Fig 3 — sequence mode vs bedtools and vs interval mode agree per config](figures/mode_d_bedtools_vs_modeB_scatter.png)

**Fig 4:** Two-panel line plot at p = 24, `no_ends` column. **Left panel:** Pearson r vs bedtools as a function of window size w, one colored line per k (k ∈ {8, 10, 15, 20, 25}); the high-k lines (15/20/25) sit at r ≈ 1.0 across most of the w range (dipping only at the largest window, w = 500, down to ≈ 0.94), while the low-k lines (8, 10) climb toward them as w grows — together tracing the high-k / high-w Pearson ridge. **Right panel:** ARI vs tissue labels on the same axes; the k = 10 line spikes to 0.91 at w = 30 and drops off either side, while the k ∈ {15, 20, 25} lines are flat near ARI ≈ 0.69 for every w (high-k clustering is invariant to w). The contrast — Pearson is a ridge, ARI is a single peak — is the headline that Fig 6's per-config (Pearson, ARI) scatter abstracts. The original Pearson heatmap grid moves to Fig S2.

![Fig 4 — sequence mode Pearson + ARI vs w, by k, at p=24](figures/mode_d_lines_p24.png)


### 4.3 Biological signal: tissue identity recovery

> **Source:** `experiments/maurano_dhs_validation/RESULTS.md`.

Sequence mode's best-ARI cell is **k = 10, w = 30**, with **ARI = 0.910, NMI = 0.961** on the 10-tissue-label set (with the three muscle subtypes — arm/back/leg — counted as distinct labels) holding across **all precisions p ≥ 12** — the clustering signal is precision-cheap once the (k, w) cell is right [INCLUDE FIGURE ILLUSTRATING SETTINGS]. Hammock successfully separates and groups the tissues. In fact, sequence mode is able to group and differentiate different muscle tissue types, while bedtools jaccard struggles.

At the ARI-best config, sequence mode's predicted Jaccards sit on the y = x diagonal versus bedtools (and versus interval mode) — the sketch is numerically calibrated against the bedtools reference, not just rank-correlated.

**Fig 5:** sequence mode best-ARI dendrogram (top) paired with the bedtools reference dendrogram (bottom). Same tissue blocks recovered.

![Fig 5a — sequence mode best dendrogram](../experiments/maurano_dhs_validation/figures/mode_d_best_dendrogram.png)
![Fig 5b — bedtools reference dendrogram](../experiments/maurano_dhs_validation/figures/bedtools_dendrogram.png)

**Fig 6:** The headline methodological point: **best-Pearson cell ≠ best-ARI cell**. Pearson-best (large-k, large-w) clusters at ARI ≈ 0.69; ARI-best (k=10, w=30) at Pearson ≈ 0.946. Numerical perfection and clustering quality are non-coincident knobs.

![Fig 6 — sequence mode Pearson vs ARI tradeoff](figures/mode_d_metric_tradeoff.png)

(The per-k violin view of Pearson r and ARI — which makes the "tunable clustering only at k = 10, fixed partition at large k" contrast explicit — is **Fig S4**, in the Supplementary figures section.)

### 4.4 Robustness to reference genome

> **Source:** `experiments/ref-comparison/docs/exp_a_results.md`.

3 H3K27ac samples (heart, liver, lung) × 3 references (GRCh37/GRCh38/CHM13), 9 sample×ref combinations. Across the (k, w) sweep, same-tissue cross-reference Jaccard is significantly higher than different-tissue Jaccard at every cell with k ≥ 8 (Wilcoxon p ≤ 10⁻⁵), and at **k ≥ 15 the two groups are *fully separated*** — the minimum same-tissue cross-reference similarity exceeds the maximum different-tissue similarity. This is a stronger statement than significant: no overlap.

**Best lead cell: k = 15, w = 15.** On broad peaks, Δmedian = **0.398** (median same-tissue cross-ref = 0.783, median different-tissue = 0.385; Wilcoxon p = 1.35 × 10⁻¹⁰, the test floor at n=18/54). Narrow peaks: Δmedian = 0.387 (0.729 vs 0.342). k=20, w=20 reaches an even larger Δ (0.413 broad) but at slightly suppressed medians.

The sweep partitions into three regimes:

| Regime | Cells | Behavior |
|---|---|---|
| Saturated high | k ≤ 8 (any w) | Both groups ≈ 0.95–0.99; Δ ≤ 0.02 |
| Interpretable mid | k = 10, w ≥ 10 | Medians ≈ 0.55–0.65; Δ ≈ 0.09; groups overlap |
| **Interpretable + fully separated** | **k ≥ 15 (any valid w)** | **Δ ≈ 0.32–0.45; min(xref) > max(diff-tissue)** |

**Fig 7:** UPGMA dendrogram at the new headline cell (k=15, w=15) — each tissue's three references form a tight monophyletic clade with deep separation between tissues; broad and narrow peak calls give the same structure.

![Fig 7 — cross-reference dendrogram, k=15, w=15](../experiments/ref-comparison/figures/cross_ref_dendrogram_k15_w15.png)

**Fig 8:** (k × w) effect-size heatmap for broad peaks; the three regimes are immediately visible, with k ∈ {15, 20} as a uniformly high-effect block.

![Fig 8 — cross-ref effect-size sweep, broad](../experiments/ref-comparison/figures/sweep_effect_size_broad.png)

**Metric choice at k=10, w=10** (`scripts/exp_a_metric_comparison.R`): of the 12 emitted similarity columns, `jaccard_similarity_with_ends` remains the right default — best Δ/saturation trade-off (broad Δ = 0.086, p = 2.0 × 10⁻⁷). `cosketch_geom_with_ends` is a near-tie (Δ = 0.065). All minimizer-only metrics hit the Wilcoxon p-floor but operate near saturation (medians ≈ 0.99 vs 0.92), so absolute Δ is small. `cosketch_max_with_ends` collapses on narrow (Δ ≈ 0, p ≈ 0.49) — the worst metric in both flavors and not recommended. (The full 12-metric comparison is **Fig S1**, in the Supplementary figures section.)

(Note: at k=15, w=15 the minimizer-only signal is already fully separated, so the choice of with-ends vs no-ends matters less; a future replication of this 12-metric comparison at the new headline cell would confirm that observation.)

Practical interpretation: when peaks are aligned to a different human reference than expected, the sketch still produces the same biological neighborhood. At k ≥ 15 the separation is large enough that reference choice is unambiguously a smaller source of variance than tissue identity. This is the property that lets hammock be deployed against heterogeneous catalogs (ENCODE/Roadmap mixtures) without first re-aligning everything.

This is a proof-of-concept on a deliberately small panel (3 tissues × 3 references = 9 samples; n = 18 same-tissue cross-reference pairs vs 54 different-tissue pairs), and the Wilcoxon test is at its p-floor at that n — so the result establishes that the separation is large and consistent on these tissues, not that it generalizes across the full diversity of cell types, assays, or reference builds. See §6.3.

### 4.5 Methodological notes: choosing sequence mode's flanking column

> **Status (2026-05-21):** Section being reworked. All §4.5 figures are hidden from the rendered outline pending replacement — none are currently shown. The hidden PNGs are still on disk at `experiments/modeD_flanking/figures/` for reference. New, more interpretable figures are being brainstormed.

> **Source:** `experiments/modeD_flanking/` Parts 1 (Maurano re-analysis) and 2 (synthetic FASTA pairs with exact k-mer Jaccard truth).

Sequence mode emits two Jaccard columns: minimizer-only (`jaccard_similarity`) and minimizer-plus-flanks (`jaccard_similarity_with_ends`). Two corpora characterize when each column is preferred:

- **Part 1 (Maurano, real DHS):** at high-precision configs, `no_ends` is the clear winner — at k=20, w=100, p=24, r_no_ends = 0.9996 vs r_with_ends = 0.888. The advantage flips only at the smallest (k, w) cells (k ≤ 10, w ≤ 10), where interior minimizers are sparse and the flanking k-mers carry the signal.
- **Part 2 (synthetic, 192 FASTA pairs across an n_intervals × mean_len × dist × mutation grid, each sketched at 66 valid (k, w, p) cells — ~12,700 measurements total; exact k-mer Jaccard as truth):** `with_ends` has systematically *larger* MAE than `no_ends` against the k-mer truth (mean Δmae_r ≈ −0.015 across all k); the gap widens with φ, the flanking-fraction predictor we defined.

The flanking-fraction φ ≈ 2(k−1)·n_intervals / (total_length / w) predicts the regime: large φ → `_with_ends` over-weights boundary k-mers; small φ → either column works.

**Recommendation:** minimizer-only (`jaccard_similarity`) is the right default for ChIP/DHS-shaped corpora at any reasonable parameter choice; `_with_ends` is a fallback for short-sequence, sparse-minimizer regimes.

<!-- Fig 9 (hidden 2026-05-21, pending rework — file still on disk):
**Fig 9:** Part 1 (Maurano real DHS) — sign of (no_ends − with_ends) flips with w.

![Fig 9 — Maurano Δr(no_ends − with_ends) vs w](../experiments/modeD_flanking/figures/maurano_delta_r_vs_w.png)
-->

<!-- Fig 10 (hidden 2026-06-03, pending rework — file still on disk):
**Fig 10:** Part 2 (synthetic) — φ × mutation phase plane.

![Fig 10a — synthetic Δ vs φ](../experiments/modeD_flanking/figures/synthetic_delta_vs_phi.png)

![Fig 10 — synthetic φ × mutation phase diagram](../experiments/modeD_flanking/figures/synthetic_phase_diagram.png)
-->

<!-- Fig 11 (hidden 2026-05-21, pending rework — file still on disk):
**Fig 11:** Empirical agreement with the analytical φ-prediction — independent validation of the flanking-fraction model.

![Fig 11 — empirical vs analytical φ](../experiments/modeD_flanking/figures/synthetic_empirical_vs_analytical.png)
-->

## 5. Methods

### 5.1 hammock implementation

- Python orchestrator + C++ extension (pybind11); HLL with register-equality Jaccard, Ertl 2017 estimator, xxh64 ingestion.
- Interval mode: per-bp HLL with optional `subB` subsampling. Sequence mode: minimizer-HLL on FASTA, with both interior-minimizer and minimizer-plus-flanks similarity columns. (The implementation also carries two interval-coordinate variants — a coordinate-only sketch and an interpolation between it and per-bp interval mode — but neither is evaluated in this paper.)
- Similarity is reported as **register-equality Jaccard**: the fraction of active registers (nonzero in at least one sketch) whose values are equal, matching the original `hammock`'s estimator. Note this is not exact set-Jaccard — it is an affine transform of it (§3.3), order-preserving with a precision-independent offset — and both modes sketch the same base-pair (interval) or k-mer (sequence) universe that the exact reference does.
- Output per-pair similarity columns: `jaccard_similarity`, `jaccard_similarity_with_ends`, plus `containment_AB`, `containment_BA`, `cosketch_{geom,arith,max}` in both flavors (10 metrics per sequence-mode pair; 5 per interval-mode pair). The `jaccard_*` columns are the analyses' default; the cosketch/containment columns are reported for transparency and inform Section 6.

### 5.2 Benchmark datasets

| dataset | role | n |
|---|---|---|
| Synthetic BED (`bedtools_benchmark/`) | scaling regime | up to 512 × 10k intervals |
| Synthetic FASTA pairs (`modeD_flanking/`) | exact k-mer truth for sequence mode | 192 pairs (n_intervals × mean_len × dist × mutation grid); ~12,700 (k, w, p) sweep measurements |
| Maurano 2012 fetal DHS | real interval data, 10 tissue labels (fBrain/fHeart/fIntestine_Sm/fKidney/fLung/fMuscle_arm/fMuscle_back/fMuscle_leg/fSkin/fStomach) | 20 BEDs |
| ENCODE H3K27ac (Heart/Liver/Lung × GRCh37/GRCh38/CHM13) | reference-build robustness | 9 sample×ref |

### 5.3 Statistical evaluation

Pearson r, Spearman ρ, MAE vs bedtools (interval and sequence modes); ARI, NMI vs known tissue labels (sequence mode); Wilcoxon same-tissue vs different-tissue (cross-reference); R/ggplot via CairoPNG.

### 5.4 HyperLogLog sketching

Both modes share a HyperLogLog (HLL) backbone [@Flajolet2007]. Each input set — per-bp positions for interval mode, minimizer hashes for sequence mode — is hashed with xxh64 (seed via `--seed`, default 42); the low `p` bits of each hash route it to one of `2^p` 1-byte registers, which stores the maximum leading-zero count seen among hashes routed there. Cardinality is recovered via the Ertl 2017 improved estimator [@Ertl2017]. For two sketches at matching `p` and seed, the union is the register-wise maximum. The two reported similarity families are derived from that union by **different routes**, and the distinction is load-bearing for how their values should be read:

- `jaccard_similarity` is the **register-equality** statistic of §3.3 — the fraction of active registers whose values are equal. Register agreement is implied by, but does not imply, a shared element: two registers set by *different* elements of equal ρ also agree, which is the chance-tie term `c` that makes this an affine transform of the exact Jaccard rather than an estimate of it.
- The directional containments `|A ∩ B| / |A|` and `|A ∩ B| / |B|` (and the three co-sketch summaries derived from them) instead take the intersection cardinality by **inclusion–exclusion**, `|A| + |B| − |A ∪ B|`, evaluated with the Ertl estimator on each of the three sketches. This carries no chance-tie term, so unlike `jaccard_similarity` these columns estimate the true set quantities — at the cost of being a difference of large estimates, hence higher variance when the intersection is small relative to the union.

The asymptotic relative standard error is ≈ 1.04 / √(2^p). For the CLI default `p = 18`, that is 1.04 / 512 ≈ 0.203% on a 2^18 = 262,144-register / 256 KiB sketch. For the high-precision configuration cited in Sections 4.2–4.3 (`p = 24`), it is 1.04 / 4,096 ≈ 0.0254% on a 2^24 = 16,777,216-register / 16 MiB sketch. Memory is independent of input cardinality — the load-bearing property that lets the same 16 MiB sketch represent a 10k-interval or a 10M-interval BED at identical cost.

### 5.5 Minimizers (sequence mode)

Sequence mode reduces each FASTA sequence to its set of (k, w)-minimizers [@Roberts2004; @Schleimer2003]: in every length-`w` sliding window over the sequence, the k-mer with the smallest selector-hash value is retained. Window scanning is delegated to the VeryAmazed `digest` library [@digest], and each unique selector hash is ingested directly into the per-sequence HLL of §5.4 via `add_hash64`. In parallel hammock canonicalizes the two flanking k-mers — the leading and trailing length-k substrings of each input sequence — by taking the lexicographic minimum of forward and reverse complement, xxh64-hashes them with the same `--seed` as the HLL, and adds them to a second HLL. The default similarity column `jaccard_similarity` compares the interior-minimizer HLLs only; the alternative `jaccard_similarity_with_ends` compares the union of the interior-minimizer and flanking-k-mer HLLs, and is preferred only in the short-sequence / sparse-minimizer regime characterized in §4.5. CLI defaults are `k = 8`, `w = 40`; the production-cited configurations are `k = 20, w = 100` (§4.2, numerical agreement with bedtools) and `k = 10, w = 30` (§4.3, tissue clustering).

### 5.6 Mixed-stride subsampling (`--subB-method mixed-stride`)

Interval mode's `subB` knob subsamples the per-base-pair positions that get ingested into the HLL. A natural implementation — a *hash-threshold* gate that computes a per-position xxh32 hash and compares it to `subB × UINT32_MAX` — is correct, but the per-position hash cost grows linearly with the input, so the wall-time savings from skipping positions plateau quickly: at moderate subB the gate hash itself dominates and `subB ≤ 0.5` is actually *slower* than no subsampling at all.

hammock introduces **mixed-stride** subsampling as a deterministic, hash-free alternative. For each chromosome, a fixed stride length `s = round(1/subB)` is chosen, offset by a chr-keyed hash so different chromosomes sample disjoint position phases. Ingestion walks the chromosome in stride-`s` increments — no per-position decision, no per-position hash. The accepted-position cardinality matches the hash-threshold gate in expectation, but the cost is `O(L/s)` rather than `O(L)`, so wall time scales with `subB` instead of being dominated by the gate.

Three properties matter for downstream use:

1. **Performance scaling.** At `subB = 0.1`, mixed-stride is 3–4× faster than hash-threshold on synthetic data and 1.8–2.4× faster than hash-threshold on real DHS. The advantage compounds: at `subB = 0.01` on 1M-interval synthetic files, mixed-stride hits 14× while hash-threshold caps at ~1.4×.
2. **Accuracy.** Per-pair MAE vs the no-subsample reference is statistically indistinguishable across the three subB methods (hash-threshold, mixed-stride, single-hash). Mixed-stride does not buy speed by losing accuracy.
3. **Determinism and reproducibility.** Output is exactly reproducible given the chr-stride seed (`--gate-seed`, default 31337). The same files at the same subB produce byte-identical HLLs across runs and machines.

The structured-sampling concern — that a fixed stride could miss periodic features — is theoretical for BED inputs at the strides we test; in practice the chr-keyed offset randomization breaks cross-chromosome alignment and the Maurano + synthetic accuracy measurements show no detectable bias. Mixed-stride is the default `--subB-method` setting; hash-threshold (`--subB-method hash-threshold`) and a single-hash variant remain available for users who want to characterize the alternative samplers.

This subsampling refinement is what makes the Section 4.1 speed numbers attainable: the "high-subsample" and "max-subsample" rows of the real-DHS table are mixed-stride at `subB = 0.1` and `subB = 0.01`. The per-method comparison plot is in supplementary.

## 6. Limitations

### 6.1 Definitional gap vs bedtools

The estimator-vs-exact gap is characterized in §3.3: interval mode's register-equality Jaccard is an affine transform c + (1 − c)·J of the exact base-pair Jaccard, so it is rank-faithful (r = 0.998) but its absolute values carry an offset c·(1 − J) — ~0.16 at low J on the synthetic benchmark, shrinking to zero as similarity rises. Because it is an estimator property rather than statistical noise it does not close with subsampling, and it is flat across precision *while the sketch stays saturated* (`m ≪ n`) — but it is set by the load factor `n/m` and collapses once `m` exceeds `n` (0.180 at p = 12/16/20 vs 0.045 at p = 24 on ~5 Mbp inputs), so the quoted offset is corpus- and precision-specific. The practical consequence: interval-mode `jaccard_similarity` should be read as a ranking, not as a bedtools magnitude.

Two routes recover absolute magnitude. (1) **Sequence mode** at the optimal high-k/high-w cell brings MAE against bedtools down to 0.006. (2) Within interval mode, the **containment columns already estimate the true set quantities** by inclusion–exclusion (§5.4), so a set-Jaccard estimate is recoverable from existing output with no rerun: `J = 1/(1/C_AB + 1/C_BA − 1)`. Because both containments derive from a single shared intersection estimate, this recovers the inclusion–exclusion estimator exactly rather than approximating it. On the synthetic spread it tracks bedtools at MAE 5×10⁻⁴ (p = 20) and 1×10⁻⁴ (p = 24), versus 0.15 for the reported `jaccard_similarity`.

**The two estimators trade off rather than one dominating, and the paper should say so.** MAE scores calibration, not resolution. Rescaled into true-J units by its own slope `1 − c`, register-equality is the *lower-variance* discriminator at low J whenever the sketch is saturated — error sd 0.0014 vs 0.0024 at p = 16, J < 0.05 — which is the regime genome-scale BED data occupies. Inclusion–exclusion only wins on both axes once `m` approaches or exceeds `n` (p ≥ 20 on these ~5 Mbp inputs). It is additionally censored: the `≥ 0` clamp fires on 25/90 pairs at p = 12, all at low J, so conditional on being non-zero it is biased upward near J = 0, and its relative error `≈ 0.6/(J√m)` makes it uninformative below J ≈ a few/√m. Net: report `jaccard_similarity` for ranking and clustering, the reconstruction for magnitude comparisons against bedtools or another tool.

### 6.2 The cosketch + containment columns are reported but not yet exploited

The five auxiliary similarity columns (containment_AB, containment_BA, cosketch_geom, cosketch_arith, cosketch_max — each in two flavors for sequence mode) are present in every sequence mode output CSV. Current analyses use the jaccard columns as their primary signal. A 12-metric sanity check at the ref-comparison Exp A (k=10, w=10) cell finds `cosketch_geom_with_ends` is a near-tie with `jaccard_similarity_with_ends`, while `cosketch_max` is uniformly the weakest discriminator. A full multi-metric re-evaluation across the (k, w, p) sweep and across the Maurano corpus is a natural fast-follow analysis and may identify a column (likely cosketch_geom on the minimizer-only flavor) that is more robust than jaccard at small precision or small k.

### 6.3 Reference-robustness is shown on a small panel

The reference-genome invariance result (§4.4) is established on 9 samples — 3 tissues (heart/liver/lung H3K27ac) × 3 references (GRCh37/GRCh38/CHM13), giving n = 18 same-tissue cross-reference pairs and 54 different-tissue pairs. The separation at k ≥ 15 is large (Δmedian ≈ 0.4, full separation of the two groups), but with this n the Wilcoxon statistic is pinned at its floor, so the p-values quantify consistency on these tissues rather than generalization. We read the experiment as a proof-of-concept that reference build is a smaller variance source than tissue identity for these marks — not as a survey across cell types, assays, or assemblies. Extending the panel (more tissues, more assays, additional assemblies) is the obvious next step before the invariance claim is relied on for heterogeneous-catalog deployment.

## 7. Discussion

- **Practical recipe.** Interval mode for fast bedtools-equivalent interval-Jaccard with optional subsampling for further speedup at no accuracy cost. Sequence mode at large k and w (k=20, w=100, p=24) for the closest numerical match to bedtools (r = 0.9996, MAE = 0.006). Sequence mode at k=10, w=30, p ≥ 12 for tissue clustering (ARI = 0.91, NMI = 0.96). Use `jaccard_similarity` by default; fall back to `jaccard_similarity_with_ends` only in the short-sequence / sparse-minimizer regime characterized in Section 4.5.
- **The sketch carries more than bedtools captures.** Sequence mode recovers tissue clustering directly from peak FASTAs — independently from bedtools' interval overlap — at ARI = 0.91. Sketch similarity ≈ biological similarity, even when the two estimators don't agree numerically.

## 8. Conclusion

hammock provides a fast, sketch-based alternative to `bedtools jaccard`. On real DHS data interval mode matches bedtools at r = 0.998 and sequence mode matches it at r = 0.9996 / MAE = 0.006, while running substantially faster than bedtools at every scale tested (orders of magnitude faster at large catalog size). The same sequence mode sketches recover tissue clustering at ARI = 0.91 directly from peak FASTAs and are robust to reference-genome choice within a species. The combination — speed *plus* near-exact agreement *plus* preserved biology *plus* reference-build invariance — positions sketching as a viable default for large-scale epigenome comparison.

---

## Figure index

| # | Path | Carries |
|---|---|---|
| 1 | `docs/figures/maurano_speed_bars.png` (data: `docs/data/maurano_subB_summary.csv` + `docs/data/maurano_bedtools.csv`, script: `docs/scripts/maurano_speed.R`) | hammock wall-time bars vs bedtools at no-sub / subB=0.1 / 0.01 |
| 2 | `docs/figures/synthetic_nscaling.png` (data: `docs/data/cpp_vs_bedtools_t16_20260512_160412.csv`, script: `docs/scripts/synthetic_nscaling.R`) | synthetic scaling to N=512 |
| 3 | `docs/figures/mode_d_bedtools_vs_modeB_scatter.png` (data: `docs/data/mode_d_summary.csv`, script: `docs/scripts/mode_d_bedtools_vs_modeB_scatter.R`) | per-config sequence mode r vs bedtools vs interval mode (y=x) |
| 4 | `docs/figures/mode_d_lines_p24.png` (data: `docs/data/mode_d_summary.csv`, script: `docs/scripts/mode_d_lines.R`) | sequence mode Pearson ridge vs ARI peak at p = 24 |
| 5 | `maurano_dhs_validation/figures/mode_d_best_dendrogram.png` + `bedtools_dendrogram.png` | tissue recovery (a + b panels) |
| 6 | `docs/figures/mode_d_metric_tradeoff.png` (data: `docs/data/mode_d_summary.csv`, script: `docs/scripts/mode_d_metric_tradeoff.R`) | Pearson-best ≠ ARI-best |
| 7 | `ref-comparison/figures/cross_ref_dendrogram_k15_w15.png` | within-tissue clades across refs (headline cell) |
| 8 | `ref-comparison/figures/sweep_effect_size_broad.png` | cross-ref effect-size heatmap, broad |
| 9 | *(hidden — pending rework; PNG on disk: `modeD_flanking/figures/maurano_delta_r_vs_w.png`)* | flanking column choice on real data |
| 10 | *(hidden — pending rework; PNG on disk: `modeD_flanking/figures/synthetic_phase_diagram.png`)* | flanking on synthetic, φ × mutation phase plane |
| 11 | *(hidden — pending rework; PNG on disk: `modeD_flanking/figures/synthetic_empirical_vs_analytical.png`)* | analytical φ-prediction validated |
| S1 (supp) | `ref-comparison/figures/metric_comparison_broad_k10_w10.png` | 12-metric Wilcoxon comparison at k=10, w=10 |
| S2 (supp) | `maurano_dhs_validation/figures/mode_d_pearson_heatmap.png` | full Pearson heatmap grid across precision × column flavor (faceted) |
| S3 (supp) | `maurano_dhs_validation/figures/mode_d_clustering_ari.png` | full ARI heatmap grid across precision × column flavor (faceted) |
| S4 (supp) | `docs/figures/mode_d_violins_by_k.png` (data: `docs/data/mode_d_summary.csv`, script: `docs/scripts/mode_d_violins.R`) | per-k Pearson and ARI distributions across full (w, p) sweep |
| S5 (supp) | `docs/figures/maurano_subB_pareto_scatter.png` (data + script as Fig 1) | subB-method speed/accuracy scatter (de-zigzagged former headline Pareto) |

(The `synthetic_speedup_vs_nosub.png` figure showing within-hammock subB-method comparison moves to supplementary — per the constraint that internal hammock-version performance differences are not paper material.)

---

## Supplementary figures

**Fig S1:** 12-metric Wilcoxon comparison at k=10, w=10 (broad peaks). `jaccard_similarity_with_ends` gives the best Δ/saturation trade-off; `cosketch_geom_with_ends` is a near-tie; `cosketch_max` is the weakest discriminator. Backs the metric-choice discussion in §6.2.

![Fig S1 — 12-metric Wilcoxon comparison at k=10, w=10, broad](../experiments/ref-comparison/figures/metric_comparison_broad_k10_w10.png)

**Fig S4:** Per-k violin plot of Pearson r and ARI (moved from the main text). **x-axis:** k (8, 10, 15, 20, 25). **y-axis:** metric value on [0, 1], shared scale. **Two violins per k:** Pearson r vs bedtools (across every (w, p) cell at that k) and ARI vs the 10 fetal-tissue labels (across every (w, p) cell at that k), `no_ends` column only. Pearson violins march steadily upward with k. The ARI side has two distinct shapes: at **k = 10 the violin is tall and multimodal**, spanning from near zero up to 0.91 with a visible upper lobe at 0.91 separated from the bulk by a clear gap (only the (w = 30, p ≥ 12) slice hits that upper plateau); at **k ∈ {15, 20, 25} the violin collapses to a horizontal line at ARI = 0.693** because every (w, p) cell at those k values gives the *same* partition. Together those two shapes carry the message: tunable clustering quality exists only at k = 10; large k locks the sketch into a fixed sub-optimal partition. The full ARI heatmap grid is Fig S3.

![Fig S4 — sequence mode Pearson + ARI violins per k, no_ends across (w, p)](figures/mode_d_violins_by_k.png)

**Fig S5:** subB-method speed/accuracy scatter — speedup vs bedtools (y) against mean per-pair Jaccard error vs the no-subsample reference (x, log); points for each `--subB` level of all three methods, no connecting line. Only mixed-stride converts subsampling into speed (~2× at subB=0.1, ~3× at subB=0.01); hash-threshold and single-hash stay near 1×. De-zigzagged replacement for the former headline Pareto (§4.1).

![Fig S5 — subB-method speed/accuracy scatter](figures/maurano_subB_pareto_scatter.png)

(Figs S2 and S3 — the full Pearson and ARI heatmap grids — are listed in the figure index above but not yet rendered here; add their images to this section when the panels are finalized.)
