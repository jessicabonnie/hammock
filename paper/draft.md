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
> from 27.61× to roughly 10.2×.** That number was itself provisional and has
> since moved again: a small-N noise bug in the same rerun was found and fixed,
> landing at **9.21× at N=512 on the register-equality (`jaccard_similarity`)
> column** — see Supplementary Figure S9. As of the 2026-08-10 figure
> simplification, main-text **Figure 3 itself reports `jaccard_similarity_ie`
> only, landing at 8.4× at N=512** — see the Figure 3 section below, which is
> current.
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

This document used for brainstorming and drafting. For current outline with all draft text see: https://www.overleaf.com/read/kggqsxtxrbdy#4685d0
## Thesis

Modern interval collections are limited by two boundaries: the computational cost of comparing every dataset and the coordinate dependence that prevents comparison across references. Hammock creates reusable sketches of either interval coordinates or interval-derived sequences, expanding comparison within and across those boundaries while preserving useful biological structure.

## Abstract

## I. Introduction

Increased availability of high-throughput sequencing has facilitated a corresponding increase in large-scale modern sequencing projects which produce vast numbers of genomic annotation datasets which then, themselves, become part of the genomic data ecosystem. Large-scale initiatives such as the ENCODE Project\cite{encode_2012}, the Roadmap Epigenomics Project\cite{roadmap_2015}, GTex[REF], and the 1000 Genomes Project\cite{10002015global} make datasets publicly available, providing unprecedented opportunities for integrative analysis. Among the most common and biologically informative outputs of these projects are genomic intervals: contiguous stretches of the genome that define transcription factor binding sites, chromatin accessibility peaks, histone modifications, structural variants, and splicing junctions. Such interval datasets have become a cornerstone of downstream analyses that seek to characterize genome function and variation across populations, tissues, cell types, and experimental conditions. In particular, pairwise comparison of interval datasets presents a measurement of biological similarity, such as identifying conserved regulatory elements between tissues. Tools such as \program{BEDTools} provide exact calculations of these measures\cite{quinlan2010bedtools}, but the computational cost of pairwise overlap calculations grows rapidly with the size of modern repositories. In practice, this creates a scalability bottleneck that hampers systematic comparisons across the interval datasets now available\cite{li2020design}.

The numbers of files in these interval databases continue to grow every year. ChIP-Atlas, for example, expanded from 37,720 accumulated experiments in 2015 to 464,655 in 2025, while repositories such as ENCODE contain thousands of additional ChIP-seq and chromatin-accessibility experiments.

Scalability is not the only obstacle to systematic interval comparison. BED coordinates are meaningful only with respect to the reference genome on which they were defined, and large public collections remain distributed across multiple genome assemblies. Standard overlap measures therefore cannot directly compare two interval files when one is defined on hg19 and the other on hg38, even when the files describe the same assay or biological feature. Figure 1 illustrates this limitation for histone-mark ChIP-seq peak files from Roadmap Epigenomics and BLUEPRINT. For histone marks represented in both collections, within-Roadmap and within-BLUEPRINT comparisons can be performed in their respective coordinate systems, whereas every Roadmap–BLUEPRINT pairing is blocked by the hg19–hg38 mismatch. Coordinate conversion can sometimes recover such comparisons, but it introduces additional preprocessing, may not map all intervals unambiguously, and continues to define similarity in terms of a selected coordinate system. An alternative is to extract the reference sequence underlying each interval and compare the resulting sequence collections directly, allowing related genomic annotations to be evaluated across references without requiring their coordinates to coincide.

![Figure 1](figures/roadmap_blueprint_top5_marker_pairwise_comparisons.png)

**Figure 1. Reference-genome fragmentation prevents direct comparison across public histone-mark collections.** Numbers of possible pairwise comparisons among processed histone-mark ChIP-seq peak BED files are shown for the five most represented marks shared by Roadmap Epigenomics and BLUEPRINT. Within-resource comparisons can be performed directly because Roadmap files share hg19 coordinates and BLUEPRINT files share hg38 coordinates. In contrast, every Roadmap–BLUEPRINT pairing is blocked for direct coordinate-overlap analysis by the hg19–hg38 reference mismatch. File counts are deduplicated at the download-URL level; included outputs are BED-, narrowPeak-, broadPeak-, or gappedPeak-like peak files. Repository counting and inclusion criteria are described in Supplementary Methods S1.

Probabilistic data structures, or sketches, offer a powerful solution to these challenges. Sketching methods compress large datasets into compact representations that allow efficient estimation of set cardinality and similarity. Techniques such as MinHash\cite{minhash} and HyperLogLog\cite{hll} have seen wide application in computer science, and their introduction to genomics has already revolutionized k-mer–based comparisons of sequencing datasets. Tools like Dashing \cite{dashing} and Mash\cite{mash} have demonstrated the value of sketching for rapid, large-scale sequence and metagenome comparisons. Extending these methods to genomic interval data, however, remains relatively unexplored. Because intervals represent continuous spans rather than discrete tokens, adapting sketches for overlap estimation presents unique methodological challenges but also significant opportunities. Moreover, representing intervals through both their coordinates and their underlying reference sequences provides complementary notions of similarity: coordinate-based sketches support rapid comparison within a shared reference, whereas sequence-based sketches can support comparisons across references.

In this study, we present \program{hammock}, a command-line tool for scalable comparison of genomic interval datasets. Within a shared reference genome, \program{hammock} applies probabilistic sketches to BED intervals to approximate overlap and similarity across large collections of files, reducing the computational burden of systematic pairwise analysis. We benchmark these estimates against exact calculations from established interval-processing tools and evaluate the resulting trade-offs among speed, memory use, and accuracy. We additionally introduce a sequence-based representation in which the reference sequences underlying BED intervals are extracted and sketched, enabling similarity measurements between annotations defined on different genome assemblies. Together, these approaches address two complementary barriers to large-scale interval analysis: the computational cost of expanding pairwise comparisons within a reference and the coordinate incompatibility that prevents direct comparison across references.

## II. Results

### 2.1 Hammock represents interval sets in complementary coordinate and sequence spaces
\program{hammock} converts each input dataset once into a reusable file-level sketch. In interval mode the sketch represents genomic positions covered by the intervals, while in sequence mode, it represents minimizer-derived sequence features. As shown in Figure \ref{fig:hammock_workflow} both representations use compact HyperLogLog (HLL) sketches which are then available for pairwise comparisons in place of the original files. The subsequent analyses evaluate whether interval sketches preserve overlap similarity at scale and whether sequence sketches preserve biological relationships across reference genomes.

![Figure 2](figures/hammock_workflow.png)

**Figure 2. Hammock provides complementary coordinate- and sequence-based representations of genomic interval sets.** (A) Public interval collections are large, comparisons scale quadratically with the number of files, and BED coordinates are tied to specific reference genomes. (B) In interval mode, each BED file is summarized as a reusable sketch of covered genomic positions, enabling fast all-pairs similarity comparisons within a shared reference. (C) In sequence mode, interval sequences are extracted from each file's native reference FASTA, representative k-mers are selected using minimizers, and the resulting sequence sketches are compared across references without requiring direct coordinate overlap. The two modes therefore answer complementary questions: whether interval sets occupy similar genomic locations and whether they contain similar underlying sequence.

### 2.2 Coordinate Space 
#### 2.2.1 Interval sketching expands feasible all-pairs BED comparison

[NOTE: find a place for this] Interval comparisons are typically conducted using files in the BED format, a simple tab-delimited representation of genomic intervals that has been in wide use since the early 2000s\cite{kent_hgbrowser}.

The Jaccard statistic, from set theory, [ETC] is one commonly used measure of similarity, offered by a number of tools. In particular, BEDtools [REF], a ubiquitous tool for genomic arithmetic, offers a `bedtools jaccard` command that calculates the exact value of intersection over union [this might now be the exact approach... subtraction is maybe used] for two BED files. \program{hammock}, which is built to compare *lists* of interval files, takes advantage of two aspects of a point-wise sketching representation approach: batch reusability of a sketch once built and point-based reproducibly random subsampling.
As can be seen in Figure 3, the majority of \program{hammock}'s time is spent summarizing each interval file into a sketch, but the pairwise calculations using and reusing the resulting bit-vectors are trivially fast. BEDtools, by contrast, must reprocess each interval file for each new comparison. As one would expect, this means the \program{hammock}'s performance advantage scales with the number of BED files included in the analysis.

\program{hammock}'s \texttt{interval} mode treats each point within an interval as an object to be hashed and added to a single sketch representation of the entire set of intervals in the file. Given the usual size of genomics intervals, this creates an opportunity for further subsampling prior to sketching.  \program{hammock}'s novel `subB` implementation uses mixed-stride methodology to sample genomic positions using deterministic chromosome-specific strides and offsets, avoiding the per-position hashing cost of other reproducibly random methods such as hash-threshold subsampling. As seen in Figure 3, this substantially reduces sketch-construction time while preserving reproducibility and producing little change relative to the unsubsampled hammock similarities in the evaluated datasets. The reduction is sublinear in the sampling fraction and corpus-dependent: a 10-fold subsample buys 4.21\times{} on the synthetic corpus at p=14 but only 1.83\times{} on Maurano at p=18.

Performance improvements are not uniform.


This section should distinguish three sources of the interval-mode speedup:

1. **Sketch reuse:** each BED file is read once and converted to a fixed-size summary that can be reused across all pairwise comparisons.
2. **Mixed-stride subsampling:** hammock's novel `subB` implementation samples genomic positions using deterministic chromosome-specific strides and offsets, avoiding the per-position hashing cost of hash-threshold subsampling. This substantially reduces sketch-construction time while preserving reproducibility and producing little change relative to the unsubsampled similarities. The reduction is sublinear in the sampling fraction and not stable across corpora — a 10-fold subsample buys 4.21x on synthetic at p=14 against 1.83x on Maurano at p=18 — so it must not be described as proportional.
3. **Workflow parallelization:** The saturation is a result, not an obstacle. It's a genuine scalability argument for sketching: an exact per-pair tool must launch N² processes and re-read each file N times, and that workflow does not parallelize — hammock launches one process and reads each file once. That's a real, defensible advantage. It just has to be attributed to batch-mode absence rather than to sketching, and reported as bedtools' achieved efficiency rather than silently labeled "t=16."
- hammock's is intra-process: OpenMP threads in one address space, each file read once, sketches shared. This is completely unaffected by everything I found today. It scales, and we can show a real thread-scaling curve.
- bedtools' is workflow-level: N² independent process launches under a GNU parallel wrapper. On this cluster that saturates at ~1.2–2.9× no matter how many cores you supply, because process creation caps near 123 exec/s.
- Concretely, three ways to put parallelization in the figure, all supported by data we already have or can cheaply get:

1. A threads panel — wall or throughput vs thread count at fixed N, both tools, efficiency annotated. hammock rises; bedtools plateaus. This makes the point vividly and is self-documenting. "One thing worth flagging for when you write the caption: this figure makes the mechanism argument, not a "we're faster" argument. The honest reading is that hammock converts cores into throughput and the bedtools pairwise workflow stops doing so almost immediately — attributable to batch-mode absence, not to sketching."
2. CPU time beside wall — already recorded in every CSV. A reader can see cores actually used, which is precisely the subtraction a skeptical referee makes.
3. Anchor speedups on bedtools t=1 as well, so the headline doesn't depend on how well the wrapper happened to parallelize on a given node.

Mixed-stride should be presented as a methodological contribution within the interval-mode scaling result, not as an unrelated implementation optimization. The main text should briefly contrast it with hash-threshold subsampling; the full comparison among mixed-stride, hash-threshold, and single-hash strategies can remain supplementary.

- **Synthetic scaling result (Figure 3A):** each benchmark configuration comprises two disjoint collections of N files compared as a full cross product, so the number of comparisons grows as N^{2}, reaching 262,144 pairs at N=512. Because hammock constructs one reusable sketch per file while the BEDTools workflow repeatedly processes the underlying interval files, the performance advantage of sketch reuse increases with collection size.
- **Sketch-comparison output cost (Figure 3A):** ~~the lower curves separate the three-column sketch-comparison output used for timing from the analysis-oriented `+IE` output... at N=512, `+IE` increases total wall time by only 1.45%... The two total-time curves are therefore indistinguishable, and only the three-column total is drawn.~~ ~~Corrected 2026-08-09, numbers refreshed 2026-08-10. First: that 1.45%/"indistinguishable" figure was measured at p=14... Second: both total-wall curves are now actually plotted (not just the three-column one)... labelled "hammock total" and "hammock total (+IE)"... with the N=512 annotation stating both speedups over BEDTools (9.2× / 8.4×).~~ **Superseded again, same day (figure simplified).** That two-curve design was itself an intermediate state. As of the 2026-08-10 simplification, Figure 3A plots `jaccard_similarity_ie` exclusively — there is only one hammock curve ("hammock total") and one N=512 annotation (8.4×). The register-equality curve and the paired 9.2×/8.4× annotation moved to Supplementary Figure S9, not deleted. The comparison-phase-cost numbers themselves are unaffected by the simplification (1.64–1.94× at N≥64, total wall +10.0% at N=512 vs +0.7% at N=32) since they describe the +IE arm, which is what's still plotted.
- **Real-data timing result (Figure 3B):** ~~on the 20 Maurano fetal-tissue DNase hypersensitivity BED files, hammock is faster than the parallelized BEDTools workflow even without subsampling.~~ ~~Superseded 2026-08-09: parity or slightly behind, 1.02× at 16 threads and 0.92× at 8 (job 29651773).~~ ~~Superseded again, same day, by the fully-fixed harness: `docs/data/maurano_bedtools.csv` + `maurano_subB_summary.csv`, unsubsampled hammock is 1.16× slower than BEDTools on this corpus... subB=0.1 is 1.58× faster, subB=0.01 is 2.49× faster.~~ **Superseded again, 2026-08-10 (figure simplified):** those 1.16×/1.58×/2.49× numbers are the register-equality (`--no-metrics`) arm (`maurano_subB_summary.csv`), which Panel B no longer plots — they are preserved in Supplementary Figure S9 instead. The current Panel B reads the `jaccard_similarity_ie` arm (`maurano_subB_ie_summary.csv`, costlier since it computes the full metrics block): unsubsampled hammock is **1.33× slower** than BEDTools, subB=0.1 is **1.41× faster**, subB=0.01 is **2.04× faster**. Still the expected shape — N=20 is far below the N≈64 crossover Panel A shows.
- **Mixed-stride result (Figure 3B):** `subB=0.1` and `subB=0.01` further reduce runtime while producing little change relative to exact BEDTools truth (not, as an earlier draft of this bullet said, relative to hammock's own unsubsampled similarity estimates — see the next bullet).
- **Definition of the plotted approximation change (Figure 3B):** ~~\Delta J is defined per pair as the difference between a subsampled run's `jaccard_similarity` and the mean `jaccard_similarity` from hammock's own unsubsampled (`subB=1.0`) runs on that same pair. Mean |\Delta J| is the mean absolute value over the 190 pairs.~~ **Superseded 2026-08-10 (figure simplified):** Panel B's accuracy label is no longer mean |ΔJ| against hammock's own baseline. It is now the MAE of `jaccard_similarity_ie` against **exact BEDTools** Jaccard, over the same 190 unique off-diagonal pairs (self-comparisons excluded). The ΔJ-against-own-baseline framing above, and its underlying `jaccard_similarity` (register-equality) numbers, are preserved as-is in Supplementary Figure S9 rather than being recomputed for the new metric; see that figure's caption for the corresponding definition. (An earlier draft reported the pair count as "950 pair-by-replicate comparisons (190 pairs \times 5 replicates)": the five replicates are byte-identical — hammock is deterministic given the seed — so the effective n is 190, not 950; this correction still applies to whichever quantity is being averaged.)

~~**NOTE: this figure is likely to change.**~~ **Resolved 2026-08-10.** Both panels are now measured on the LD_LIBRARY_PATH-fixed, rotated-arm harness (see `docs/bedtools-baseline-retraction.md`'s "Re-measured and gated" section for the full gate history superseding the process-creation-cap story this note used to tell), and Panel A's small-N cells (N≤32) — which an n=3 draw had gotten backwards, showing hammock faster than bedtools at N=2 — are now corrected by a same-node, 20-replicate rerun (job 29671317). bedtools reliably wins below N≈64; the crossover falls between N=32 and N=64.

![Figure 3](figures/pairwise_scaling.png)

**Figure 3. Hammock expands feasible all-pairs comparison as interval collections grow.** (A) Wall time for hammock and BEDTools across synthetic collections containing 10,000 intervals per BED file, using HyperLogLog precision p=18 (the hammock CLI default — corrected from p=14, used nowhere else in the paper), 16 threads, and three (twenty at N ≤ 32) runs per configuration; N files per side gives N^{2} ordered pairs, 262,144 at N=512. Curves show total hammock wall time on `jaccard_similarity_ie` — **the column we recommend reading for any comparison of magnitude, and the only one this figure reports** — plus its sketch-comparison-phase component in isolation ("hammock sketch comparison"). BEDTools is faster than hammock below N≈64; the two total-wall curves (BEDTools and hammock) cross between N=32 and N=64. At N=512, hammock is 8.4× faster than BEDTools. (The register-equality (`jaccard_similarity`) variant of this panel, plotted alongside `jaccard_similarity_ie` rather than in its place, is Supplementary Figure S9; there the equivalent N=512 numbers are 9.2×/8.4×.) (B) Wall time on 20 Maurano fetal-tissue DNase hypersensitivity BED files using interval mode with p=18 and eight threads, at the `jaccard_similarity_ie` arm (the full metrics block — costlier than a register-equality-only run, see Figure S9). Bars compare BEDTools with hammock at `subB=1.0`, `0.1`, and `0.01`; labels report wall time, speedup relative to BEDTools, and MAE of `jaccard_similarity_ie` relative to exact BEDTools Jaccard over the 190 unique off-diagonal pairs (self-comparisons excluded).

Together, the synthetic and real-data benchmarks show that hammock improves the feasibility of exhaustive pairwise analysis through both reusable sketches and efficient optional subsampling during sketch construction.

**Figure 3 generation.** The combined figure is generated by `paper/pairwise_scaling/plot_pairwise_scaling.R` and written to `paper/figures/pairwise_scaling.png`. Panel A uses `docs/data/cpp_vs_bedtools_t16_p18.csv`, job **29671317**, node **c529**, one exclusive allocation, two passes merged (N∈{2,4,8,16,32} at 20 replicates, N∈{64,128,256,512} at 3) — supersedes job 29656140 (node c531), whose N=2/N=4 cells were noise-corrupted at n=3, which itself replaced the p=14 `cpp_vs_bedtools_t16_20260804_172242.csv`/job 29552415 cited in an earlier draft of this note; Panel B uses `docs/data/maurano_subB_ie_summary.csv` (median wall time and MAE-vs-exact-BEDTools of the `+IE` arm at each subB level, 5 replicates, p=18, t=8) and `docs/data/maurano_bedtools.csv`. **Do not compare Panel A's absolute times to any earlier run.** The p=14 predecessor (job 29552415, node sr14) was itself a genuine allocated SLURM job; the "1.27–1.53× slower, no SLURM allocation" finding belongs to a *different*, unrelated 2026-08-04 file (`pairwise_cost_by_precision_20260804_164807.csv`, the metrics-block/fusion-cost experiment) compared against its 08-08 successor — not to job 29552415 (see `docs/seed-benchmark-methodology.md`, which distinguishes the two). The two share a date but are otherwise unconnected; conflating them was an error in an earlier draft of this note. Regardless of mechanism, do not compare absolutes across any of these runs — within-run ratios survive cross-run confounds; absolutes do not. The script rebuilds both panels with harmonized typography, panel labels, tool colors, wall-time units, and the pair count N^{2}.

#### 2.2.2 Interval sketches recover exact-overlap similarity structure and values

![Figure 4](figures/interval_accuracy.png)

**Figure 4. Inclusion–exclusion estimates closely reproduce exact BEDTools Jaccard.** Pairwise comparisons were performed on 20 Maurano fetal-tissue DNase hypersensitivity BED files, yielding 190 unique off-diagonal pairs; self-comparisons were excluded. Hammock's `jaccard_similarity_ie` uses HyperLogLog cardinality estimates with inclusion–exclusion to approximate the same base-pair Jaccard statistic computed exactly by BEDTools. At HyperLogLog precision p=21, the estimated values closely track exact BEDTools Jaccard over the observed off-diagonal range (MAE =4.3\times10^{-4}, Pearson r=0.99999, Kendall \tau=0.9947). Across precisions p=18, p=21, and p=23, absolute deviation from BEDTools decreases as precision increases, from approximately 1\times10^{-3} to approximately 2\times10^{-4}. The figure therefore evaluates the numerical accuracy of the inclusion–exclusion estimate relative to exact interval Jaccard; it does not compare inclusion–exclusion with hammock's register-equality statistic.

**Figure 4 generation.** The figure is generated by `paper/interval_accuracy/plot_interval_accuracy.R` and written to `paper/figures/interval_accuracy.png`. Inputs are `docs/data/maurano_bedtools_ref.tsv` and `docs/data/hammock_hll_p{18,21,23}_jaccB.csv`. Those three CSVs were regenerated with hammock 0.5.0 on 2026-07-31 to obtain the `jaccard_similarity_ie` column; earlier copies carried the placeholder `containment` column and could not supply it. The regeneration reproduces the archived `jaccard_similarity` values **exactly** (0 of 400 pairs differ at any precision), so panel A's register-equality series is unchanged from the previously published figure.

### 2.3 Sequence Space
#### 2.3.1 Biological identity is preserved across references in sequence space

![Figure 5](figures/cross_reference_identity.png)

**Figure 5. Sequence sketches group samples by tissue across genome references.** Three ENCODE H3K27ac ChIP-seq samples—heart ENCSR175ABH, liver ENCSR864OOO, and lung ENCSR954JMZ—were independently aligned and peak-called against GRCh37, GRCh38, and CHM13, yielding nine peak sets. Each peak set was converted to sequence against its native reference and sketched in sequence mode using broad peaks, k=10, w=10, HyperLogLog precision p=24, and the minimizer-only Jaccard metric. UPGMA clustering on 1-J groups the nine peak sets by tissue rather than by reference genome: heart, liver, and lung each form a distinct clade containing the GRCh37-, GRCh38-, and CHM13-derived representation of that tissue. In this proof-of-concept panel, reference choice therefore behaves as a within-tissue perturbation rather than as the dominant axis of separation, showing that sequence sketches can retain biological identity across genome references without requiring direct coordinate overlap.

**Figure 5 generation.** The figure is generated by `paper/cross_reference_identity/plot_cross_reference_identity.R` and written to `paper/figures/cross_reference_identity.png`. Inputs are `docs/data/exp_a_broad_k10_w10.csv` (staged from `experiments/ref-comparison/results/exp_a/broad/k10_w10/exp_a_mnmzr_p24_jaccD_k10_w10.csv`) and `docs/data/exp_a_metadata.tsv` (staged from the experiment configuration). The script constructs the nine-sample similarity matrix, converts similarity to 1-J, performs average-linkage hierarchical clustering, and renders the tissue-colored dendrogram.

### 2.3.2 Sequence sketches recover tissue organization

- **Quantitative clustering result for the Results prose:** cutting the average-linkage dendrogram into the ten annotated tissue classes gives adjusted Rand index (ARI) = 0.9102 and normalized mutual information (NMI) = 0.9610. These values quantify a separate 10-cluster evaluation and should not be presented as properties of the visual tissue annotations in Figure 6.
- **Figure interpretation:** Figure 6 itself does not display a discrete cluster cut. Its labels are colored by the known tissue annotations to show how biological identity is organized along the dendrogram.

![Figure 6](figures/sequence_tissue_clustering.png)

**Figure 6. Sequence sketches recover fetal-tissue organization in the Maurano DNase hypersensitivity collection.** Twenty fetal-tissue DNase-seq interval sets spanning ten annotated tissue labels were converted to interval-derived sequence and compared using the minimizer-only `jaccard_similarity` score at k=10, w=30, and HyperLogLog precision p=24. Average-linkage hierarchical clustering was performed on 1-J. Leaf labels show dataset accessions and are colored by annotated tissue. The tissue coloring is an annotation of the dendrogram and does not represent a discrete cluster cut.

**Figure 6 generation.** The figure is generated by `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` and written to `paper/figures/sequence_tissue_clustering.png`. By default, the script reads `experiments/maurano_dhs_validation/results/raw_d/hammock_mnmzr_p24_jaccD_k10_w30.csv`, uses the `jaccard_similarity` column, and obtains tissue annotations from `experiments/maurano_dhs_validation/data/maurano_filenames_key.tsv`. Optional command-line arguments can substitute another similarity CSV, tissue key, output path, or similarity column. Passing `jaccard_similarity_ie` as the fourth argument produces Supplementary Figure S5a; the column is derived from the containment columns the archived CSV already carries, since the sweep predates it. Two further renderings of the same configuration at different precisions, `paper/figures/sequence_tissue_clustering_p12.png` and `sequence_tissue_clustering_p24.png`, record that the clustering is unchanged across the p ≥ 12 plateau and are retained as supporting material rather than manuscript figures.

### 2.3.3 Parameterization separates numerical agreement from biological resolution

- **Figure 6 parameter choice for the Results prose:** k=10, w=30 is the ARI-optimal sequence-mode setting under both similarity estimators at p=24; the parameter sweep, rather than Figure 6 itself, is the evidence for that choice.

![Figure 7](figures/parameter_objective_tradeoff_estimators.png)

**Figure 7. Sequence parameters determine biological recovery, while similarity estimator choice primarily changes numerical calibration.** Sequence-mode parameter configurations were evaluated on 20 Maurano fetal-tissue DNase hypersensitivity interval sets at HyperLogLog precision p=24 using the same minimizer-derived HLL sketches but two alternative similarity estimators. Panel A shows Hammock's register-equality statistic (`jaccard_similarity`), and Panel B shows the inclusion–exclusion Jaccard estimate (`jaccard_similarity_ie`). Each point represents one combination of k-mer size (k) and minimizer-window size (w); the x-axis gives Pearson correlation with exact BEDTools Jaccard across pairwise comparisons, and the y-axis gives adjusted Rand index (ARI) for recovery of the ten annotated tissue classes after clustering. Point color encodes k and point size encodes w. In each panel, the estimator-specific configuration with greatest numerical agreement is emphasized by an orange-red outline, and the shared k=10, w=30 biological optimum used in Figure 6 is marked by a hollow orange-red diamond; dotted leader lines connect these highlighted points to their callouts. The numerical optimum shifts from k=20, w=100 for register equality (r=0.9996) to k=20, w=20 for inclusion–exclusion (r=0.99997), whereas the biological optimum remains k=10, w=30 with ARI=0.9102 under both estimators. At p=24, the estimators give the same ARI at 38 of 39 parameter cells. The lone exception, k=8, w=500, is identified by a grey outline and is biologically uninformative under both estimators (ARI 0.059 versus 0.025). At the Figure 6 setting, the two estimators produce the identical ten-cluster partition, not merely the same summary score. Thus estimator choice changes numerical calibration and the location of the BEDTools-matching optimum, while the biologically useful parameter region and the selected tissue-recovery optimum are stable at the precision used in the manuscript.

**Figure 7 generation.** The paired figure is generated by `paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R` and written to `paper/figures/parameter_objective_tradeoff_estimators.png`. The existing `plot_parameter_objective_tradeoff.R` is retained unchanged as the single-estimator diagnostic workflow. The new script reads `docs/data/mode_d_summary.csv`, supplementing from `experiments/maurano_dhs_validation/results/mode_d_summary.csv` if the staged paper copy lacks a current estimator, and filters to p=24, BEDTools reference, k≥8, and non-missing Pearson/ARI values. It computes the estimator-specific Pearson optima, verifies the k=10, w=30 ARI optimum under both estimators, and derives the p=24 ARI-disagreement count directly from matched (k,w) cells.

## III. Discussion

### 3.1 Scaling comparison within references

Discuss both layers of the computational contribution: reusable sketches reduce the cost of repeated all-pairs comparisons, while mixed-stride makes optional base-pair subsampling computationally effective by avoiding a per-position hashing gate. Mixed-stride should be described as a novel deterministic subsampling strategy with chromosome-specific phase offsets, not merely as a command-line tuning option.

### 3.2 Comparing interval-derived sequence across references

### 3.3 Coordinate similarity and sequence similarity are complementary

### 3.4 Practical recommendations

### 3.5 Limitations and future work

## IV. Methods

### 4.1 Software implementation

\program{hammock} is a command-line tool and importable Python package with a C++17 extension that uses OpenMP for parallel execution. It accepts two positional inputs, each a plain-text file containing paths to a collection of genomic datasets to be compared. Hammock selects interval or sequence mode based on the input dataset formats and supplied reference arguments unless a mode is explicitly requested. In all modes, one reusable sketch is constructed for each listed dataset, and pairwise comparisons between the two collections are written to a comma-separated table.

Hammock implements four representations that can be selected with `--mode`.

| CLI mode                       | Input                                | Representation                                       |
| ------------------------------ | ------------------------------------ | ---------------------------------------------------- |
| `interval-string`              | BED                                  | Exact chromosome–start–end interval strings          |
| `interval` / `interval-points` | BED                                  | Base-level genomic positions covered by intervals    |
| `interval-hybrid`              | BED                                  | Combined interval-string and position representation |
| `sequence`                     | FASTA or BED with a reference genome | Sliding-window minimizers derived from sequence      |

This study evaluates the base-level `interval` representation for within-reference comparisons and the `sequence` representation for comparisons of interval-derived sequence. The interval-string and interval-hybrid representations remain available in the software but are not evaluated as primary methods here.

For sequence comparisons beginning from BED files, hammock extracts the nucleotide sequence underlying each interval using the reference genome associated with that collection. `hammock fetch-ref` can be used to facilitate the download of references. Sequence extraction produces one FASTA file per BED file, which is then processed using the same sequence-mode workflow as directly supplied FASTA input.

For each pair of input datasets, hammock reports the following similarity summaries:

| Output field            | Interpretation                                                                 |
| ----------------------- | ------------------------------------------------------------------------------ |
| `jaccard_similarity_ie` | Set Jaccard estimated from HyperLogLog cardinalities using inclusion–exclusion |
| `jaccard_similarity`    | Fraction of active HyperLogLog registers with equal values                     |
| `containment_AB`        | Estimated fraction of dataset A contained in dataset B                         |
| `containment_BA`        | Estimated fraction of dataset B contained in dataset A                         |
| `cosketch_geom`         | Geometric mean of the two directional containment estimates                    |
| `cosketch_arith`        | Arithmetic mean of the two directional containment estimates                   |
| `cosketch_max`          | Maximum of the two directional containment estimates                           |

The analyses in this study use the register-equality Jaccard estimate or the inclusion–exclusion Jaccard estimate when comparison with exact set Jaccard is required. Sketch construction and interval subsampling are reproducible for fixed parameter values and seeds.

### 4.2 Interval-mode sketch construction

#### 4.2.1 Hyperloglog Sketching

#### 4.2.2 Mixed-stride subsampling

Interval mode passes the points of a base-pair expansion of each interval to a HLL sketch. There is a certain computational excess involved in adding every point in an interval, so \program{hammock} includes a \texttt{--subB} command that allows for a subsampling of the points to be added to the sketch, where the value of \texttt{subB} is a proportion between $0$ and $1$. In order to maintain reproducible randomness of point selection, mixed-stride subsampling is used.

In mixed-stride subsampling, a deterministic chromosome-keyed offset is used to select the sampling phase, interpolating between two strides to exactly match a target probability $p$, with deterministic reproducibility.

ADD MATH HERE

Define `subB`, the sampled base-pair universe, and the mixed-stride algorithm. For sampling fraction `subB`, set the approximate stride to s=\mathrm{round}(1/\mathrm{subB}), use a deterministic chromosome-keyed offset to select the sampling phase, and advance directly between accepted positions rather than evaluating every covered position.

Explain determinism, seed handling, expected sampling density, computational cost, and the rationale for chromosome-specific offsets. Contrast the algorithm with hash-threshold subsampling, whose inclusion test requires a hash computation for every covered position regardless of the requested sampling fraction.

### 4.3 Sequence-mode sketch construction

#### 4.3.1 Minimizers

### 4.4 Similarity estimators and computational complexity

*Draft outline. Key sentences are given; connecting prose is not written.*

- **Framing.** Both estimators are computed from the same pair of sketches. *"Hammock's two similarity columns differ in what they estimate and in what they cost, not in the data they see."*
- **Register-equality (**`jaccard_similarity`**).** Define as the fraction of active registers whose values agree.
  - *"This is not set Jaccard: two registers agree whenever the largest ρ observed in that bucket happens to coincide, which occurs at some rate even for disjoint inputs, so the statistic carries a chance-agreement floor c."*
  - c is set by the sketch load factor λ = n/m **and** by the cardinality ratio |A|/|B| — not by precision as such. Report 0.1699 at equal cardinality as λ→∞, 0.058 at ratio 10, 0.045 at p = 24 once m > n.
  - *"Because c differs between pairs of differing size, register-equality is order-preserving within a cardinality ratio but not across ratios."* Ties to the 2.45% Maurano inversions in §2.3.
- **Inclusion–exclusion (**`jaccard_similarity_ie`**).** Define |A ∩ B| = |A| + |B| − |A ∪ B| with each term an Ertl (2017) estimate, union by register-wise maximum, intersection clamped to ≥ 0; J = |A∩B|/|A∪B|, equivalently 1/(1/C_AB + 1/C_BA − 1).
  - *"An exact 0.0 in this column means the intersection estimate hit the non-negativity clamp or the inputs are genuinely disjoint; it never means a measured zero."*
  - Relative error ≈ 0.6/(J√m), hence uninformative below J ≈ a few/√m — state as the column's domain of validity.
- **Cost.**
  - *"Register-equality is a ratio of two register counts and costs a single pass over the register arrays; inclusion–exclusion additionally requires a union sketch and cardinality estimates per pair, so its cost scales as O(N²·2^p) against O(N·M) for sketch construction."*
  - Measured overhead, same job and binary, 16 threads, p = 14 (`docs/data/cpp_vs_bedtools_t16_20260804_172242.csv`): the comparison *phase* costs 2.9–3.8× more for N ≥ 32, but that phase is only 0.4–0.5% of wall time, so **total wall rises 1.45% at N = 512** (54.15 s against 53.38 s) and stays below 0.6% at every smaller N. Peak RSS 34.8 MB against 22.7 MB.
    - *(Corrected 2026-08-09. This read "not measurable at N = 20, +2.4% at N = 100, +10.4% at N = 512". None of the three was sourced — that CSV has no N=20 or N=100 row — and the N=512 value contradicted the 1.45% recorded above under "Sketch-comparison output cost" and in `docs/figure3-panel-a-rebuild.md`. The old number is consistent with a pre-fusion measurement: the 2026-08-08 fused register pass bought 2.2–5.4× on exactly this path.)*
  - *"The overhead is therefore negligible in wall time across the entire measured range, and the choice between the two columns can be made on statistical grounds rather than cost."* Note the phase ratio does grow with N, so the phase must eventually dominate; it has not done so by N = 512. (This bullet previously concluded the overhead was "material at the catalog scale", which followed from the +10.4% number corrected above.)
  - ~~State explicitly which configuration Figure 3 benchmarks, so the speed claim and the reporting recommendation refer to the same thing.~~ **Done, 2026-08-09**: see the Figure 3 caption above.
- **Why not the joint maximum-likelihood estimator.** *"Ertl's joint MLE is the lower-variance choice for HyperLogLog Jaccard, but it requires the joint register histogram of the pair, which hammock's sketch interface does not expose; the inclusion–exclusion form is recoverable from the containment estimates already computed and needs no additional interface."* Pre-empts the obvious reviewer question.
- **Which column to report.** *"We report* `jaccard_similarity_ie` *for any comparison of magnitude, and for any comparison spanning different set sizes;* `jaccard_similarity` *is retained as the low-cost default."*
  - One clause on the exception (register-equality ranks better below J ≈ 0.05 at p ≤ 20 among comparably-sized pairs, and one step of precision removes the advantage), with the quantification deferred to a supplementary note rather than carried in Methods.
  - **Figure 6 estimator note:** the current Figure 6 is generated from `jaccard_similarity` because that is the column emitted by the original sequence-mode sweep. At Figure 6's configuration (k=10, w=30, p=24) the two columns induce the *identical* 10-cluster partition — not merely equal scores, the same member sets — with ARI and NMI agreeing to sixteen digits. The equivalence is **local, not global**: across the 235-cell sweep the two columns give different ARI at 48 cells (max |ΔARI| 0.301). The disagreement decays with precision — 16 of 17 cells differ at p=10, but only 1 of 39 at p=24, and that one (k=8, w=500) has ARI 0.059 against 0.025, i.e. both clusterings are degenerate. So the defensible statement is "identical at p=24 except one degenerate cell, and at every cell the manuscript cites", not "identical everywhere". Generated by `paper/sequence_tissue_clustering/estimator_agreement_stats.csv` and `experiments/maurano_dhs_validation/results/estimator_ari_by_config.csv`. This is methodological context and should not be carried in the Figure 6 caption.
- **Complexity.** Retain the original stub: cost of mixed-stride sketch construction as a function of covered length and stride, alongside full interval ingestion and pairwise sketch comparison.

### 4.5 Datasets

### 4.6 Performance benchmarking

Include the direct comparison of no subsampling and mixed-stride `subB = 0.1` and `0.01` used in Figure 3B. The main benchmark establishes the practical speed–approximation trade-off; exhaustive comparison against hash-threshold and single-hash alternatives can be reported in Supplementary Results.

### 4.7 Interval-mode accuracy evaluation

### 4.8 Sequence-mode biological evaluation

## Data and code availability

## Author contributions

## Acknowledgments

## References

## Supplementary Methods

### S1. Quantification of public interval collections

### S2. Alternative `subB` sampling strategies and implementation details

Document hash-threshold and single-hash alternatives, the full mixed-stride derivation and reproducibility checks not needed in the main text, and any structured-sampling sensitivity analyses.

## Supplementary Results

### S3. Comparison of mixed-stride, hash-threshold, and single-hash subsampling

Report runtime, similarity deviation, and scaling across sampling fractions and file sizes. This supports the mixed-stride contribution without requiring an additional main-text figure unless its comparative advantage becomes a headline result.

### S4. Supplementary note for Section 2.3: register-equality and inclusion–exclusion estimate different quantities

Figure 4 evaluates only the numerical accuracy of `jaccard_similarity_ie` relative to exact BEDTools Jaccard. Hammock also emits `jaccard_similarity`, a register-equality statistic defined as the fraction of active HyperLogLog registers whose values agree. This statistic is not set Jaccard: registers can tie by chance, producing a positive chance-agreement floor whose magnitude depends on sketch load and on the cardinality ratio |A|/|B|.

In separate analyses on the same 20-file Maurano corpus, register-equality preserved the broad ordering of pairwise similarities but did not reproduce exact Jaccard magnitudes. At p=21, register-equality had MAE =0.1378, Pearson r=0.99720, and Kendall \tau=0.9511 relative to BEDTools Jaccard, whereas inclusion–exclusion had MAE =4.3\times10^{-4}, Pearson r=0.99999, and Kendall \tau=0.9947. The register-equality offset was approximately 0.14 and remained largely unchanged across p=18, p=21, and p=23, consistent with an estimator-specific chance-agreement floor rather than ordinary sketch sampling error. By contrast, inclusion–exclusion error decreased with precision.

The ordering differences under register-equality were small but systematic. Kendall \tau=0.951 corresponds to 439 of 17,955 pair-of-pair comparisons (2.45%) in which register-equality and BEDTools disagreed on ordering; all involved BEDTools Jaccard differences below 0.025. The same residual pattern recurred across independently generated sketch sets at p=18 and p=21 (r=0.994), consistent with the chance floor changing with the relative cardinalities of the compared interval sets. Under inclusion–exclusion at p=21, the same corpus produced 48 ordering inversions (0.27%). These analyses explain why the main-text Figure 4 focuses on inclusion–exclusion when the goal is to recapitulate exact Jaccard values.

### S5. Sequence-mode results are unchanged when read on `jaccard_similarity_ie`

Every sequence-mode figure (Figures 5, 6, 7) is drawn on `jaccard_similarity`, because that is the column the sequence-mode sweeps emitted; `jaccard_similarity_ie` shipped afterwards. Because the archived CSVs carry `containment_AB` and `containment_BA`, the inclusion–exclusion column is exactly recoverable from them (`J = 1/(1/C_AB + 1/C_BA − 1)`), so this check required no re-sketching.

**Tissue clustering (Figure 6).** At k=10, w=30, p=24 the two columns induce the identical ten-cluster partition (ARI = 0.9102, NMI = 0.9610 under both, agreeing to sixteen digits). Supplementary Figure S5a (`paper/figures/sequence_tissue_clustering_k10_w30_p24_ie.png`) is the Figure 6 dendrogram redrawn on `jaccard_similarity_ie`; every annotated tissue occupies a contiguous, monophyletic block, as in Figure 6. Across the full 235-cell (k, w, p) sweep the columns are *not* interchangeable — ARI differs at 48 cells, maximum |ΔARI| = 0.301 — but the disagreement is concentrated at low precision (16 of 17 cells at p=10 against 1 of 39 at p=24) and does not touch any configuration the manuscript reports.

**Parameter selection (Figure 7).** The ARI-optimal configuration is k=10, w=30 under both columns, at the same ARI of 0.9102, and the best-ARI value is reached on a plateau spanning every precision p ≥ 12. The Pearson-optimal configuration differs between columns — k=20, w=100 on `jaccard_similarity` (r = 0.9996) against k=20, w=20 on `jaccard_similarity_ie` (r = 0.99997) — so the two panels of Figure 7 annotate different best-Pearson points. Figure 7's claim is unaffected: under either column the Pearson optimum has ARI 0.693 against the ARI optimum's 0.910, so numerical agreement and biological resolution still select different settings.

> **Caveat to resolve before publishing the Pearson optima (added 2026-08-09).** The BEDTools reference behind these correlations was computed with **bedtools 2.27.1**, which has an **order-dependent union**: 93 of the 190 unordered Maurano pairs had J(A,B) ≠ J(B,A) in what is supposed to be the exact reference. The magnitude is small (max |Δ| = 1.3×10⁻⁵, mean 5.6×10⁻⁷), so no clustering conclusion moves. But the r = 0.9996 versus r = 0.99997 distinction above lives at 1−r ≈ 4×10⁻⁴ against 3×10⁻⁵, which is within an order of magnitude of the reference's own inconsistency — so the *identity* of the Pearson optimum, and the claim that it shifts between columns, is not yet safe to assert at that resolution. `bedtools` is now pinned to 2.30.0; re-derive the reference before this paragraph goes in the manuscript. The ARI results are unaffected, being rank-and-partition quantities. See CLAUDE.md, benchmark-methodology seed.

**Cross-reference identity (Figure 5).** Redrawn on `jaccard_similarity_ie` (`paper/figures/cross_reference_identity_ie.png`), the nine peak sets still form three monophyletic tissue clades, each containing all three reference-genome representations. Ranking the twenty (k, w) cells by the separation statistic Δ = median(same-tissue cross-reference) − median(different-tissue) gives Spearman ρ = 0.9925 between the two columns on both broad and narrow peaks, with rank changes confined to the top five cells, which lie within 0.06 of one another. The k ≥ 15 and k ≤ 10 groups remain disjoint under both columns, with a slightly wider gap on `jaccard_similarity_ie` (broad +0.415 against +0.389).

Caveat worth carrying: `jaccard_similarity_ie` is censored at zero by the non-negativity clamp on the inclusion–exclusion intersection, so an exact 0.0 means "clamped or empty", never "measured zero". This does not arise at the configurations above but would in a low-similarity corpus.

## Supplementary Figures and Tables

- **Figure S5a — sequence-mode tissue dendrogram on `jaccard_similarity_ie`:** `paper/figures/sequence_tissue_clustering_k10_w30_p24_ie.png`, generated by `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` with `jaccard_similarity_ie` as the fourth positional argument. Companion to Figure 6; see S5.
- **Figure S5b — cross-reference dendrogram on `jaccard_similarity_ie`:** `paper/figures/cross_reference_identity_ie.png`, generated by `paper/cross_reference_identity/plot_cross_reference_identity.R` with the output path and `jaccard_similarity_ie` as its two positional arguments. Companion to Figure 5; see S5.
- **Table S5 — estimator agreement:** `paper/sequence_tissue_clustering/estimator_agreement_stats.csv` (per-column ARI/NMI and partition identity at Figure 6's cell), `experiments/maurano_dhs_validation/results/estimator_ari_by_config.csv` (all 235 cells), `experiments/maurano_dhs_validation/results/best_config_by_column.csv` (per-column optimum with tie sets), and `experiments/ref-comparison/results/exp_a_estimator_delta.csv` (cross-reference Δ per cell).
- **TODO — sequence-mode parameter heatmaps:** include the full k × w parameter-response heatmaps for numerical agreement with BEDTools and biological clustering performance (ARI/NMI), potentially stratified by HyperLogLog precision and/or estimator. These heatmaps should provide the complete response surface underlying the compact objective-space summary in Figure 7 without burdening the main-text figure.

### Figure S9 — Figure 3 companion: register-equality vs +IE, both scored against exact BEDTools

<img src="figures/pairwise_scaling_supplement.png" alt="Figure S9" width="850">

**Figure S9.** Companion to Figure 3, for the supplement / advisor review rather than the main text. Same two panels and same underlying runs as Figure 3, but Panel A keeps both hammock variants (`jaccard_similarity` and `jaccard_similarity_ie`, i.e. "hammock total" and "hammock total (+IE)", plus their sketch-comparison-phase components), and Panel B's register-equality bars are now scored the same way as the `+IE` bars — mean absolute error against exact BEDTools Jaccard, not drift from hammock's own unsubsampled run (which is what a pre-2026-08-10 version of this panel reported, and which is tautologically small at `subB=1.0` since that run *is* the baseline it would be diffed against). On that shared footing, register-equality's MAE is ~0.137–0.138 across all three subB levels — roughly 100× the `+IE` MAE (0.0012–0.0016) — and it is essentially flat across subsampling levels, because it is dominated by the register-equality chance-agreement floor (CLAUDE.md divergence #2, `c` ≈ 0.17), not by anything subsampling does. This is the reason Figure 3 itself reports `+IE` only: putting a bedtools-uncalibrated column next to a bedtools-calibrated one at face value, without this framing, invites exactly the miscomparison this figure is designed to make legible instead.

**Figure S9 generation.** `paper/pairwise_scaling/plot_pairwise_scaling_supplement.R`, forked 2026-08-10 from the script that produces Figure 3 at the commit where it still plotted both hammock variants. Panel A shares Figure 3's `docs/data/cpp_vs_bedtools_t16_p18.csv` (job 29671317). Panel B shares Figure 3's `docs/data/maurano_bedtools.csv` and `docs/data/maurano_subB_ie_summary.csv`, plus a new `docs/data/maurano_subB_re_vs_bedtools.csv` (register-equality MAE vs exact BEDTools, computed from `experiments/subB_mixed_stride/results/sweep_maurano_ie_20260809_200658.csv` — mixed-stride, p=18, t=8, `rep=0`, against `docs/data/maurano_bedtools_ref.tsv`; reps are byte-identical, see `experiments/subB_mixed_stride/RESULTS.md` "Harness validation") and the existing register-equality wall times in `docs/data/maurano_subB_summary.csv`.

### Figure S6 — Thread scaling

<img src="figures/threading_supplement.png" alt="Figure S6" width="850">

**Figure S6.** Wall time, speedup vs each tool's own *t*=1, and parallel efficiency across thread counts 1–48, on a synthetic 64-files-per-side corpus (10k intervals/file, p=18, subB=1.0, 3 replicates, means). Deliberately separate from Figure 3, which is about cost vs. collection size; this is about what each tool does with the cores it is given, and shares no axis with Figure 3. BEDTools plateaus by t=8 (4.22× over its own t=1) and stays flat through t=48 (4.20×), converting at most 6.5 of its available cores into throughput because `bedtools jaccard` has no batch mode — a pairwise workflow is N² separate process launches, each re-reading its input files, wrapped in GNU `parallel`. hammock keeps improving to t=32 (12.11×, 24.0 of 32 cores busy) before easing back slightly at t=48 (11.59×, 28.1 of 48 — an oversubscription effect, not a new mechanism; see `docs/seed-mode-d-threading.md` for the analogous Mode D finding). Every BEDTools row in the benchmark CSVs carries `mean_bedtools_parallel_eff`, so this achieved-parallelism gap is checkable rather than assumed on any figure that quotes a BEDTools speedup at a fixed thread count. Caveat: one node, one corpus size, one precision — the *level* of BEDTools' achieved parallelism is node-dependent (measured 1.17–2.86× across three otherwise-identical runs), but the *shape* (early saturation, then a plateau rather than a steep regression) is the part to trust as general.

**Figure S6 generation.** `paper/pairwise_scaling/plot_threading_supplement.R` from `docs/data/sweep_threads_p18.csv` (job 29670792, node c531, one exclusive allocation).

### Figure S7 — Catalog scale

<img src="figures/largeN_supplement.png" alt="Figure S7" width="850">

**Figure S7.** hammock measured to N=2048 (the ChIP-Atlas hg38 CTCF manifest holds 2,206 verified files, so this is a real corpus size) against a **projected**, not measured, BEDTools curve — drawn dashed, grey, and with open points to keep a projection from being mistaken for a third measurement. Projecting was the only practical option: BEDTools measures 653.4 s at N=512 (the same corrected, job-29671317 measurement Figure 3A plots — do not substitute the retracted pre-LD_LIBRARY_PATH-fix figure of 1978 s from the top-of-document notice), growing ~3.97× per doubling, so one replicate is ~0.7 h at N=1024 and ~2.9 h at N=2048 — roughly 11 h of node time for three replicates at each of N=1024 and N=2048, to extend a comparison Figure 3 already settles by N=256. The fitted exponent is 2.000 (theory 2.0, R²=0.9999, fit on N≥32 only, since smaller N is startup-dominated). At N=512, 1024, and 2048 hammock measures 71.35 s, 164.68 s, and 417.27 s against BEDTools' 653.4 s (measured) and 2574 s / 10294 s (projected), for speedups of 9.2×, 15.6× (projected), and 24.7× (projected). hammock's own curve is not extrapolated — every hammock point shown is measured — and its cost per doubling is itself rising (2.30× then 2.53×) as the Θ(N²) pairwise phase begins to overtake linear sketching, which is why only the BEDTools segment is dashed. The 9-column metrics block costs 1.27× at N=2048 (529.6 s vs 417.3 s) against **~1.15×** at N=20 on Maurano (`maurano_subB_ie_summary.csv` wall 9.404 s vs `maurano_subB_summary.csv` wall 8.188 s at subB=1.0 — not directly readable off Figure 3B, which as of the 2026-08-10 simplification plots only the +IE arm; see Figure S9 for the register-equality comparison point) — the block is a per-*pair* cost, so its share still grows markedly with N, and "the metrics are nearly free" from the small-N comparison does not generalize.

**Figure S7 generation.** `paper/pairwise_scaling/plot_largeN_supplement.R` from `docs/data/cpp_vs_bedtools_t16_p18.csv` (job 29671317, node c529 — same run as Figure 3A) joined to `docs/data/cpp_vs_bedtools_t16_p18_largeN.csv` (job 29652432, node c594, hammock-only, N ∈ {512, 1024, 2048}, same corpus generator and seed schedule as Figure 3A). The script checks and prints two internal consistency gates rather than assuming them: the N=512 cell measured by both jobs on different nodes agrees to 1.1% (70.97 s vs 71.73 s), and the fitted BEDTools exponent (2.000) is reported alongside its R².

### Figure S8 — Precision frontier

<img src="figures/precision_frontier.png" alt="Figure S8" width="850">

**Figure S8.** x = MAE of `jaccard_similarity_ie` vs exact BEDTools (log), y = speedup vs BEDTools (log), one point per HyperLogLog precision p ∈ {12,14,16,18,20,22,24}, on the 20-file Maurano corpus (380 ordered pairs after self-pair exclusion, subB=1.0, 3 replicates, means). Deliberately not twin y-axes against p: twin log axes have no canonical registration, so an apparent crossing would be an artifact of axis placement rather than a property of the data. Plotted this way, iso-accuracy is a vertical line and "slower than BEDTools" is the region below y=1 — which is every precision tested on this corpus (0.77× at p=12 down to 0.08× at p=24). Accuracy improves 57× from p=12 to p=24 (MAE 8.8e-3 → 1.5e-4) while speed falls ~10×, so there is no accuracy/speed knee to trade toward on this corpus — the choice of p here is purely an accuracy decision. This panel must be read together with Figure 3A, not against it: at N=20 (this corpus) hammock is slower than BEDTools at every precision, while Figure 3A reports 8.4× at N=512 — both are true, because Maurano sits well below the N≈64 crossover, so sketching cost has not yet been amortized over enough pairs. `jaccard_similarity` (register-equality) cannot usefully appear on this x-axis: its MAE is 0.1374–0.1383 across the whole precision sweep (varying 1.01× while `_ie`'s MAE moves 57×), because it is dominated by the register-equality chance-agreement floor rather than by sampling noise — the same floor Figure S9 makes visible directly.

**Figure S8 generation.** `paper/pairwise_scaling/plot_precision_frontier.R` from `docs/data/sweep_precision_maurano_p18_t16.csv` (job 29670793, node c432; the t=8 companion `sweep_precision_maurano_p18_t8.csv` from the same job is a cross-check, not a second plotted series — accuracy at t=8 is bit-identical to t=16, since thread count changes only how fast a deterministic sketch is built, not its content, and this was verified exactly: max |ΔMAE| = 0 at every precision). **Only the x-axis (MAE) is guaranteed to match Figure 3B's Maurano bar; the y-axis (speedup) is not**, since this panel's t=16 is a different thread count from Figure 3B's t=8 — the p=18 speedup reads 0.69× here against 0.86× (1.16× slower) on Figure 3B, a real ~25% gap from the thread-count difference, not an inconsistency to resolve.