## Thesis

Modern interval collections are limited by two boundaries: the computational cost of comparing every dataset and the coordinate dependence that prevents comparison across references. Hammock creates reusable sketches of either interval coordinates or interval-derived sequences, expanding comparison within and across those boundaries while preserving useful biological structure.

## Abstract

## I. Introduction

Increased availability of high-throughput sequencing has facilitated a corresponding increase in large-scale modern sequencing projects which produce vast numbers of genomic annotation datasets which then, themselves, become part of the genomic data ecosystem. Large-scale initiatives such as the ENCODE Project\cite{encode_2012}, the Roadmap Epigenomics Project\cite{roadmap_2015}, GTex[REF], and the 1000 Genomes Project\cite{10002015global} make datasets publicly available, providing unprecedented opportunities for integrative analysis. Among the most common and biologically informative outputs of these projects are genomic intervals: contiguous stretches of the genome that define transcription factor binding sites, chromatin accessibility peaks, histone modifications, structural variants, and splicing junctions. Such interval datasets have become a cornerstone of downstream analyses that seek to characterize genome function and variation across populations, tissues, cell types, and experimental conditions. In particular, pairwise comparison of interval datasets presents a measurement of biological similarity, such as identifying conserved regulatory elements between tissues. Tools such as \program{BEDTools} provide exact calculations of these measures\cite{quinlan2010bedtools}, but the computational cost of pairwise overlap calculations grows rapidly with the size of modern repositories. In practice, this creates a scalability bottleneck that hampers systematic comparisons across the interval datasets now available\cite{li2020design}.

The numbers of files in these interval databases continue to grow every year. ChIP-Atlas, for example, expanded from 37,720 accumulated experiments in 2015 to 464,655 in 2025, while repositories such as ENCODE contain thousands of additional ChIP-seq and chromatin-accessibility experiments.

Scalability is not the only obstacle to systematic interval comparison. BED coordinates are meaningful only with respect to the reference genome on which they were defined, and large public collections remain distributed across multiple genome assemblies. Standard overlap measures therefore cannot directly compare two interval files when one is defined on hg19 and the other on hg38, even when the files describe the same assay or biological feature. Figure 1 illustrates this limitation for histone-mark ChIP-seq peak files from Roadmap Epigenomics and BLUEPRINT. For histone marks represented in both collections, within-Roadmap and within-BLUEPRINT comparisons can be performed in their respective coordinate systems, whereas every Roadmap–BLUEPRINT pairing is blocked by the hg19–hg38 mismatch. Coordinate conversion can sometimes recover such comparisons, but it introduces additional preprocessing, may not map all intervals unambiguously, and continues to define similarity in terms of a selected coordinate system. An alternative is to extract the reference sequence underlying each interval and compare the resulting sequence collections directly, allowing related genomic annotations to be evaluated across references without requiring their coordinates to coincide.

![Figure 1](figures/roadmap_blueprint_top5_marker_pairwise_comparisons.png)

**Figure 1. Reference-genome fragmentation prevents direct comparison across public histone-mark collections.** Numbers of possible pairwise comparisons among processed histone-mark ChIP-seq peak BED files are shown for the five most represented marks shared by Roadmap Epigenomics and BLUEPRINT. Labels above the bars give pair counts; labels above the blocked cross-reference bars also give the percentage of all possible pairs for that mark. Within-resource comparisons can be performed directly because Roadmap files share hg19 coordinates and BLUEPRINT files share hg38 coordinates. In contrast, every Roadmap–BLUEPRINT pairing is blocked for direct coordinate-overlap analysis by the hg19–hg38 reference mismatch. File counts are deduplicated at the download-URL level; included outputs are BED-, narrowPeak-, broadPeak-, or gappedPeak-like peak files.

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

The Jaccard statistic, from set theory, [ETC] is one commonly used measure of similarity, offered by a number of tools. In particular, BEDtools [REF], a ubiquitous tool for genomic arithmetic, offers a `bedtools jaccard` command that calculates the exact value of intersection over union for two BED files. \program{hammock}, which is built to compare *lists* of interval files, takes advantage of two aspects of a point-wise sketching representation approach: batch reusability of a sketch once built and point-based reproducibly random subsampling.
As can be seen in Figure 3, the majority of \program{hammock}'s time is spent summarizing each interval file into a sketch, but the pairwise calculations using and reusing the resulting bit-vectors are trivially fast. BEDtools, by contrast, must reprocess each interval file for each new comparison. As one would expect, this means the \program{hammock}'s performance advantage scales with the number of BED files included in the analysis.

\program{hammock}'s \texttt{interval} mode treats each point within an interval as an object to be hashed and added to a single sketch representation of the entire set of intervals in the file. Given the usual size of genomics intervals, this creates an opportunity for further subsampling prior to sketching.  \program{hammock}'s novel `subB` implementation uses mixed-stride methodology to sample genomic positions using deterministic chromosome-specific strides and offsets, avoiding the per-position hashing cost of other reproducibly random methods such as hash-threshold subsampling. As seen in Figure 3, this substantially reduces sketch-construction time while preserving reproducibility and producing little change relative to the unsubsampled hammock similarities in the evaluated datasets. The reduction is sublinear in the sampling fraction and corpus-dependent: a 10-fold subsample buys 2.12\times{} on a synthetic corpus at N=512 (p=18) but only 1.87\times{} on Maurano at p=18 — and the synthetic-corpus benefit itself shrinks with catalog scale, from a peak of 2.62\times{} at N=16 down to 1.57\times{} at N=2048, as the pairwise-comparison phase (which `subB` does not touch) comes to dominate wall time (Supplementary Figure S10).




This section should distinguish three sources of the interval-mode speedup:

1. **Sketch reuse:** each BED file is read once and converted to a fixed-size summary that can be reused across all pairwise comparisons.
2. **Mixed-stride subsampling:** hammock's novel `subB` implementation samples genomic positions using deterministic chromosome-specific strides and offsets, avoiding the per-position hashing cost of hash-threshold subsampling. This substantially reduces sketch-construction time while preserving reproducibility and producing little change relative to the unsubsampled similarities. The reduction is sublinear in the sampling fraction and corpus-dependent — a 10-fold subsample buys only 1.87x on Maurano at p=18 but 2.12x on a synthetic corpus at N=512, and even the synthetic-corpus benefit shrinks with scale (2.62x at N=16 down to 1.57x at N=2048, Supplementary Figure S10) as the subB-invariant pairwise-comparison phase comes to dominate wall time — so it must not be described as proportional, or as a fixed multiplier independent of catalog size.
3. **Workflow parallelization:** The saturation is a result, not an obstacle. It's a genuine scalability argument for the evaluated workflows: the BEDTools process-per-pair workflow launches N² comparisons and re-reads each file N times, parallelizes poorly, and saturates early, whereas hammock launches one process and reads each file once. That's a real, defensible advantage. It just has to be attributed to batch-mode absence rather than to sketching, and reported as BEDTools' achieved efficiency rather than silently labeled "t=16."

![Figure 3](figures/pairwise_scaling.png)

**Figure 3. Hammock expands feasible all-pairs comparison as interval collections grow.** (A) **Sketch reuse increases the advantage as collections grow.** Wall time for hammock and BEDTools across synthetic collections containing 10,000 intervals per BED file, using HyperLogLog precision p=18 (the hammock CLI default), 16 threads, and three runs per configuration (twenty at N ≤ 32 — see generation note); N files per side gives N^{2} ordered pairs, 262,144 at N=512. Timings are of the standalone \texttt{hammock-cpp} binary, which the wheel installs but does not place on \texttt{PATH}. "BEDTools at 16 threads" is 16 concurrent single-threaded \texttt{bedtools jaccard} processes under GNU \texttt{parallel}, since \texttt{bedtools jaccard} has no batch or threaded mode; each BEDTools row records the parallelism it actually achieved, which on this corpus is 1.6 of 16 cores (Supplementary Fig S6). Curves show total hammock wall time on `jaccard_similarity_ie` — **the column we recommend reading for any comparison of magnitude**, and the only one this figure now reports — plus its sketch-comparison-phase component in isolation ("hammock sketch comparison"). BEDTools is faster than hammock below N≈64 (e.g. 0.29× at N=8) because hammock pays a roughly fixed per-invocation cost regardless of collection size while BEDTools' pairwise work is still cheap at this scale; the two total-wall curves cross between N=32 and N=64 as sketch reuse begins to dominate. At N=512, hammock is 8.4× faster than BEDTools. (B) **Subsampling further reduces runtime.** Wall time on 20 Maurano fetal-tissue DNase hypersensitivity BED files using interval mode with p=18 and eight threads, at the `jaccard_similarity_ie` arm (the full metrics block, costlier than a legacy register-equality-similarity-only run — see Fig S9). Bars compare BEDTools with hammock at `subB=1.0`, `0.1`, and `0.01`; labels report wall time, speedup relative to BEDTools, and MAE of `jaccard_similarity_ie` relative to exact BEDTools truth over the 190 unique off-diagonal pairs (self-comparisons excluded).

Together, the synthetic and real-data benchmarks show that hammock improves the feasibility of exhaustive pairwise analysis through both reusable sketches and efficient optional subsampling during sketch construction.

**Figure 3 generation.** Generated by `paper/pairwise_scaling/plot_pairwise_scaling.R`. Panel A: `docs/data/cpp_vs_bedtools_t16_p18.csv` (job 29671317). Panel B: `docs/data/maurano_bedtools.csv` and `docs/data/maurano_subB_ie_summary.csv`. Full provenance and re-measurement history: `docs/bedtools-baseline-retraction.md`.

#### 2.2.2 Interval sketches recover exact-overlap similarity structure and values

![Figure 4](figures/interval_accuracy.png)

**Figure 4. Inclusion–exclusion estimates closely reproduce exact interval Jaccard, with a tunable precision–runtime tradeoff.** (A) Pairwise comparisons were performed on 20 Maurano fetal-tissue DNase hypersensitivity BED files, yielding 190 unique off-diagonal pairs; self-comparisons were excluded. At the hammock CLI default, HyperLogLog precision p=18, the inclusion–exclusion Jaccard estimate closely tracks exact BEDTools base-pair Jaccard over the observed off-diagonal range (MAE=0.00115, Pearson r=0.99991, Spearman ρ=0.99966, Kendall τ=0.9883). (B) Accuracy and wall time were evaluated across p=12–24 on the same corpus using 16 threads and `subB=1.0`. The x-axis is MAE against exact BEDTools over the same 190 pairs on a base-10 logarithmic scale; lower values are more accurate and successive labeled ticks differ tenfold. The y-axis is median hammock wall time divided by median BEDTools wall time on a base-2 logarithmic scale; y=1 denotes equal runtime and each labeled tick doubles. Vertical ranges show the observed range among three hammock runs. Increasing precision reduced MAE from 0.00880 at p=12 to 0.000154 at p=24, while relative hammock wall time increased from 1.34× to 13.8×. Every point lies above y=1 because this 20-file corpus is below the approximately 64-file crossover at which reusable sketches amortize construction cost in the synthetic scaling benchmark. The magenta circle marks p=18, the CLI default.

**Figure 4 generation.** Both panels are generated by `paper/interval_accuracy/plot_interval_accuracy.R` with reference precision 18 and written to `paper/figures/interval_accuracy.png`. Panel A uses `docs/data/maurano_bedtools_ref.tsv` and `docs/data/hammock_hll_p18_jaccB_full.csv`. Panel B uses `docs/data/sweep_precision_maurano_p18_t16.csv`; accuracy is the MAE of `jaccard_similarity_ie`, and relative wall time is median hammock `wall_time` divided by median BEDTools wall time at each precision.

#### 2.3.1 Biological identity is preserved across references in sequence space

![Figure 5](figures/cross_reference_identity.png)

**Figure 5. Sequence sketches group samples by tissue across genome references.** Three ENCODE H3K27ac ChIP-seq samples—heart ENCSR175ABH, liver ENCSR864OOO, and lung ENCSR954JMZ—were independently aligned and peak-called against GRCh37, GRCh38, and CHM13, yielding nine peak sets. (A) Each broad-peak set was converted to sequence against its native reference and sketched in sequence mode at k=20, w=20, and HyperLogLog precision p=24. Average-linkage agglomerative hierarchical clustering on one minus the inclusion–exclusion Jaccard estimate groups the nine peak sets by tissue rather than by reference genome. (B) Median inclusion–exclusion Jaccard estimates for same-tissue comparisons across references and different-tissue comparisons are shown for selected cells from the expanded parameter sweep. Each vertical segment spans the two medians. The largest observed separation, Δ=0.563, occurred at k=20, w=20; its segment is highlighted in magenta.

**Figure 5 generation.** Panel A is generated by `paper/cross_reference_identity/plot_cross_reference_identity.R` from `docs/data/exp_a_broad_k20_w20_full.csv` and `docs/data/exp_a_metadata.tsv`. Panel B is generated by `paper/cross_reference_identity/plot_cross_reference_parameter_plateau.R` from the 63-cell-per-peak-type table `docs/data/exp_a_estimator_delta_expanded.csv`; the displayed broad-peak subset retains k=5–40 and selected tested window sizes. Both panels use `jaccard_similarity_ie`. `paper/cross_reference_identity/assemble_cross_reference_figure.R` stacks the title-free panels and writes `paper/figures/cross_reference_identity.png`.

#### 2.3.2 Sequence sketches recover tissue organization

- **Quantitative clustering result for the Results prose:** cutting the average-linkage dendrogram into the ten annotated tissue classes gives adjusted Rand index (ARI) = 0.9102 and normalized mutual information (NMI) = 0.9610. These values quantify a separate 10-cluster evaluation and should not be presented as properties of the visual tissue annotations in Figure 6.
- **Figure interpretation:** Figure 6A does not display a discrete cluster cut; its labels are colored and boxed by the known tissue annotations. Figure 6B shows the separate ten-cluster evaluation as lime-green outlines on the corresponding distance matrix.

![Figure 6](figures/sequence_tissue_clustering.png)

![Figure 6B](figures/sequence_tissue_distance_heatmap.png)

**Figure 6. Sequence sketches recover fetal-tissue organization in the Maurano DNase hypersensitivity collection.** Twenty fetal-tissue DNase-seq interval sets spanning ten annotated tissue labels were converted to interval-derived sequence and compared using the inclusion–exclusion Jaccard estimate, `jaccard_similarity_ie`, at k=10, w=30, and HyperLogLog precision p=18, the hammock CLI default. **(A)** Average-linkage agglomerative hierarchical clustering was performed on 1−J. Leaf labels show dataset accessions and are colored and boxed by annotated tissue; these annotations do not represent a discrete cluster cut. **(B)** The complete pairwise 1−J distance matrix follows the leaf order of the average-linkage dendrogram in A for direct comparison; this display ordering does not affect distances or cluster assignments. Accession-label colors and marginal bars denote annotated tissues, with arm, back, and leg grouped under muscle. Lime-green diagonal outlines show the result of cutting the average-linkage agglomerative hierarchy into ten clusters. They expose the departures from the annotated tissue partition: the two brain samples occupy separate singleton clusters, whereas the back sample and two leg samples occupy one cluster. The brain samples remain adjacent in the dendrogram leaf order, but their pairwise distance (1−J=0.136) exceeds the tightest within-tissue distances for heart (~0.069) and lung (~0.073). The outlined ten-cluster cut separates the brain pair while retaining those tighter groups, yielding ARI=0.910 rather than 1.0. The ten-cluster partition, ARI, and NMI are unchanged across every tested precision from p=12 through p=24.

**Figure 6 generation.** Panel A is generated by `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` and written to `paper/figures/sequence_tissue_clustering.png`; Panel B is generated by `paper/sequence_tissue_clustering/plot_distance_heatmap.R` and written to `paper/figures/sequence_tissue_distance_heatmap.png`. Both scripts read `experiments/maurano_dhs_validation/results/raw_d/hammock_mnmzr_p18_jaccD_k10_w30.csv`, use `jaccard_similarity_ie`, and obtain tissue annotations from `experiments/maurano_dhs_validation/data/maurano_filenames_key.tsv`. The archived sweep predates the native inclusion–exclusion output column, so the scripts derive it from `containment_AB` and `containment_BA`.

### 2.3.3 Parameterization separates numerical agreement from biological resolution

- **Figure 6 parameter choice for the Results prose:** Under inclusion–exclusion, k=10, w=30 attains the maximum ARI at the displayed default precision p=18 and at every tested precision p=12–24; at p=12, k=10, w=20 ties that maximum. The k=10, w=30 ten-cluster partition persists across the tested precision range; the parameter sweep, rather than Figure 6 itself, is the evidence for that choice.

![Figure 7](figures/parameter_objective_tradeoff_estimators.png)

**Figure 7. Sequence parameters separate numerical agreement from biological recovery.** Sequence-mode parameter configurations were evaluated on 20 Maurano fetal-tissue DNase hypersensitivity interval sets at HyperLogLog precision p=24 using the inclusion–exclusion Jaccard estimate (`jaccard_similarity_ie`). Each point represents one combination of k-mer size (k) and minimizer-window size (w); the x-axis gives Pearson correlation with exact BEDTools Jaccard across pairwise comparisons, and the y-axis gives adjusted Rand index (ARI) for recovery of the ten annotated tissue classes after clustering. Point color encodes k and point size encodes w. The configuration with greatest numerical agreement, k=20 and w=20 (r=0.99997; ARI=0.693), is emphasized by a green circle. The k=10, w=30 biological optimum used in Figure 6 (ARI=0.910; r=0.946) is marked by a yellow circle. Thus the parameters that best reproduce exact pairwise Jaccard values differ from those that best recover the annotated tissue organization. The exact BEDTools reference was independently regenerated with BEDTools 2.30.0 over all 400 ordered pairs and matched the checked-in reference byte-for-byte, with 0 of 190 unique pairs asymmetric.

**Figure 7 generation.** The figure is generated by `paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R` and written to `paper/figures/parameter_objective_tradeoff_estimators.png`. The script reads `docs/data/mode_d_summary.csv`, supplementing the inclusion–exclusion rows from `experiments/maurano_dhs_validation/results/mode_d_summary.csv` when they are absent from the staged paper copy, and filters to `jaccard_similarity_ie`, p=24, BEDTools reference, k≥8, and non-missing Pearson/ARI values. It identifies the greatest-Pearson configuration and verifies the k=10, w=30 ARI optimum.

### Supplementary Figure S1 — interval accuracy at higher precision

![Figure S1](figures/interval_accuracy_p21.png)

**Figure S1. Inclusion–exclusion estimates closely reproduce exact interval Jaccard at higher HyperLogLog precision.** (A) Pairwise comparisons on 20 Maurano et al. fetal-tissue DNase hypersensitivity BED files yielded 190 unique off-diagonal pairs after self-comparisons and reciprocal duplicates were excluded. At p=21, `jaccard_similarity_ie` closely tracked exact BEDTools base-pair Jaccard (MAE=0.000430, Pearson r=0.99999, Spearman ρ=0.99987, and Kendall τ=0.9947). (B) The precision–runtime analysis from Figure 4 is repeated for context. The x-axis is MAE on a base-10 logarithmic scale, where lower is more accurate and successive labeled ticks differ tenfold. The y-axis is hammock/BEDTools relative wall time on a base-2 logarithmic scale, where y=1 denotes equal runtime and each labeled tick doubles; vertical ranges show the three-run hammock range. The magenta circle marks p=18, the CLI default.

**Figure S1 generation.** `paper/interval_accuracy/plot_interval_accuracy.R` with its optional reference-precision argument set to 21, using `docs/data/maurano_bedtools_ref.tsv`, `docs/data/hammock_hll_p21_jaccB_full.csv`, and `docs/data/sweep_precision_maurano_p18_t16.csv`; output is `paper/figures/interval_accuracy_p21.png`.

### Supplementary Figure S5a — tissue organization at the numerical optimum

![Figure S5a, panel A](figures/sequence_tissue_clustering_k20_w20.png)

![Figure S5a, panel B](figures/sequence_tissue_distance_heatmap_k20_w20.png)

**Figure S5a. Numerical agreement with exact interval Jaccard does not maximize recovery of fetal-tissue organization.** The same 20 fetal-tissue DNase-seq interval sets shown in Figure 6 were compared using the inclusion–exclusion Jaccard estimate, `jaccard_similarity_ie`, at k=20, w=20, and HyperLogLog precision p=24—the parameter combination with the greatest Pearson correlation with exact BEDTools Jaccard within the p=24 comparison shown in Figure 7 (r=0.99997). **(A)** Average-linkage agglomerative hierarchical clustering was performed on 1−J. Leaf labels show dataset accessions and are colored and boxed by annotated tissue; these annotations do not represent a discrete cluster cut. **(B)** The complete pairwise 1−J distance matrix follows the dendrogram leaf order in A; this display ordering does not affect distances or cluster assignments. Accession-label colors and marginal bars denote annotated tissues, with arm, back, and leg grouped under muscle. Lime-green diagonal outlines show the result of cutting the average-linkage agglomerative hierarchy into ten clusters. Relative to the annotated tissue partition, the cut separates the two brain samples, separates one of the five small-intestine samples, and combines the arm, back, and two leg samples into one cluster, yielding ARI=0.693 and NMI=0.905. In contrast, the biologically selected k=10, w=30 setting in Figure 6 at the default p=18 yields ARI=0.910 and NMI=0.961, demonstrating that numerical agreement with exact pairwise Jaccard and recovery of tissue organization select different sequence-sketch parameters.

**Figure S5a generation.** Panel A is generated by `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` using `experiments/maurano_dhs_validation/results/raw_d/hammock_mnmzr_p24_jaccD_k20_w20.csv`, the Maurano tissue key, `paper/figures/sequence_tissue_clustering_k20_w20.png`, and `jaccard_similarity_ie` as its four positional arguments. Panel B is generated by `paper/sequence_tissue_clustering/plot_distance_heatmap.R` with the same input CSV, tissue key, `paper/figures/sequence_tissue_distance_heatmap_k20_w20.png`, and similarity-column argument.

## III. Discussion

### 3.1 Scaling comparison within references

Discuss both layers of the computational contribution: reusable sketches reduce the cost of repeated all-pairs comparisons, while mixed-stride makes optional base-pair subsampling computationally effective by avoiding a per-position hashing gate. Mixed-stride should be described as a novel deterministic subsampling strategy with chromosome-specific phase offsets, not merely as a command-line tuning option.

### 3.2 Comparing interval-derived sequence across references

### 3.3 Coordinate similarity and sequence similarity are complementary

### 3.4 Practical recommendations

### 3.5 Limitations and future work

Every wall-time result reported here, including Figure 3 and its supplements, times the standalone \texttt{hammock-cpp} binary rather than the \program{hammock} Python command-line interface most users invoke directly (the installed wheel provides \texttt{hammock-cpp} but does not place it on \texttt{PATH}). We measured the gap between the two front-ends directly, at the same precision and thread count used throughout this study (p=18, 16 threads, synthetic catalogs of N=2 to N=2048 files per side, \texttt{--exclusive} SLURM allocation). The gap is not a fixed overhead: the two front-ends differ in how they dispatch per-file sketch construction — \texttt{hammock-cpp} sketches files serially, one dedicated thread team at a time, while the Python CLI dispatches multiple files concurrently through a thread pool. This makes the CLI faster than \texttt{hammock-cpp} across a middle range of catalog sizes (N≈8–1024 files per side, peaking near 28\% faster at N=32–128), converging back toward parity by N=1024 and reversing by N=2048, where \texttt{hammock-cpp} is roughly 13\% faster — the catalog size closest to the scaling regime Figure 3 addresses. Because \texttt{hammock-cpp} is the faster of the two tools specifically in the N≥1024 regime this study's scaling claims target, the wall times reported here are not shown to understate what the CLI most users run would achieve at that scale. At smaller-to-moderate catalog sizes, users invoking the CLI directly may see somewhat faster wall times than the \texttt{hammock-cpp} results presented here would suggest. A candidate explanation for the reversal — that the pairwise-comparison phase, identical compiled code in both front-ends and quadratic in N, comes to dominate wall time at large N and dilutes a dispatch advantage confined to the sketching phase — is consistent with the data but has not been confirmed by timing the two phases separately. Combining both front-ends' advantages, by parallelizing \texttt{hammock-cpp}'s per-file sketch dispatch to match the CLI's, is a natural direction for future work.

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
| `reg_eq_similarity`    | Legacy register-equality similarity: fraction of active HyperLogLog registers with equal values; not set Jaccard |
| `containment_AB`        | Estimated fraction of dataset A contained in dataset B                         |
| `containment_BA`        | Estimated fraction of dataset B contained in dataset A                         |
| `cosketch_geom`         | Geometric mean of the two directional containment estimates                    |
| `cosketch_arith`        | Arithmetic mean of the two directional containment estimates                   |
| `cosketch_max`          | Maximum of the two directional containment estimates                           |

Analyses that estimate set Jaccard magnitudes or compare hammock with exact BEDTools Jaccard use the inclusion–exclusion estimate, `jaccard_similarity_ie`; Figure 6 and its companion heatmap also use this estimator. Sequence-mode biological analyses that use the register-equality score, `reg_eq_similarity`, identify it explicitly; despite its legacy column name, this score is not set Jaccard and should not be interpreted on the same numerical scale. The legacy name is retained for output-schema and downstream-workflow compatibility. Sketch construction and interval subsampling are reproducible for fixed parameter values and seeds.

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
- **Register equality (`reg_eq_similarity`; legacy column name retained for compatibility).** Define as the fraction of active registers whose values agree.
  - *"This is not set Jaccard: two registers agree whenever the largest ρ observed in that bucket happens to coincide, which occurs at some rate even for disjoint inputs, so the statistic carries a chance-agreement floor c."*
  - c is set by the sketch load factor λ = n/m **and** by the cardinality ratio |A|/|B| — not by precision as such. Report 0.1699 at equal cardinality as λ→∞, 0.058 at ratio 10, and 0.045 in the ~5-Mbp synthetic benchmark at p=24, where m exceeded n.
  - *"Because c differs between pairs of differing size, register-equality is order-preserving within a cardinality ratio but not across ratios."* Ties to the 2.45% Maurano inversions in §2.3.
- **Inclusion–exclusion (**`jaccard_similarity_ie`**).** Define |A ∩ B| = |A| + |B| − |A ∪ B| with each term an Ertl (2017) estimate, union by register-wise maximum, intersection clamped to ≥ 0; J = |A∩B|/|A∪B|, equivalently 1/(1/C_AB + 1/C_BA − 1).
  - *"An exact 0.0 in this column means the intersection estimate hit the non-negativity clamp or the inputs are genuinely disjoint; it never means a measured zero."*
  - Relative error ≈ 0.6/(J√m), hence uninformative below J ≈ a few/√m — state as the column's domain of validity.
- **Cost.**
  - *"Register-equality is a ratio of two register counts and costs a single pass over the register arrays; inclusion–exclusion additionally requires a union sketch and cardinality estimates per pair, so its cost scales as O(N²·2^p) against O(N·M) for sketch construction."*
  - Measured overhead, 16 threads, p = 18 (`docs/data/cpp_vs_bedtools_t16_p18.csv`, job 29671317): the inclusion–exclusion comparison phase costs **1.64–1.94× as much as** the register-equality comparison phase for N ≥ 64. Its share of inclusion–exclusion total wall rises from **1.6% at N = 32 to 22.1% at N = 512**, while total inclusion–exclusion wall is 0.7%, 2.6%, 5.3%, and 10.0% higher at N = 32, 128, 256, and 512, respectively (78.08 s against 70.97 s at N = 512). Peak RSS is 278 MB against 266 MB.
  - *"The overhead is modest but no longer negligible at the default precision: the inclusion–exclusion comparison phase grows from 1.6% of wall at N = 32 to 22.1% at N = 512 because the phase is Θ(N²) against Θ(N) sketching."* Its growing share is visible in these data, but it does not yet dominate total wall at the largest measured N.
  - Figure 3 plots `jaccard_similarity_ie` only, so the speed claim and the reporting recommendation refer to literally the same run. The register-equality comparison is Supplementary Figure S9.
- **Why not the joint maximum-likelihood estimator.** *"Ertl's joint MLE is the lower-variance choice for HyperLogLog Jaccard, but it requires the joint register histogram of the pair, which hammock's sketch interface does not expose; the inclusion–exclusion form is recoverable from the containment estimates already computed and needs no additional interface."* Pre-empts the obvious reviewer question.
- **Which column to report.** *"We report* `jaccard_similarity_ie` *for any comparison of magnitude, and for any comparison spanning different set sizes;* `reg_eq_similarity` *is retained as the low-cost default."*
  - One clause on the exception (register-equality ranks better below J ≈ 0.05 at p ≤ 20 among comparably-sized pairs, and one step of precision removes the advantage); the quantification is deferred rather than carried in Methods.
  - **Figure 6 estimator note:** Figure 6 uses `jaccard_similarity_ie`, derived from the archived sweep's containment columns. At Figure 6's displayed configuration (k=10, w=30, p=18), inclusion–exclusion and register equality induce the *identical* 10-cluster partition, with ARI and NMI agreeing to sixteen digits. Under inclusion–exclusion, the same partition and scores occur at every tested precision p=12–24, although the complete topology and leaf order match p=24 only at p≥22. The estimator equivalence is **local, not global**: across the 235-cell sweep the two columns give different ARI at 48 cells (max |ΔARI| 0.301). Generated by `paper/sequence_tissue_clustering/estimator_agreement_stats.csv` and `experiments/maurano_dhs_validation/results/estimator_ari_by_config.csv`.
- **Complexity.** Retain the original stub: cost of mixed-stride sketch construction as a function of covered length and stride, alongside full interval ingestion and pairwise sketch comparison.

### 4.5 Datasets

### 4.6 Performance benchmarking

Include the direct comparison of no subsampling and mixed-stride `subB = 0.1` and `0.01` used in Figure 3B. The main benchmark establishes the practical speed–approximation trade-off; exhaustive comparison against the hash-threshold and single-hash alternatives is deferred.

### 4.7 Interval-mode accuracy evaluation

### 4.8 Sequence-mode biological evaluation

## Data and code availability

## Author contributions

## Acknowledgments

## References
- **TODO — sequence-mode parameter heatmaps:** include the full k × w parameter-response heatmaps for numerical agreement with BEDTools and biological clustering performance (ARI/NMI), potentially stratified by HyperLogLog precision and/or estimator. These heatmaps should provide the complete response surface underlying the compact objective-space summary in Figure 7 without burdening the main-text figure.
