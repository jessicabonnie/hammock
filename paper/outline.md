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

### 2.2 Interval sketching expands feasible all-pairs BED comparison

[NOTE: find a place for this] Interval comparisons are typically conducted using files in the BED format, a simple tab-delimited representation of genomic intervals that has been in wide use since the early 2000s\cite{kent_hgbrowser}.

The Jaccard statistic, from set theory, [ETC] is one commonly used measure of similarity, offered by a number of tools. In particular, BEDtools [REF], a ubiquitous tool for genomic arithmetic, offers a `bedtools jaccard` command that calculates the exact value of intersection over union for two BED files. \program{hammock}, which is built to compare *lists* of interval files, takes advantage of two aspects of a point-wise sketching representation approach: batch reusability of a sketch once built and point-based reproducibly random subsampling.
As can be seen in Figure 3, the majority of \program{hammock}'s time is spent summarizing each interval file into a sketch, but the pairwise calculations using and reusing the resulting bit-vectors are trivially fast. BEDtools, by contrast, must reprocess each interval file for each new comparison. As one would expect, this means the \program{hammock}'s performance advantage scales with the number of BED files included in the analysis.

\program{hammock}'s \texttt{interval} mode treats each point within an interval as an object to be hashed and added to a single sketch representation of the entire set of intervals in the file. Given the usual size of genomics intervals, this creates an opportunity for further subsampling prior to sketching.  \program{hammock}'s novel `subB` implementation uses mixed-stride methodology to sample genomic positions using deterministic chromosome-specific strides and offsets, avoiding the per-position hashing cost of other reproducibly random methods such as hash-threshold subsampling. As seen in Figure 3, this reduces sketch-construction time in proportion to the requested sampling fraction while preserving reproducibility and producing little change relative to the unsubsampled hammock similarities in the evaluated datasets.

Performance improvements are not uniform.


This section should distinguish two sources of the interval-mode speedup:

1. **Sketch reuse:** each BED file is read once and converted to a fixed-size summary that can be reused across all pairwise comparisons.
2. **Mixed-stride subsampling:** hammock's novel `subB` implementation samples genomic positions using deterministic chromosome-specific strides and offsets, avoiding the per-position hashing cost of hash-threshold subsampling. This reduces sketch-construction time in proportion to the requested sampling fraction while preserving reproducibility and producing little change relative to the unsubsampled hammock similarities in the evaluated datasets.
3. **Workflow parallelization:** The saturation is a result, not an obstacle. It's a genuine scalability argument for sketching: an exact per-pair tool must launch N² processes and re-read each file N times, and that workflow does not parallelize — hammock launches one process and reads each file once. That's a real, defensible advantage. It just has to be attributed to batch-mode absence rather than to sketching, and reported as bedtools' achieved efficiency rather than silently labeled "t=16."
- hammock's is intra-process: OpenMP threads in one address space, each file read once, sketches shared. This is completely unaffected by everything I found today. It scales, and we can show a real thread-scaling curve.
- bedtools' is workflow-level: N² independent process launches under a GNU parallel wrapper. On this cluster that saturates at ~1.2–2.9× no matter how many cores you supply, because process creation caps near 123 exec/s.
- Concretely, three ways to put parallelization in the figure, all supported by data we already have or can cheaply get:

1. A threads panel — wall or throughput vs thread count at fixed N, both tools, efficiency annotated. hammock rises; bedtools plateaus. This makes the point vividly and is self-documenting. "One thing worth flagging for when you write the caption: this figure makes the mechanism argument, not a "we're faster" argument. The honest reading is that hammock converts cores into throughput and the bedtools pairwise workflow stops doing so almost immediately — attributable to batch-mode absence, not to sketching."
2. CPU time beside wall — already recorded in every CSV. A reader can see cores actually used, which is precisely the subtraction a skeptical referee makes.
3. Anchor speedups on bedtools t=1 as well, so the headline doesn't depend on how well the wrapper happened to parallelize on a given node.

Mixed-stride should be presented as a methodological contribution within the interval-mode scaling result, not as an unrelated implementation optimization. The main text should briefly contrast it with hash-threshold subsampling; the full comparison among mixed-stride, hash-threshold, and single-hash strategies can remain supplementary.

- **Synthetic scaling result (Figure 3A):** each benchmark configuration comprises two disjoint collections of N files compared as a full cross product, so the number of comparisons grows as N^{2}, reaching 262,144 pairs at N=512. Because hammock constructs one reusable sketch per file while the BEDTools workflow repeatedly processes the underlying interval files, the performance advantage of sketch reuse increases with collection size.
- **Sketch-comparison output cost (Figure 3A):** the lower curves separate the three-column sketch-comparison output used for timing from the analysis-oriented `+IE` output, which additionally emits `jaccard_similarity_ie`, containment, and co-sketch columns. The additional columns cost a factor of 3.4–3.8 within the comparison phase at N\ge64, but this phase is small relative to sketch construction: at N=512, `+IE` increases total wall time by only 1.45% (54.15 s versus 53.38 s). The two total-time curves are therefore indistinguishable, and only the three-column total is drawn.
- **Real-data timing result (Figure 3B):** on the 20 Maurano fetal-tissue DNase hypersensitivity BED files, hammock is faster than the parallelized BEDTools workflow even without subsampling.
- **Mixed-stride result (Figure 3B):** `subB=0.1` and `subB=0.01` further reduce runtime while producing little change relative to hammock's own unsubsampled similarity estimates. Mixed-stride deterministically advances through genomic positions at chromosome-specific strides and offsets rather than applying a hash-based inclusion test at every covered position.
- **Definition of the plotted approximation change (Figure 3B):** \Delta J is defined per pair as the difference between a subsampled run's `jaccard_similarity` and the mean `jaccard_similarity` from hammock's own unsubsampled (`subB=1.0`) runs on that same pair. Mean |\Delta J| is the mean absolute value over the 190 pairs. (An earlier draft reported this as "950 pair-by-replicate comparisons (190 pairs \times 5 replicates)". The five replicates are byte-identical --- hammock is deterministic given the seed, and all 3420 archived (method, subB, pair) cells contain exactly one distinct Jaccard value across the five --- so the replicates carry no additional information and the effective n is 190. The reported mean is unchanged; only the stated n was wrong.)

![Figure 3](figures/pairwise_scaling.png)

**Figure 3. Hammock expands feasible all-pairs comparison as interval collections grow.** (A) Wall time for hammock and BEDTools across synthetic collections containing 10,000 intervals per BED file, using HyperLogLog precision p=14, 16 threads, and three runs per configuration. Curves show total hammock wall time as well as sketch-construction and sketch-comparison components; the `+IE` comparison curve includes `jaccard_similarity_ie`, containment, and co-sketch output. (B) Wall time on 20 Maurano fetal-tissue DNase hypersensitivity BED files (190 unique pairs) using interval mode with p=18 and eight threads. Bars compare BEDTools with hammock at `subB=1.0`, `0.1`, and `0.01`; labels report wall time, speedup relative to BEDTools, and mean |\Delta J| relative to hammock's unsubsampled similarity estimates.

Together, the synthetic and real-data benchmarks show that hammock improves the feasibility of exhaustive pairwise analysis through both reusable sketches and efficient optional subsampling during sketch construction.

**Figure 3 generation.** The combined figure is generated by `paper/pairwise_scaling/plot_pairwise_scaling.R` and written to `paper/figures/pairwise_scaling.png`. Panel A uses `docs/data/cpp_vs_bedtools_t16_20260804_172242.csv` (job 29552415); Panel B uses `docs/data/maurano_subB_summary.csv` and `docs/data/maurano_bedtools.csv`. The script rebuilds both panels with harmonized typography, panel labels, tool colors, wall-time units, and the pair count N^{2}.

### 2.3 Interval sketches recover exact-overlap similarity structure and values

![Figure 4](figures/interval_accuracy.png)

**Figure 4. Inclusion–exclusion estimates closely reproduce exact BEDTools Jaccard.** Pairwise comparisons were performed on 20 Maurano fetal-tissue DNase hypersensitivity BED files, yielding 190 unique off-diagonal pairs; self-comparisons were excluded. Hammock's `jaccard_similarity_ie` uses HyperLogLog cardinality estimates with inclusion–exclusion to approximate the same base-pair Jaccard statistic computed exactly by BEDTools. At HyperLogLog precision p=21, the estimated values closely track exact BEDTools Jaccard over the observed off-diagonal range (MAE =4.3\times10^{-4}, Pearson r=0.99999, Kendall \tau=0.9947). Across precisions p=18, p=21, and p=23, absolute deviation from BEDTools decreases as precision increases, from approximately 1\times10^{-3} to approximately 2\times10^{-4}. The figure therefore evaluates the numerical accuracy of the inclusion–exclusion estimate relative to exact interval Jaccard; it does not compare inclusion–exclusion with hammock's register-equality statistic.

**Figure 4 generation.** The figure is generated by `paper/interval_accuracy/plot_interval_accuracy.R` and written to `paper/figures/interval_accuracy.png`. Inputs are `docs/data/maurano_bedtools_ref.tsv` and `docs/data/hammock_hll_p{18,21,23}_jaccB.csv`. Those three CSVs were regenerated with hammock 0.5.0 on 2026-07-31 to obtain the `jaccard_similarity_ie` column; earlier copies carried the placeholder `containment` column and could not supply it. The regeneration reproduces the archived `jaccard_similarity` values **exactly** (0 of 400 pairs differ at any precision), so panel A's register-equality series is unchanged from the previously published figure.

### 2.4 Biological identity is preserved across references

![Figure 5](figures/cross_reference_identity.png)

**Figure 5. Sequence sketches group samples by tissue across genome references.** Three ENCODE H3K27ac ChIP-seq samples—heart ENCSR175ABH, liver ENCSR864OOO, and lung ENCSR954JMZ—were independently aligned and peak-called against GRCh37, GRCh38, and CHM13, yielding nine peak sets. Each peak set was converted to sequence against its native reference and sketched in sequence mode using broad peaks, k=10, w=10, HyperLogLog precision p=24, and the minimizer-only Jaccard metric. UPGMA clustering on 1-J groups the nine peak sets by tissue rather than by reference genome: heart, liver, and lung each form a distinct clade containing the GRCh37-, GRCh38-, and CHM13-derived representation of that tissue. In this proof-of-concept panel, reference choice therefore behaves as a within-tissue perturbation rather than as the dominant axis of separation, showing that sequence sketches can retain biological identity across genome references without requiring direct coordinate overlap.

**Figure 5 generation.** The figure is generated by `paper/cross_reference_identity/plot_cross_reference_identity.R` and written to `paper/figures/cross_reference_identity.png`. Inputs are `docs/data/exp_a_broad_k10_w10.csv` (staged from `experiments/ref-comparison/results/exp_a/broad/k10_w10/exp_a_mnmzr_p24_jaccD_k10_w10.csv`) and `docs/data/exp_a_metadata.tsv` (staged from the experiment configuration). The script constructs the nine-sample similarity matrix, converts similarity to 1-J, performs average-linkage hierarchical clustering, and renders the tissue-colored dendrogram.

### 2.5 Sequence sketches recover tissue organization

- **Quantitative clustering result for the Results prose:** cutting the average-linkage dendrogram into the ten annotated tissue classes gives adjusted Rand index (ARI) = 0.9102 and normalized mutual information (NMI) = 0.9610. These values quantify a separate 10-cluster evaluation and should not be presented as properties of the visual tissue annotations in Figure 6.
- **Figure interpretation:** Figure 6 itself does not display a discrete cluster cut. Its labels are colored by the known tissue annotations to show how biological identity is organized along the dendrogram.

![Figure 6](figures/sequence_tissue_clustering.png)

**Figure 6. Sequence sketches recover fetal-tissue organization in the Maurano DNase hypersensitivity collection.** Twenty fetal-tissue DNase-seq interval sets spanning ten annotated tissue labels were converted to interval-derived sequence and compared using the minimizer-only `jaccard_similarity` score at k=10, w=30, and HyperLogLog precision p=24. Average-linkage hierarchical clustering was performed on 1-J. Leaf labels show dataset accessions and are colored by annotated tissue. The tissue coloring is an annotation of the dendrogram and does not represent a discrete cluster cut.

**Figure 6 generation.** The figure is generated by `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` and written to `paper/figures/sequence_tissue_clustering.png`. By default, the script reads `experiments/maurano_dhs_validation/results/raw_d/hammock_mnmzr_p24_jaccD_k10_w30.csv`, uses the `jaccard_similarity` column, and obtains tissue annotations from `experiments/maurano_dhs_validation/data/maurano_filenames_key.tsv`. Optional command-line arguments can substitute another similarity CSV, tissue key, output path, or similarity column. Passing `jaccard_similarity_ie` as the fourth argument produces Supplementary Figure S5a; the column is derived from the containment columns the archived CSV already carries, since the sweep predates it. Two further renderings of the same configuration at different precisions, `paper/figures/sequence_tissue_clustering_p12.png` and `sequence_tissue_clustering_p24.png`, record that the clustering is unchanged across the p ≥ 12 plateau and are retained as supporting material rather than manuscript figures.

### 2.6 Parameterization separates numerical agreement from biological resolution

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
- **Cost, and why the two columns are not interchangeable operationally.**
  - *"Register-equality is a ratio of two register counts and costs a single pass over the register arrays; inclusion–exclusion additionally requires a union sketch and cardinality estimates per pair, so its cost scales as O(N²·2^p) against O(N·M) for sketch construction."*
  - Measured overhead, same binary, 16 threads, p = 14: not measurable at N = 20 (p = 18), +2.4% at N = 100, **+10.4% at N = 512** — the largest collection in Figure 3A.
  - *"The overhead is therefore negligible for small collections and material at the catalog scale the method targets."*
  - State explicitly which configuration Figure 3 benchmarks, so the speed claim and the reporting recommendation refer to the same thing.
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

**Cross-reference identity (Figure 5).** Redrawn on `jaccard_similarity_ie` (`paper/figures/cross_reference_identity_ie.png`), the nine peak sets still form three monophyletic tissue clades, each containing all three reference-genome representations. Ranking the twenty (k, w) cells by the separation statistic Δ = median(same-tissue cross-reference) − median(different-tissue) gives Spearman ρ = 0.9925 between the two columns on both broad and narrow peaks, with rank changes confined to the top five cells, which lie within 0.06 of one another. The k ≥ 15 and k ≤ 10 groups remain disjoint under both columns, with a slightly wider gap on `jaccard_similarity_ie` (broad +0.415 against +0.389).

Caveat worth carrying: `jaccard_similarity_ie` is censored at zero by the non-negativity clamp on the inclusion–exclusion intersection, so an exact 0.0 means "clamped or empty", never "measured zero". This does not arise at the configurations above but would in a low-similarity corpus.

## Supplementary Figures and Tables

- **Figure S5a — sequence-mode tissue dendrogram on `jaccard_similarity_ie`:** `paper/figures/sequence_tissue_clustering_k10_w30_p24_ie.png`, generated by `paper/sequence_tissue_clustering/plot_sequence_tissue_clustering.R` with `jaccard_similarity_ie` as the fourth positional argument. Companion to Figure 6; see S5.
- **Figure S5b — cross-reference dendrogram on `jaccard_similarity_ie`:** `paper/figures/cross_reference_identity_ie.png`, generated by `paper/cross_reference_identity/plot_cross_reference_identity.R` with the output path and `jaccard_similarity_ie` as its two positional arguments. Companion to Figure 5; see S5.
- **Table S5 — estimator agreement:** `paper/sequence_tissue_clustering/estimator_agreement_stats.csv` (per-column ARI/NMI and partition identity at Figure 6's cell), `experiments/maurano_dhs_validation/results/estimator_ari_by_config.csv` (all 235 cells), `experiments/maurano_dhs_validation/results/best_config_by_column.csv` (per-column optimum with tie sets), and `experiments/ref-comparison/results/exp_a_estimator_delta.csv` (cross-reference Δ per cell).
- **TODO — sequence-mode parameter heatmaps:** include the full k × w parameter-response heatmaps for numerical agreement with BEDTools and biological clustering performance (ARI/NMI), potentially stratified by HyperLogLog precision and/or estimator. These heatmaps should provide the complete response surface underlying the compact objective-space summary in Figure 7 without burdening the main-text figure.