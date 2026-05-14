# Experiment Design: Minimizer-Based Sketching on ENCODE Epigenomic Data

**Project:** Validation of minimizer-based sequence sketching across reference genomes and biological contexts  
**Status:** Workflow built; pending accession confirmation and parameter selection  
**Last Updated:** 2026-04 (Exp A scoping); 2026-05 (Exp B split out — see note below)

> **Note (2026-05-12):** This document still contains the original Experiment B
> (tissue-over-species clustering) sections for historical reference, but
> Exp B has been split into separate experiments under
> `hammock_claude/experiments/`:
> - `mus-homo/` — tissue-over-species in human + mouse (probably H3K4me3)
> - `primate-phylogeny/` — phylogeny recovery across mammalian liver (Villar 2015)
> - `man-monkey/` — tissue-over-species at primate distance, non-ENCODE data
>
> This directory (`ref-comparison`) is now **Experiment A only** (cross-reference
> robustness). All Exp B sections below are kept for reference but the Snakefile
> here does not run them.

---

## Overview

This document describes the experimental design for two complementary validation experiments using minimizer-based sequence sketching applied to ENCODE and Roadmap Epigenomics ChIP-seq data:

1. **Experiment A — Reference Robustness:** Demonstrate that minimizer sketches of the same biological sample aligned to GRCh37 vs. GRCh38 are more similar to each other than to sketches from a different sample aligned to either reference.
2. **Experiment B — Tissue-over-Species Clustering:** Reproduce the known result that epigenomic profiles cluster by tissue type rather than species when comparing human and mouse samples, using minimizer sketching as the similarity metric.

Both experiments sketch from aligned reads (BAM → FASTA). Experiment B additionally runs a parallel peak-called track (BAM → peak BED → FASTA) as the primary expected path to recovery of tissue-over-species clustering; the alignment-only track is retained as a secondary check. Experiment A uses alignment-only (peak calling is not meaningful for the cross-reference robustness question).

---

## Biological Motivation and Prior Work

### Tissue-over-Species Clustering (Experiment B)

The clustering of epigenomic samples by tissue type across species rather than by species is a well-established finding. The key references establishing this hierarchy are:

| Reference | Key Finding | Data Used |
|---|---|---|
| Yue et al. (2014) *Nature* [@Yue2014] | Mouse and human chromatin state landscapes cluster by tissue; RNA-seq and DNase-seq confirm tissue identity dominates species identity | Mouse ENCODE (GEO: see Supplementary) and human ENCODE/REMC |
| Barbosa-Morais et al. (2012) *Science* [@BarbosaMorais2012] | Alternative splicing patterns cluster by tissue across vertebrates more strongly than by species | RNA-seq across 9 species |
| Merkin et al. (2012) *Science* [@Merkin2012] | Conserved and species-specific splicing across mammals; tissue dominates over species in clustering | RNA-seq across mammals |
| Lin et al. (2014) *Nature* [@Lin2014] | Transcriptome-level comparison across 15 tissues in human and mouse: inter-tissue differences exceed inter-species differences for tissue-specific genes | RNA-seq; histone ChIP-seq (REMC + Mouse ENCODE) |
| Roadmap Epigenomics Consortium (2015) *Nature* [@Roadmap2015] | 111 human reference epigenomes cluster by tissue/lineage; H3K4me1, H3K27ac most informative for cell type identity | GEO: GSE16256; ENCODE accessions E001–E113 |
| ENCODE Phase III (2020) *Nature* [@ENCODE2020] | Cross-species comparison of cis-regulatory elements confirms tissue-specific clustering in both human and mouse | ENCODE portal: https://www.encodeproject.org |

The histone marks **H3K4me3** (active promoters), **H3K27ac** (active enhancers), and **H3K4me1** (poised/active enhancers) are the most informative for distinguishing tissue types and have been consistently used across all major studies. These are therefore the primary targets for Experiment B.

### Reference Genome Robustness (Experiment A)

No prior work has explicitly tested sketch-based similarity across reference versions, but the experiment is motivated by:

- GRCh37 (hg19) → GRCh38 (hg38) differences are predominantly in unplaced scaffolds, alternate loci, and refined centromere/telomere sequences; the core euchromatic sequence is >99.9% identical
- Minimizers derived from reads aligned to either reference capture the same underlying biological k-mer content, modulo multi-mapping edge cases
- This robustness property is important for users wishing to compare datasets across legacy (GRCh37) and current (GRCh38) alignments

---

## Experiment A: Cross-Reference Robustness

### Design

**Input:** Paired BAM files for the same biological sample, one aligned to GRCh37 and one to GRCh38.

**Samples:** Select 4–6 ENCODE ChIP-seq experiments for H3K27ac in human tissues where both GRCh37 and GRCh38 alignments are available via the ENCODE portal. Select 2–3 distinct tissue types to serve as "negative controls" (different-tissue pairs should have lower sketch similarity than same-tissue cross-reference pairs).

**Tissue targets (Experiment A):**

| Tissue | Roadmap ID | Mark |
|---|---|---|
| Heart left ventricle | E095 | H3K27ac |
| Liver | E066 | H3K27ac |
| Lung | E096 | H3K27ac |

> **TODO:** Confirm SRA run IDs for Roadmap E095/E066/E096 via `scripts/fetch_accessions.py` and update `sra_map` in `config/config.yaml`.

### Metric

For each BAM file, extract minimizer sketches using the project's sketching tool. Compute pairwise Jaccard similarity (or containment index) between all sketch pairs. The experiment succeeds if:

```
sim(sample_i_hg19, sample_i_hg38) > sim(sample_i_hg19, sample_j_hg19)
```

for all i ≠ j tissue pairs.

### Controls

- Replicate pairs (same sample, same reference): should have highest similarity
- Same tissue, same reference: high similarity
- Same tissue, different reference: high similarity (this is what we demonstrate)
- Different tissue, either reference: lower similarity

---

## Experiment B: Tissue-over-Species Clustering

### Design

Reproduce the clustering from Yue et al. (2014) and Lin et al. (2014) using minimizer sketches instead of alignment-based signal tracks or count matrices.

### Tissue and Assay Selection

Following Yue et al. (2014) and Lin et al. (2014), use the overlapping tissue types available in both human (ENCODE/Roadmap) and mouse (Mouse ENCODE):

| Tissue | Human (Roadmap) | Human GEO | Mouse (LICR) | Mouse GEO |
|---|---|---|---|---|
| Heart | E095 | GSM1127115 / GSE16256 | GSM1000093 | GSE49847 |
| Liver | E066 | GSM1127068 / GSE16256 | GSM1000140 | GSE49847 |
| Brain (cortex) | E067 | GSM1127069 / GSE16256 | GSM1000100 | GSE49847 |
| Spleen | E113 | GSM1127132 / GSE16256 | GSM1000138 | GSE49847 |

> **Note:** Lung is excluded from the initial H3K27ac analysis — H3K27ac ChIP-seq is absent from GSE49847 for lung (only H3K4me3/H3K4me1/Pol2/CTCF are available). Lung will be added in the planned H3K4me3 extension (see below). GSE49847 is the GEO Series for Mouse ENCODE LICR ChIP-seq data (Yue et al. 2014). Roadmap human data are under GEO Superseries GSE16256.

**Histone marks to include (in priority order):**

1. H3K4me3 — active promoters; used in all key papers
2. H3K27ac — active enhancers; strongest tissue discriminator (Roadmap 2015)
3. H3K4me1 — poised/active enhancers; used in Yue et al. and Roadmap 2015

**Planned extension — H3K4me3 across all tissues including lung:**
The initial analysis uses H3K27ac across 4 tissues (heart, liver, brain, spleen). A subsequent run will switch to H3K4me3, which is available in GSE49847 for all target tissues including lung, allowing a 5-tissue comparison and a cross-mark consistency check. H3K4me3 GSM accessions for Roadmap human samples (E095/E066) are already confirmed (GSM1127114, GSM1127067); mouse lung H3K4me3 GSMs need to be pulled from GSE49847.

### Expected Clustering Result

Based on prior published results, minimizer sketch similarity should recover:

```
(human_heart, mouse_heart) cluster together
(human_liver, mouse_liver) cluster together
...
All human tissues cluster together at a higher level than with their mouse counterparts
```

The species barrier should be visible at the top level of a dendrogram, with tissue identity dominating at finer resolution — matching the "tissue-over-species" hierarchy described in Yue et al. (2014) Figure 2 and Lin et al. (2014) Figure 1.

### Sketching tracks

Two parallel tracks will be run and compared:

| Track | Input to sketcher | Rationale |
|---|---|---|
| **Peak-called (primary)** | FASTA extracted from peak BED regions | Isolates tissue-specific signal; matches what prior papers used; expected to be necessary for clean clustering |
| **Alignment-only (secondary)** | FASTA from full BAM | Tests whether sketching alone is sufficient without peak calling; retained as a check, not the expected success path |

Peak calling will use MACS2 (narrow peaks for H3K4me3; broad peaks for H3K27ac and H3K4me1) with the corresponding input control where available. Sequences are then extracted with `bedtools getfasta`.

If both tracks recover the expected clustering, that would be a positive result for the alignment-only approach. If only the peak-called track succeeds, that informs the recommended usage of the sketching tool.

### Metric

Pairwise minimizer sketch similarity (Jaccard / containment) computed across all human and mouse samples, then hierarchical clustering (complete linkage) used to generate a dendrogram. Cluster recovery assessed by comparing to ground-truth tissue labels. Metrics computed independently for each sketching track.

---

## Data Acquisition Plan

### Step 1: Download raw or aligned reads from ENCODE portal

```bash
# Example: ENCODE portal metadata query for H3K27ac human tissues
curl -L "https://www.encodeproject.org/metadata/?type=Experiment&assay_title=Histone+ChIP-seq&target.label=H3K27ac&biosample_ontology.organ_slims=heart&files.file_type=fastq" \
  -o metadata_h3k27ac_heart.tsv
```

For Roadmap data, use the EBI ENA or NCBI SRA entries linked from GEO Superseries GSE16256.

### Step 2: Align (if using FASTQ)

- Human: align to GRCh38 (primary) and GRCh37 (Experiment A only) using BWA-MEM2
- Mouse: align to mm10 (GRCm38; mm39 is not available on this cluster)

### Step 3: Peak calling (Experiment B only)

Call peaks with MACS2 on each aligned BAM:

```bash
# Narrow peaks (H3K4me3)
macs2 callpeak -t sample.bam -c input.bam -f BAM -g hs --nomodel --extsize 200 \
  -n sample_h3k4me3 --outdir peaks/

# Broad peaks (H3K27ac, H3K4me1)
macs2 callpeak -t sample.bam -c input.bam -f BAM -g hs --broad \
  -n sample_h3k27ac --outdir peaks/
```

Use `-g mm` for mouse samples. Input controls should be matched by experiment where available; fall back to pooled input if not.

Peak BED files from this step are used as regions for `bedtools getfasta` in the peak-called track.

### Step 4: Sketch

Two tracks for Experiment B; one track for Experiment A:

- **Alignment-only track:** Extract FASTA from full BAM (e.g., via `samtools fasta`), then sketch with hammock mode D.
- **Peak-called track (Experiment B):** Extract FASTA from peak regions via `bedtools getfasta`, then sketch with hammock mode D.

Sketching parameters (k, w) are being selected via the companion `mode-d-optimality` sweep experiment; see that experiment's results before finalising values in `config/config.yaml`.

### Step 5: Compute pairwise similarities

Generate all-vs-all similarity matrix, separately for each track.

### Step 6: Cluster and visualize

Hierarchical clustering, PCA, and/or t-SNE/UMAP to visualize groupings.

---

## Reproducibility Infrastructure

- **Workflow manager:** Snakemake (primary) with Nextflow alternative
- **Containerization:** Singularity (cluster-compatible)
- **Cluster:** SLURM; bigmem partition for large BAM operations
- **Version control:** Git; all documents in Markdown (LaTeX migration planned)
- **Code access:** Claude Code extension on Cursor

See `workflow/Snakefile` for pipeline definitions (Snakemake + SLURM; no Nextflow config).

---

## Empirical Findings (2026-05-11) & Redesign Options for Experiment B

### Observation

After running the full `(k, w)` sweep (19 valid pairs after dropping `(5, 5)` for runtime; both `broad` and `narrow` peaks; refactored `hammock` from `hammock_claude` installed in the `claude-ref-comparison` conda env), we evaluated each `cluster_assignments.tsv` against tissue and species labels using Adjusted Rand Index:

| Metric | Result |
|---|---|
| Sweep cells evaluated | 38 (19 k,w × 2 peak types) |
| Cells where tissue ARI > species ARI | **0** |
| Typical species ARI | 0.706 (humans form one cluster; mouse splits 2-way) |
| Typical tissue ARI | **−0.296** (clustering is *anti-correlated* with tissue) |
| Best tissue ARI observed | −0.061 (at k=10, w≥20; still negative) |

**The exp_b hypothesis as currently designed cannot be validated by this pipeline.** No parameter combination recovers tissue-over-species clustering. This is consistent across both broad and narrow peak calling.

### Why this happens

Peak FASTAs are extracted as `bedtools getfasta -fi {ref}.fa -bed {peaks}.bed`. The output sequence is whatever sits under the called peaks **in the source genome**. Two consequences:

1. Human samples produce only human-genome k-mers; mouse samples produce only mouse-genome k-mers.
2. Human ↔ mouse 15-mer overlap is ~0 even at orthologous, evolutionarily-conserved enhancers (~80 Mya divergence; sequence has drifted past exact-match threshold).

Hammock's k-mer-set similarity therefore sees species composition, not tissue identity. A minhash/HLL sketch over raw peak FASTA **cannot see orthology** — it only sees k-mer-set membership.

The known biological result of tissue-over-species clustering in the literature ([@Yue2014], [@Lin2014], [@Roadmap2015]) does not contradict this — those analyses used peak *positions* in a lifted/orthologous coordinate system, signal *intensity* at orthologous regions, or hand-curated conserved enhancer sets. None of those representations is what a raw k-mer sketch can capture.

### Three redesign paths for Experiment B

Each redesign factors species composition out of the input, so the residual signal is tissue identity.

| Option | Mechanism | Pipeline change | Trade-offs |
|---|---|---|---|
| **B1. Liftover-based** | Map mouse peaks to GRCh38 coordinates (e.g., UCSC `liftOver` with `mm10ToHg38.over.chain`), then extract FASTA from GRCh38 at the lifted intervals. All samples' FASTAs now come from one genome. | New rule before `peaks_to_fasta`: `liftOver mouse_peaks.bed mm10ToHg38.over.chain mouse_peaks.lifted.bed unmapped.bed`. Use lifted BED for mouse track. | Liftover drops 15–30% of mouse peaks that lack a 1:1 GRCh38 ortholog. Lifted peaks span the *human* genome's version of the conserved region — sequence is fully human. Most direct way to test the original hypothesis, but mouse-specific peaks are lost from the analysis. |
| **B2. Orthologous-peak filter** | Restrict both species' peaks to a curated set of mammalian-conserved enhancer regions (e.g., VISTA enhancers, phastCons-conserved blocks, or ENCODE cCRE conservation track). | New rule: intersect peaks with conservation BED before `peaks_to_fasta`. Peaks become a subset; downstream unchanged. | Filters out the most species-specific peaks (which would have dominated the sketch). Both species' inputs are at conserved genomic regions but still each in its own genome — sequence still divergent. Likely partial fix, may need to combine with B1. Easiest to implement (one `bedtools intersect` rule). |
| **B3. Compositional control** | Subtract genome-wide background k-mer frequency from each sketch before computing similarity. Removes species composition signal at the sketch level. | Requires C++-side changes in `hammock_claude` to support background subtraction — significant scope. Or: post-process the all-vs-all matrix to regress out within-species mean similarity before clustering. | Doesn't lose any peaks. Mathematically cleanest. But the post-processing variant is a one-off statistical correction (residualize after similarity), not a property of the sketch itself — fine for one figure but doesn't generalize. Either form needs care to avoid removing too much signal. |

### Decision steps

1. **Within-species tissue signal — sanity check completed 2026-05-11.** Across all 38 sweep cells, within-species pair rankings are highly consistent:
   - **Human**: heart–lung is the most similar pair in 34/38 cells (89%); the residual 4 cells flip heart-lung ↔ liver-lung.
   - **Mouse**: liver–lung is the most similar pair in 34/38 cells (89%); the rest are minor variants.
   - Within-species similarity spread is 5–10% of max (broad: human 0.088, mouse 0.052; narrow: human 0.110, mouse 0.063) — comparable in scale to the cross-reference effect size from exp_a. So **the sketch sees real tissue-axis signal within species**; it's just numerically dominated by the cross-species composition difference. Greenlights pursuing the redesign options below.
   - **Caveat:** humans pick `heart–lung` while mice pick `liver–lung`, despite tissue-tissue biology being conserved across mammals. Some of the "tissue signal" we're picking up is likely technical (lab/depth/peak-calling differences between ENCSR samples), not purely biological. Sample provenance is worth a sanity audit before drawing biological conclusions, but the orthogonal demonstration that within-species rankings *are stable across (k, w)* still validates that the sketch is sensitive to non-species axes.

2. **If within-species tissue signal is present (it is), pursue B2 (orthologous-peak filter) as the cheapest first attempt.** Single `bedtools intersect` rule; everything downstream is unchanged. If B2 alone recovers tissue-over-species clustering, ship it. If B2 helps but isn't enough, escalate to B1 (liftover) which factors species composition out completely.

3. **B3 is a longer-horizon option** — useful if both B1 and B2 fail, or if the project later wants a sketch-level solution that works without peak-call preprocessing. Not the right first move.

4. **Report exp_b as-is alongside the redesign.** The current negative result *is* a publishable finding ("k-mer sketches of peak FASTAs cannot do cross-species tissue identification without genome-aware preprocessing") — it characterizes what hammock can and can't do. The redesign result then shows what additional preprocessing buys.

---

## Open Questions / TODOs

- [ ] Run `scripts/fetch_accessions.py` to confirm SRA run IDs for Roadmap E095/E066/E096 and Mouse ENCODE GSE49847; update `sra_map` in `config/config.yaml`
- [ ] Select final (k, w) parameters from `mode-d-optimality` sweep results and update `kmer_sizes` / `window_sizes` in `config/config.yaml`
- [ ] Add MACS2 peak calling rules to `workflow/Snakefile` (narrow for H3K4me3, broad for H3K27ac/H3K4me1); wire peak-called track as primary for Experiment B
- [ ] Confirm input control BAM availability for each sample; update `config/config.yaml` with input accessions
- [ ] Establish storage requirements for full FASTQ + BAM + FASTA dataset
- [x] **Exp B redesign step 1 (DONE 2026-05-11):** within-species tissue-axis sanity check confirmed tissue signal exists (humans 34/38 cells: heart-lung closest; mice 34/38 cells: liver-lung closest; spread ~5–10% of max within species). The differing human/mouse "closest pair" suggests technical (lab/depth) contributions on top of biology — audit sample provenance before publishing biological conclusions.
- [x] **Exp B redesign step 2 (DONE 2026-05-11, reframed 2026-05-12):** implemented B2 (orthologous-peak filter via `bedtools intersect` with UCSC phastCons elements: hg38 100-way, mm10 60-way). Added rules: `prep_conservation_bed`, `filter_peaks_to_conserved`, `peaks_conserved_to_fasta`, `exp_b_fasta_list_conserved`, `exp_b_hammock_conserved`, `exp_b_cluster_and_plot_conserved`. Outputs under `results/exp_b_conserved/`. **Result at N=12 (2026-05-12): null.** 0 of 40 (k, w) × peak-type cells achieve tissue-dominant clustering. The earlier 6-sample result (1/38 at broad k=5 w=8 with ARI tissue 0.444 vs species −0.176) **did not replicate** at 12 samples — that cell is now species-dominant. Strong evidence the prior hit was an underpowered fluke. Interpretation: conservation filtering reduces species-specific *sequence content* but does not align orthologous peaks to each other, so the sketch still scores composition rather than orthology.
- [ ] **Exp B redesign step 3 (IN PROGRESS 2026-05-12):** implement B1 (coordinate-space liftover). Plan: download UCSC `mm10ToHg38.over.chain.gz`; liftover mouse mm10 peaks → GRCh38 with UCSC `liftOver` (drop unmapped peaks, ~30% per published rates); extract FASTA from GRCh38 at the lifted intervals; re-run hammock all-vs-all on the same 12-sample design. Outputs under `results/exp_b_lifted/`. Both species' inputs come from a single genome (GRCh38), so coordinate-space species composition is entirely factored out — only orthology + peak landscape remain in the input. Decisive test: if B1 also produces null result, the literature's "tissue beats species" claim is properly localized to *intensity* information (peak height / coverage), not sequence content of peaks.
- [ ] **Exp A finalization:** pick a single `(k, w)` operating point from `results/exp_a/{broad,narrow}/sweep_effect_size.png` based on effect size (`median(cross-ref) − median(diff-tissue)`); record choice + rationale here.
- [ ] **`(5, 5)` revisit:** broad-peak hammock at `k=5, w=5` exceeded 4h runtime and was dropped from the sweep. Confirm it's not worth profiling/optimizing — adjacent `(5, 8)` had effect size ≈ 0.03, well below the `k=10` peak (~0.10), so `(5, 5)` was unlikely to be the optimum anyway.
- [x] **Sample expansion for Exp B (DONE 2026-05-11, phase 1):** added 3 adult-matched tissues (spleen, small intestine, testis), bringing exp_b to 6 tissues × 2 species = 12 samples (6 cross-species pairs vs prior 3). Added accessions:
  - human: `ENCSR436JNB` (spleen, Broad, 76bp SE, 2 repls, 61F), `ENCSR655XLM` (sm.int., UCSD, 50bp SE, 1 repl, 30F — flag in methods), `ENCSR136ZQZ` (testis, Broad, 76bp SE, 2 repls, 37M)
  - mouse (all Bing Ren UCSD, 36bp SE, 2 repls, 2mo B6NCrl): `ENCSR000CDJ` (spleen), `ENCSR000CCQ` (sm.int.), `ENCSR000CCU` (testis)
  - All have `possible_controls` set in ENCODE → nf-core/chipseq pairs them automatically.
- [ ] **Sample expansion for Exp B (phase 2 — brain regions, deferred):** mouse adult H3K27ac is available for cerebellum (`ENCSR000CDC`), cortical plate (`ENCSR000CDD`), olfactory bulb (`ENCSR000CCE`) — all UCSD B6NCrl 2mo, matching existing mouse protocol. **Blocker:** need matched-age human brain H3K27ac (ENCODE portal not yet queried; many brain ChIP-seq experiments are fetal/embryonic). Re-run sample-identification queries restricting to `life_stage=adult` for human brain anatomical terms (cerebellum, frontal cortex, etc.).
- [ ] **Sample expansion for Exp B (phase 3 — replicates):** current design is still 1 ENCSR per tissue per species. Adding 2–3 reps per cell (different donors / mouse strains) would push to 24–36 samples and sharpen tissue-vs-species ARI.
- [ ] **Sample-identification helper (paired with sample expansion):** query the ENCODE portal API for matched-tissue H3K27ac ChIP-seq across human + Mouse ENCODE; filter on data-quality flags (audit status, read depth, peak count), confirm matched input controls; emit an updated `config/config.yaml` `samples:` block. Likely scripted, not manual.

---

## References

See `references.bib` for full bibliographic data.

Key references: [@Yue2014], [@Lin2014], [@Roadmap2015], [@ENCODE2020], [@BarbosaMorais2012], [@Merkin2012]

