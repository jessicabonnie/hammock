# Experiment Design: Minimizer-Based Sketching on ENCODE Epigenomic Data

**Project:** Validation of minimizer-based sequence sketching across reference genomes and biological contexts  
**Status:** Workflow built; pending accession confirmation and parameter selection  
**Last Updated:** 2026-04

---

## Overview

This document describes the experimental design for two complementary validation experiments using minimizer-based sequence sketching applied to ENCODE and Roadmap Epigenomics ChIP-seq data:

1. **Experiment A — Reference Robustness:** Demonstrate that minimizer sketches of the same biological sample aligned to GRCh37 vs. GRCh38 are more similar to each other than to sketches from a different sample aligned to either reference.
2. **Experiment B — Tissue-over-Species Clustering:** Reproduce the known result that epigenomic profiles cluster by tissue type rather than species when comparing human and mouse samples, using minimizer sketching as the similarity metric.

Both experiments operate directly on aligned sequence files (BAM/FASTQ), not on peak calls or signal tracks, making them a clean test of the sketching approach independent of any downstream peak-calling pipeline.

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

| Tissue | Human Source | Human GEO/Accession | Mouse Source | Mouse GEO/Accession |
|---|---|---|---|---|
| Heart | Roadmap E095 (left ventricle) | GSE16256 | Mouse ENCODE | GSE49847 (LICR) |
| Liver | Roadmap E066 | GSE16256 | Mouse ENCODE | GSE49847 (LICR) |
| Lung | Roadmap E096 | GSE16256 | Mouse ENCODE | GSE49847 |
| Brain (frontal cortex / cerebellum) | Roadmap E067/E071 | GSE16256 | Mouse ENCODE | GSE49847 |
| Spleen | Roadmap E113 | GSE16256 | Mouse ENCODE | GSE49847 |
| Small intestine | Roadmap E109 | GSE16256 | Mouse ENCODE | GSE49847 |

> **Note:** GSE49847 is the GEO Series for Mouse ENCODE LICR ChIP-seq data used in the Lin et al. (2014) analysis. The Roadmap human data are under GEO Superseries GSE16256. See bibliography for full citations.

**Histone marks to include (in priority order):**

1. H3K4me3 — active promoters; used in all key papers
2. H3K27ac — active enhancers; strongest tissue discriminator (Roadmap 2015)
3. H3K4me1 — poised/active enhancers; used in Yue et al. and Roadmap 2015

### Expected Clustering Result

Based on prior published results, minimizer sketch similarity should recover:

```
(human_heart, mouse_heart) cluster together
(human_liver, mouse_liver) cluster together
...
All human tissues cluster together at a higher level than with their mouse counterparts
```

The species barrier should be visible at the top level of a dendrogram, with tissue identity dominating at finer resolution — matching the "tissue-over-species" hierarchy described in Yue et al. (2014) Figure 2 and Lin et al. (2014) Figure 1.

### Metric

Pairwise minimizer sketch similarity (Jaccard / containment) computed across all human and mouse samples, then hierarchical clustering (complete linkage) used to generate a dendrogram. Cluster recovery assessed by comparing to ground-truth tissue labels.

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

### Step 3: Sketch

Extract sequences from aligned BAM files as FASTA (via `bedtools getfasta`), then apply minimizer sketching (hammock mode D) to the FASTA files. Sketching parameters (k, w) are being selected via the companion `mode-d-optimality` sweep experiment; see that experiment's results before finalising values in `config/config.yaml`.

### Step 4: Compute pairwise similarities

Generate all-vs-all similarity matrix.

### Step 5: Cluster and visualize

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

## Open Questions / TODOs

- [ ] Run `scripts/fetch_accessions.py` to confirm SRA run IDs for Roadmap E095/E066/E096 and Mouse ENCODE GSE49847; update `sra_map` in `config/config.yaml`
- [ ] Select final (k, w) parameters from `mode-d-optimality` sweep results and update `kmer_sizes` / `window_sizes` in `config/config.yaml`
- [ ] Decide whether to sketch from raw FASTQ, aligned BAM, or peak-called BED (current plan: BAM → FASTA via bedtools)
- [ ] Establish storage requirements for full FASTQ + BAM + FASTA dataset

---

## References

See `references.bib` for full bibliographic data.

Key references: [@Yue2014], [@Lin2014], [@Roadmap2015], [@ENCODE2020], [@BarbosaMorais2012], [@Merkin2012]

