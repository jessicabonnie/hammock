# Slide Deck Prompt: Data Selection for Minimizer Sketch Validation

---

## Slide 1 — Overview: Two Complementary Validation Experiments

**Layout:** Two-column with a title banner at top.

**Title:** "Validation Strategy: Two Independent Experiments"

**Left column — Experiment A:**
- Heading: "Reference Robustness"
- Body: Same biological sample, two reference genomes (GRCh37 vs. GRCh38). Does the sketch "see through" reference version differences?
- Icon: two reference genome icons (e.g., chromosome diagrams) with a double-headed arrow labeled "same sample"

**Right column — Experiment B:**
- Heading: "Tissue-over-Species Clustering"
- Body: Human and mouse samples from matched tissues. Does sketch similarity reproduce the known biological hierarchy — tissue identity stronger than species identity?
- Icon: a small 2×2 grid of tissue icons (heart, liver, brain, spleen) with human and mouse silhouettes as row/column headers

**Bottom bar:** A single sentence tying them together — "Both experiments operate directly on aligned reads (BAM/FASTQ), not on peak calls, for a clean test of the sketching method."

---

## Slide 2 — Experiment A: Data Source and Logic

**Layout:** Top title, then a 3-row comparison table in the center, with an explanatory annotation panel on the right.

**Title:** "Experiment A — Cross-Reference Data: Roadmap Epigenomics H3K27ac"

**Central figure — comparison matrix (hand-drawn style or simple table):**

```
              GRCh37 (hg19)    GRCh38 (hg38)
Heart (E095)      ●                ●
Liver (E066)      ●                ●
Lung  (E096)      ●                ●
```

Dots are filled circles. Arrows connect same-tissue dots across columns (solid, labeled "should be similar"). Arrows connecting different-tissue dots within a column (dashed, labeled "should be less similar"). This is the core prediction diagram.

**Right annotation panel:**
- "Why H3K27ac?" → Roadmap 2015 identified it as the strongest tissue-type discriminator among histone marks — maximizes contrast between positive and negative pairs
- "Why these tissues?" → Biologically distinct enough to serve as tight negative controls; all three have both hg19 and hg38 alignments available via ENCODE portal

**Bottom source note:** Data from Roadmap Epigenomics Consortium (2015, *Nature*); GEO Superseries GSE16256; ENCODE portal accessions E095, E066, E096. Downloaded as FASTQ via SRA; aligned to GRCh37 and GRCh38 with BWA-MEM2.

---

## Slide 3 — Experiment B: Tissue × Species Data Grid

**Layout:** Large central grid figure, with a sidebar listing mark selection rationale.

**Title:** "Experiment B — Cross-Species Data: Human (Roadmap) + Mouse ENCODE (LICR)"

**Central figure — filled data availability grid:**

A 4-row × 2-column grid (tissues as rows, species as columns). Each cell contains a small icon for the histone mark available.

```
              Human (Roadmap/ENCODE)      Mouse (ENCODE LICR)
Heart              H3K27ac ✓                  H3K27ac ✓
Liver              H3K27ac ✓                  H3K27ac ✓
Brain (cortex)     H3K27ac ✓                  H3K27ac ✓
Spleen             H3K27ac ✓                  H3K27ac ✓
```

Below the grid, add a grayed-out row:
```
Lung               H3K27ac ✓                  H3K27ac ✗  ← excluded (unavailable in GSE49847)
                   H3K4me3 ✓                  H3K4me3 ✓  ← planned extension
```

**Right sidebar — Why these marks?**
- H3K27ac: active enhancers; strongest tissue discriminator (Roadmap 2015); primary mark for this analysis
- H3K4me3: active promoters; used in Yue et al. 2014 and Lin et al. 2014; planned cross-mark consistency check
- H3K4me1: poised/active enhancers; secondary validation

**Bottom source note:** Human data — Roadmap Epigenomics Consortium; GEO GSE16256; ENCODE accessions E095, E066, E067, E113. Mouse data — Mouse ENCODE LICR (Yue et al. 2014); GEO GSE49847; GSM accessions GSM1000093, GSM1000140, GSM1000100, GSM1000138.

---

## Slide 4 — Why This Data? Grounding in Prior Literature

**Layout:** Timeline / evidence chain running left to right, with the key publications as nodes.

**Title:** "Data Choices Are Anchored to Published Ground Truth"

**Main figure — horizontal evidence chain:**

```
Yue et al. 2014          Lin et al. 2014          Roadmap 2015
(Nature)                 (Nature)                  (Nature)
Mouse ENCODE             H3K27ac/RNA-seq           111 human epigenomes
tissue clustering        15 tissues H+M            cluster by tissue/lineage
        │                       │                        │
        └───────────────────────┴────────────────────────┘
                                │
                    "Tissue identity > species identity"
                    established in matched tissues:
                    heart, liver, brain, spleen
                                │
                                ▼
              THIS EXPERIMENT: reproduce with minimizer sketches
              using the same datasets → direct replication test
```

**Below the chain:** Two bullet points
- "Using the exact datasets from these papers makes our result a direct replication, not a loose analogy"
- "If minimizer sketching fails to recover the known hierarchy, it is a falsifiable failure — not an ambiguous result"

---

## Slide 5 — Data Acquisition Pipeline

**Layout:** Vertical flowchart, center-aligned, with annotations on the left for human data and right for mouse data.

**Title:** "How the Data Was Obtained"

**Central flowchart (5 steps):**

```
[ENCODE Portal / GEO SRA]
         │
         ▼
  Download FASTQ via SRA
  (prefetch + fasterq-dump)
         │
         ▼
  Align to reference genome
  (BWA-MEM2; human→hg38+hg19,
   mouse→mm10)
         │
         ▼
  Extract sequences from BAM
  (bedtools getfasta → FASTA)
         │
         ▼
  Apply minimizer sketching
  (hammock mode D; k,w from
   mode-d-optimality sweep)
         │
         ▼
  Pairwise sketch similarity
  (Jaccard / containment)
```

**Left annotation (human):**
- ENCODE portal metadata API query for H3K27ac + tissue + file type
- GEO Superseries GSE16256 (Roadmap)
- Roadmap IDs: E095, E066, E067, E096, E113

**Right annotation (mouse):**
- GEO Series GSE49847 (Mouse ENCODE LICR)
- GSM accessions confirmed via `scripts/fetch_accessions.py`
- Reference: mm10 (GRCm38); mm39 not available on cluster

**Bottom note:** Workflow managed by Snakemake + Singularity + SLURM (bigmem partition for alignment); all steps are fully reproducible from the repo.

---

## General Design Notes

- Use a consistent color scheme: one color for human samples, a contrasting color for mouse samples, a third for "planned/future" items
- Keep tissue icons consistent across slides (same icon for heart on every slide, etc.)
- Grayed-out items (lung exclusion, planned extensions) should use a visually distinct muted style
- Source citations can be footnote-sized at the bottom of each slide
- Avoid dense text — each slide should have at most one short explanatory paragraph; the rest should be visual
