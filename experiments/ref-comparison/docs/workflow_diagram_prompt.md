Create a detailed workflow figure suitable for a computational biology methods paper (Nature/Genome Biology style).

The figure shows a two-branch pipeline for validating minimizer-based sequence sketching on ChIP-seq data.

## Layout
Vertical flow, top to bottom. Shared preprocessing steps in a center column, then the pipeline forks left (Experiment A) and right (Experiment B) after the FASTA extraction step. Use a clean sans-serif font (Helvetica or Arial). Color scheme: muted blues and grays for data/process boxes, warm orange for the hammock sketching step (the key method being validated), green for output/result boxes.

## Data sources (top of figure)
Two labeled data source boxes with dataset tags:
- "Human ChIP-seq (H3K27ac)" — Roadmap Epigenomics GSE16256; 6 tissues: heart, liver, lung, brain, spleen, small intestine
- "Mouse ChIP-seq (H3K27ac)" — Mouse ENCODE LICR GSE49847; 4 tissues: heart, liver, brain, spleen
Both feed downward into the shared preprocessing column.

## Shared preprocessing steps (center column)
Sequential boxes with labeled arrows showing data transformation:
1. **Download** — fasterq-dump (SRA Tools 3.0); output: paired-end FASTQ
2. **Align** — BWA-MEM2 v2.2 → SAMtools sort/index; output: sorted BAM
   - Show three reference genome labels branching off: GRCh37, GRCh38 (human), mm10 (mouse)
3. **Peak calling** — MACS2 v2.2 (--nomodel, BAMPE mode); output: narrowPeak BED
4. **Sequence extraction** — bedtools getfasta; output: peak FASTA
   - This is the fork point — add a horizontal dashed line here labeled "per-experiment sketching"

## Left branch — Experiment A: Cross-Reference Robustness
Title label: "Exp A: Reference Robustness"
Steps:
- Hammock sketching box (orange): minimizer sketching, k∈{3,5,7}, w∈{5,10,15,20}, precision=24; inputs: all sample FASTAs from GRCh37 AND GRCh38 combined
- All-vs-all Jaccard similarity matrix (heatmap icon, 6×6 grid schematic)
- Output (green): boxplot icon — "Same-tissue cross-ref vs. different-tissue similarity; Wilcoxon test"

## Right branch — Experiment B: Tissue-over-Species Clustering
Title label: "Exp B: Tissue Identity"
Steps:
- Hammock sketching box (orange): same parameters; inputs: human (GRCh38) + mouse (mm10) FASTAs
- All-vs-all Jaccard similarity matrix (9×9 grid schematic with faint tissue grouping lines)
- Output (green): dendrogram icon — "Hierarchical clustering (complete linkage); PCA"; label citing Yue et al. 2014, Lin et al. 2014

## Style details
- Process steps: rounded rectangles with a subtle drop shadow
- Data objects (FASTQ, BAM, BED, FASTA, CSV): parallelograms or document icons
- Arrows: directed, labeled with the data type being passed
- Parameter sweep annotation: small bracket near the hammock boxes labeled "11 (k,w) pairs"
- Container annotations: small pill badges on each process step reading "Singularity" with the tool version
- Figure panel labels: (A) and (B) in bold upper-left of each branch
- No decorative gradients; flat design with thin borders
- Resolution: 300 dpi, suitable for print; aspect ratio approximately 2:3 (wide portrait)
