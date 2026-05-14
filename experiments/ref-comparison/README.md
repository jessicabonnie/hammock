# ref-comparison — Cross-Reference Robustness

Validates minimizer-based sequence sketching (hammock Mode D) on ENCODE H3K27ac
ChIP-seq by asking whether **the same biological sample aligned to different
human reference genomes** (GRCh37, GRCh38, CHM13) sketches more similarly across
references than across tissues on any one reference.

This experiment was originally part of `claude-ref-comparison` together with a
tissue-over-species experiment ("Exp B"). The tissue-over-species work has been
split into separate experiments (`mus-homo`, `primate-phylogeny`, `man-monkey`)
because it doesn't actually need the multi-reference alignment pipeline that
justifies this directory's existence.

---

## Repository structure

```
ref-comparison/
├── docs/
│   ├── experiment_design.md   ← design rationale, accessions, success criteria
│   ├── paper_outline.md       ← paper outline (Exp A only)
│   └── references.bib         ← BibTeX bibliography
├── workflow/
│   ├── Snakefile              ← downstream pipeline (peaks → FASTA → hammock → plots)
│   └── slurm_profile/
│       └── config.yaml        ← SLURM cluster profile for Snakemake
├── config/
│   ├── config.yaml            ← samples, refs, sketch parameters
│   ├── exp_a_metadata.tsv     ← sample × ref metadata for plotting
│   ├── confirmed_accessions.tsv
│   └── nfcore_samplesheet_human.csv
├── scripts/
│   ├── fetch_accessions.py    ← resolve GEO/SRA accession numbers
│   ├── run_nfcore.sh          ← launches nf-core/chipseq per reference
│   ├── exp_a_validate_plot.R  ← Wilcoxon + 2-panel boxplot + heatmap
│   └── exp_a_sweep_summary.R  ← (k × w) effect-size heatmap across the sweep
└── README.md
```

Hammock + nf-core outputs are stored at `/vast/blangme2/jbonnie/hammock/claude-ref-comparison/`
(kept on the legacy path so existing CSVs / plots remain valid without re-running).

---

## Two-stage pipeline

1. **nf-core/chipseq** (upstream) — download FASTQs, align to each reference, peak-call. Multi-reference re-alignment is the part of the workflow that **needs** a full pipeline; the other experiments in `hammock_claude/experiments/` can use pre-called peak BEDs from ENCODE directly.
2. **Snakemake** (downstream, this directory) — peaks-to-FASTA → hammock all-vs-all → (k, w) sweep → plots.

---

## Quickstart

```bash
# Step 1: Upstream — run nf-core/chipseq across references
bash scripts/run_nfcore.sh                       # all references (GRCh38, GRCh37, CHM13)
bash scripts/run_nfcore.sh GRCh38                # one reference only

# Step 2: Downstream — Snakemake (peaks → FASTA → hammock → plots)
conda activate claude-ref-comparison
snakemake --dryrun                               # check job plan
snakemake --profile workflow/slurm_profile/      # submit to SLURM
```

---

## Sample design

| Tissue | ENCODE accession | References |
|--------|------------------|------------|
| Heart | ENCSR175ABH | GRCh37, GRCh38, CHM13 |
| Liver | ENCSR864OOO | GRCh37, GRCh38, CHM13 |
| Lung  | ENCSR954JMZ | GRCh37, GRCh38, CHM13 |

3 samples × 3 references = 9 sample × ref combinations. Two peak types
(MACS3 broad + narrow). Parameter sweep over `k ∈ {5, 8, 10, 15, 20}`
× `w ∈ {5, 8, 10, 15, 20, 30}` filtered to `w ≥ k` (20 cells).

---

## Compute requirements

| Step | CPUs | RAM | Est. walltime |
|---|---|---|---|
| nf-core download + alignment (per sample × ref) | 16 | 32 GB | ~4 h |
| hammock per (k, w) | 8 | 16 GB | 10 min – 12 h |
| R plotting (per cell, sweep summary) | 1 | 8 GB | <2 min |
