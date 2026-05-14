# ref-comparison

Validates minimizer-based sketching (Mode D) on ENCODE H3K27ac ChIP-seq for
**cross-reference robustness**: the same biological sample aligned to GRCh37,
GRCh38, and CHM13 should sketch more similarly across references than across
tissues on any one reference.

Migrated from `hammock/experiments/claude-ref-comparison/` (2026-05-12). The
former Experiment B (tissue-over-species) was split off into separate
experiments (`mus-homo`, `primate-phylogeny`, `man-monkey`) because it
doesn't actually need the cross-reference alignment work that justifies
this directory's existence.

## Two-stage pipeline

1. **nf-core/chipseq** (upstream) — download FASTQs → align to each reference → peak-call. Run via `scripts/run_nfcore.sh`. This is the only piece of the project that genuinely needs a multi-reference alignment pipeline; everything else downstream of peak calling can use pre-called BEDs from ENCODE directly.
2. **Snakemake** (downstream) — peaks-to-FASTA → hammock (k, w) sweep → plots.

## Key files

- `docs/experiment_design.md` — full design rationale, accession tables, success criteria
- `docs/paper_outline.md` — current paper outline (Exp A only)
- `config/config.yaml` — sample metadata, nf-core output paths, sketch parameters
- `config/nfcore_samplesheet_human.csv` — nf-core input: 3 ChIP-seq samples + matched controls × 3 references
- `workflow/Snakefile` — downstream pipeline (peaks → FASTA → hammock → plots)
- `workflow/slurm_profile/config.yaml` — SLURM cluster profile for Snakemake
- `scripts/run_nfcore.sh` — launches nf-core/chipseq for GRCh38, GRCh37, CHM13

## Running

```bash
# Step 1: Run nf-core/chipseq (handles download, alignment, peak calling)
bash scripts/run_nfcore.sh

# Step 2: Run downstream Snakemake (peaks → FASTA → hammock → plots)
conda activate claude-ref-comparison
snakemake --dryrun                              # check job plan
snakemake --profile workflow/slurm_profile/     # submit to SLURM
```

## Sample design

| Tissue | Sample | References |
|--------|--------|------------|
| Heart | ENCSR175ABH | GRCh37, GRCh38, CHM13 |
| Liver | ENCSR864OOO | GRCh37, GRCh38, CHM13 |
| Lung  | ENCSR954JMZ | GRCh37, GRCh38, CHM13 |

The human samplesheet (`config/nfcore_samplesheet_human.csv`) contains additional
biosamples (spleen, sm. intestine, testis, etc.) that were used in the former
Experiment B; they're orphaned here but kept on disk so nf-core's `-resume` cache
stays warm. They're not referenced by the Snakefile.
