# claude-ref-comparison

Validates minimizer-based sketching on ENCODE/Roadmap ChIP-seq data via two experiments:

- **Experiment A** — same biological sample aligned to GRCh37 vs GRCh38 should sketch more similarly than different tissues on either reference
- **Experiment B** — minimizer sketch similarity should reproduce tissue-over-species clustering (Yue et al. 2014, Lin et al. 2014, Roadmap 2015)

## Key files

- `docs/experiment_design.md` — full design rationale, accession tables, success criteria
- `docs/references.bib` — BibTeX bibliography for all cited papers
- `config/config.yaml` — sample metadata, SRA IDs, reference paths, sketch parameters
- `workflow/Snakefile` — full pipeline (download → align → sketch → similarity → plots)
- `workflow/slurm_profile/config.yaml` — SLURM cluster profile

## Before running the workflow

1. Run `python scripts/fetch_accessions.py` to confirm SRA run IDs for Mouse ENCODE samples (GSE49847)
2. Update `sra_map` in `config/config.yaml` with confirmed SRR IDs
3. Update reference FASTA paths in `config/config.yaml` for your cluster

## Running

```bash
# Dry run
snakemake --dryrun --forceall

# Local
snakemake --cores 8 --use-singularity

# SLURM cluster
snakemake --profile workflow/slurm_profile/ --use-singularity
```

## Gotchas

- Alignment jobs need `bigmem` partition — already set in `slurm_profile/config.yaml`
- `sketch_bin` in `config/config.yaml` must point to the built sketching binary before any sketch rules will run
- Mouse ENCODE LICR sample GSM accessions (GSE49847) need manual confirmation — placeholders are in the config
