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

## Environment

```bash
ml anaconda && conda activate hammock   # always run this at the start of a terminal session
conda env create -f environment.yaml    # first time only (env name: claude-ref-comparison)
```

`environment.yaml` is the source of truth — update it and re-run `conda env update -f environment.yaml --prune` when dependencies change.

## Hammock output columns

When using minimizer (mode D) output CSVs, always use:
- `hash_with_ends_similarity` — primary metric (minimizers + flanking end k-mers)
- `hash_similarity` — secondary metric (minimizers only)

**Never use `jaccard_similarity`** — that column is legacy and should be ignored even if present.

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
