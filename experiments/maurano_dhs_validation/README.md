# maurano_dhs_validation — Modes A/B/C/D vs bedtools on fetal-tissue DHS

How accurately do Hammock's four modes (with HLL backing) recover bedtools
pairwise Jaccard on the Maurano et al. 2012 fetal-tissue DNase-seq corpus,
and what parameters are optimal?

This is the `hammock_claude` recreation of `dnase1-hypersensitivity` from the
original `hammock` repo. The headline result we expect to reproduce: Mode D
at k=10, w=20–50, p=24 hits Pearson r ≈ 0.998 against bedtools with perfect
clustering reproduction (ARI=NMI=1, RF=0) on the 20-sample x 20-sample
matrix.

## Corpus

Maurano et al. 2012 ("Systematic localization of common disease-associated
variation in regulatory DNA", *Science*) — 20 human fetal-tissue DNase-seq
samples across 8 tissues (Brain, Heart, Small Intestine, Kidney, Lung,
Muscle [arm/back/leg], Skin, Stomach), hotspot-called at FDR < 0.05 and
merged. hg19/GRCh37 coordinates. Sample list:
`data/maurano_filenames_key.tsv`.

The BEDs and pre-built hg19 FASTAs are symlinked from the old hammock
experiment dir by `prepare_data.sh`; nothing is re-downloaded.

## Design

### Modes A/B/C
- **Mode A** at precision ∈ {18, 21, 23}
- **Mode B** at precision ∈ {18, 21, 23} (mixed-stride default, subB=1.0)
- **Mode C** at p=21 with two parameter sweeps:
  - `expA ∈ {0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0}` (subB=1.0)
  - `subB ∈ {0.01, 0.05, 0.1, 0.25, 0.5, 0.8}` (expA=0)

### Mode D
- `k ∈ {8, 10, 15, 20, 25}` × `w ∈ {8, 10, 20, 30, 50, 100, 200}` (w ≥ k)
  × `p ∈ {20, 22, 23, 24}` — ~140 valid configs after filtering.
- Analysis uses `jaccard_similarity_with_ends` (boundary k-mers included);
  the boundary-less `jaccard_similarity` column is known to be noisy.

### Ground truth
`bedtools jaccard -a A -b B` over all 20×20 BED pairs. Built by
`prepare_data.sh` into `data/maurano_bedtools_ref.tsv`.

### Per-config metrics (computed by `analyze.R`)
- Pearson r, MAE, RMSE, Frobenius, max-error vs bedtools (over all 400
  pairs).
- Best Mode D config gets a scatter plot and a tissue-labelled dendrogram
  side-by-side with the bedtools reference dendrogram.

ARI/NMI/RF clustering agreement are not yet implemented here (would need
`mclust` + `phangorn`); the original experiment's value was r=0.998 and
perfect clustering at k=10, w=20–50, p=24, so the visual dendrogram
comparison plus the Pearson r heatmap should suffice for a first pass.

## Layout

```
experiments/maurano_dhs_validation/
├── README.md                 (this file)
├── prepare_data.sh           idempotent data prep (symlinks + bedtools ref)
├── run_sweep_abc.py          [optional] local A/B/C sweep driver (no Snakemake)
├── run_sweep_d.py            [optional] local Mode D sweep driver  (no Snakemake)
├── analyze.R                 plots + summaries (CairoPNG)
├── sbatch_workflow.sh        SLURM wrapper that drives the Snakemake DAG  ← preferred
├── sbatch_sweep.sh           single-job fallback (sequential, no DAG)
├── workflow/
│   ├── Snakefile             per-config rules; parallelizes across SLURM
│   ├── config.yaml           sweep grid (precisions, k/w, expA/subB)
│   └── slurm_profile/
│       └── config.yaml       Snakemake 7.x cluster profile
├── data/    → /vast/.../maurano_dhs_validation/data
│   ├── maurano.dnaseI/       20 BEDs (symlinks into old hammock dir)
│   ├── fastas/               20 hg19 FASTAs (symlinks)
│   ├── maurano_files.txt
│   ├── maurano_fastas.txt
│   ├── maurano_filenames_key.tsv
│   └── maurano_bedtools_ref.tsv
├── results/ → /vast/.../maurano_dhs_validation/results
│   ├── raw_abc/              one CSV per A/B/C config
│   ├── raw_d/                one CSV per Mode D config
│   ├── abc_summary.csv       per-config accuracy (Modes A/B/C)
│   └── mode_d_summary.csv    per-config accuracy (Mode D)
├── figures/                  Cairo PNGs (lives in-repo; small)
└── logs/    → /vast/.../maurano_dhs_validation/logs
    └── slurm/                per-job Snakemake submission logs
```

`data/`, `results/`, and `logs/` are `/vast` symlinks created by
`prepare_data.sh`. Recreate the destination tree per machine:

```bash
mkdir -p /vast/blangme2/jbonnie/hammock_claude_experiments/maurano_dhs_validation/{data,results,logs}
```

## Environment

The driver scripts call the refactored hammock installed in the
`claude-ref-comparison` conda env:
`/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock`
(both Python drivers default to this path; pass `--hammock <path>` to
override).

## How to run

The primary entry point is the Snakemake workflow in `workflow/`. Per-config
rules fan out as parallel SLURM jobs; the `run_mode_d` rule alone can have
~140 jobs in flight.

```bash
cd experiments/maurano_dhs_validation/

# 1. Dry-run shows the full DAG; sanity-check before committing cluster time.
snakemake -s workflow/Snakefile -n

# 2. Cluster execution (preferred). The wrapper sbatch runs the Snakemake
#    controller itself on a tiny node; it submits per-rule child jobs via
#    workflow/slurm_profile/.
sbatch sbatch_workflow.sh

# 3. Local execution (no SLURM): pick a thread budget.
snakemake -s workflow/Snakefile --cores 8

# 4. Partial targets:
snakemake -s workflow/Snakefile sweep_abc          # A/B/C raw CSVs only
snakemake -s workflow/Snakefile sweep_d            # Mode D raw CSVs only
snakemake -s workflow/Snakefile results/abc_summary.csv   # add analyze step
snakemake -s workflow/Snakefile -R analyze         # re-run just analyze.R

# 5. Force re-run of one config (Snakemake skips outputs that already exist):
snakemake -s workflow/Snakefile --forcerun run_mode_d \
    results/raw_d/hammock_mnmzr_p24_jaccD_k10_w20.csv
```

### Fallback: single-job sequential sbatch

`sbatch_sweep.sh` still exists as a no-Snakemake fallback — it runs the
old `run_sweep_abc.py` + `run_sweep_d.py` drivers serially in one job
slot. Use it only when you don't want Snakemake's dependency on cluster
scheduling (e.g., a dev workstation).

```bash
sbatch experiments/maurano_dhs_validation/sbatch_sweep.sh
```

## What to look for

- **Mode D Pearson heatmap** (`figures/mode_d_pearson_heatmap.png`): a band
  of r ≈ 0.99 across k=10, w∈{20,30,50}, p=24. That's the published optimum.
- **Mode D best dendrogram**: should reproduce tissue clustering exactly
  (Brain/Heart/Intestine/Kidney/Lung/Muscle/Skin/Stomach blocks).
- **Mode C expA sweep**: original hammock noted higher expA dramatically
  suppresses Jaccard — confirm at p=21.
- **Mode C subB sweep**: at p=21 we expect MAE to stay below ~1e-3 for
  subB ≥ 0.05 (consistent with the `subB_mixed_stride` experiment).

## Divergences from the original

1. **HLL-only.** The new repo doesn't ship `--minhash` or `--exact` for
   Modes A/B/C, so the old `parameter_scan/` exact/minhash comparisons are
   not part of this sweep.
2. **Mixed-stride is the new default for `--subB`.** Mode B/C subB runs
   here use mixed-stride; the original used hash-threshold. See
   `CLAUDE.md` divergence #3.
3. **Containment column is well-defined** in the new repo (formula in
   `CLAUDE.md` divergence #2). It is not currently compared against
   anything here.

## Citation
Maurano MT, Humbert R, Rynes E, et al. *Systematic localization of common
disease-associated variation in regulatory DNA.* Science.
2012;337(6099):1190–1195. doi:10.1126/science.1222794
