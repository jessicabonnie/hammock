# maurano_dhs_validation — Modes A/B/C/D vs bedtools on fetal-tissue DHS

How accurately do Hammock's four modes (with HLL backing) recover bedtools
pairwise Jaccard on the Maurano et al. 2012 fetal-tissue DNase-seq corpus,
and what parameters are optimal?

This is the `hammock_claude` recreation of `dnase1-hypersensitivity` from the
original `hammock` repo. The headline result we expected to reproduce: Mode D
at k=10, w=20–50, p=24 hits Pearson r ≈ 0.998 against bedtools with perfect
clustering reproduction (ARI=NMI=1, RF=0) on the 20-sample x 20-sample
matrix.

**What actually happened is in [RESULTS.md](RESULTS.md)** — the Pearson
optimum landed at high k + high w (r = 0.9996) and the clustering optimum at
a different cell (k=10, w=30, ARI = 0.910), i.e. the two optima are not
coincident. Plot-generation notes for the two paper figures are in
[PLOT_GENERATION.md](PLOT_GENERATION.md).

## Corpus

Maurano et al. 2012 ("Systematic localization of common disease-associated
variation in regulatory DNA", *Science*) — 20 human fetal-tissue DNase-seq
samples across 8 tissues (Brain, Heart, Small Intestine, Kidney, Lung,
Muscle [arm/back/leg], Skin, Stomach), hotspot-called at FDR < 0.05 and
merged. The clustering metrics score against **10** labels, since the three
muscle sub-tissues are kept distinct. hg19/GRCh37 coordinates. Sample list:
`data/maurano_filenames_key.tsv`.

The BEDs and pre-built hg19 FASTAs are symlinked from the old hammock
experiment dir by `prepare_data.sh`; nothing is re-downloaded.

## Design

The authoritative grid is `workflow/config.yaml`; the summary below tracks it.

### Modes A/B/C
- **Mode A** at precision ∈ {18, 21, 23}
- **Mode B** at precision ∈ {18, 21, 23} (mixed-stride default, subB=1.0)
- **Mode C** at p=21 with two parameter sweeps:
  - `expA ∈ {0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0}`
    (subB=1.0)
  - `subB ∈ {0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.80}`
    (expA=0)
- 29 A/B/C configs in total (3 + 3 + 23).

### Mode D
- Base grid: `k ∈ {8, 10, 15, 20, 25}` × `w ∈ {8, 10, 20, 30, 50, 100, 200}`
  (w ≥ k) × `p ∈ {20, 22, 23, 24}`.
- Extended grid (added after run 1): `k ∈ {8, 10}` × `w ∈ {…, 300, 500}` ×
  `p ∈ {10, 12, 14, 16, 18, 20, 22, 23, 24}`, unioned with the base grid.
- Fill-in jobs outside `config.yaml`: `sbatch_fill_highk_w.sh`
  (k ∈ {15,20,25} × w ∈ {300,500} × p ∈ {20,22,23,24}) and the two one-off
  k=5 scripts. `results/mode_d_summary.csv` currently holds **235** configs.
- Analysis automatically summarizes the seven current minimizer metrics:
  `jaccard_similarity`, `jaccard_similarity_ie`, `containment_AB/BA`, and
  `cosketch_geom/arith/max`. Published figures remain explicitly filtered to
  `jaccard_similarity`.

> **`jaccard_similarity_with_ends` was removed in v0.6.0** (`CLAUDE.md`
> divergence #8). The archived CSVs under `results/raw_d/` still carry it and
> `analyze.R` intentionally ignores it; a fresh Mode D run will not emit it.
> The earlier
> claim here that the boundary-less column "is known to be noisy" was an
> artifact of the `add_string` bug fixed 2026-05-14; post-fix it is the better
> of the two on every metric.
>
> **Parse Mode D CSVs by column name, not position.** The Mode D header has
> changed four times since this sweep ran: containment/cosketch block
> (2026-05-14), `jaccard_similarity_ie` (v0.5.0), `with_ends` family dropped
> (v0.6.0), trailing `ref1`/`ref2` for bed2fasta runs.
>
> **`jaccard_similarity` is register-equality, not set Jaccard.** A high
> Pearson r against bedtools means the two covary affinely, not that the
> values agree; see `docs/jaccard-definitional-gap.md` and the caveat box in
> `RESULTS.md`. `jaccard_similarity_ie` (v0.5.0) is the calibrated column and
> is now scored alongside it in `mode_d_summary.csv`.

### Ground truth
`bedtools jaccard -a A -b B` over all 20×20 BED pairs. Built by
`prepare_data.sh` into `data/maurano_bedtools_ref.tsv`.

### Per-config metrics (computed by `analyze.R`)
- Pearson r, Spearman ρ, MAE, RMSE, Frobenius, max-error vs bedtools **and**
  vs Mode B (over all 400 pairs). Both references land in
  `results/mode_d_summary.csv` under the `reference` column.
- ARI and NMI against the 10 fetal-tissue labels (implemented since the first
  draft of this file; `mode_d_summary.csv` has `ari`/`nmi` columns and
  `figures/mode_d_clustering_{ari,nmi}.png`).
- Best Mode D config gets a scatter plot and a tissue-labelled dendrogram
  side-by-side with the bedtools reference dendrogram, plus a cluster
  contingency table (`results/best_cluster_{assignment,contingency}.csv`).

Robinson–Foulds tree distance is still not computed (would need `phangorn`);
the dendrogram comparison is visual.

## Layout

```
experiments/maurano_dhs_validation/
├── README.md                 (this file — question, design, how to run)
├── RESULTS.md                findings + figures
├── PLOT_GENERATION.md        spec for the two paper metric figures (executed)
├── prepare_data.sh           idempotent data prep (symlinks + bedtools ref)
├── run_sweep_abc.py          [optional] local A/B/C sweep driver (no Snakemake)
├── run_sweep_d.py            [optional] local Mode D sweep driver  (no Snakemake)
├── analyze.R                 plots + summaries (CairoPNG)
├── mode_c_interpolation.R    Mode C ↔ Mode A/B interpolation figures
├── render_dendrogram.R       standalone dendrogram renderer: <csv> <png>
├── scripts/
│   └── make_metric_plots.R   the three PLOT_GENERATION.md figures
├── sbatch_workflow.sh        SLURM wrapper that drives the Snakemake DAG  ← preferred
├── sbatch_sweep.sh           single-job fallback (sequential, no DAG)
├── sbatch_fill_highk_w.sh    fill-in: k∈{15,20,25} × w∈{300,500} × p∈{20,22,23,24}
├── sbatch_k5_w5_p24.sh       one-off k=5 w=5  p=24 (saturates: all jaccards 1.0)
├── sbatch_k5_w20_p24.sh      one-off k=5 w=20 p=24 (companion to the above)
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
│   ├── raw_d_buggy_pre_fix/  pre-2026-05-14 Mode D CSVs; do not analyse
│   ├── abc_summary.csv       per-config accuracy (Modes A/B/C)
│   ├── mode_d_summary.csv    per-config accuracy (Mode D)
│   ├── best_cluster_assignment.csv / best_cluster_contingency.csv
│   └── sweep_d_*.csv         Mode D sweep index (also cited in CLAUDE.md as
│                             the cpu_s/wall_s ≈ 0.8 threading evidence)
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
override). R steps need `ml r/4.3.0` and render via `Cairo::CairoPNG()` —
system R 4.3.0 has no native PNG device.

> **Mode D wall times from this experiment are inflated (noted 2026-08-06).**
> Every Mode D invocation here ran at `--threads 8`
> (`config.yaml: threads_per_run: 8`, and the same in the sbatch scripts).
> v0.6.1 changed the Mode D default to `--threads 1` because the thread pool
> was a GIL convoy — measured 2–4.5× *slower* than single-threaded. The
> accuracy numbers are unaffected (the CSVs are byte-identical either way),
> but do not compare any timing from this experiment against a post-0.6.1
> run. See `docs/seed-mode-d-threading.md`.

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

The two Python drivers take: `--threads`, `--hammock`, `--dry-run`, plus
`--ks/--ws/--ps` (Mode D) or `--subset` (A/B/C). Neither takes a grid file —
the Snakemake path is the one that reads `config.yaml`, so a driver run and a
Snakemake run can disagree about the grid unless you pass it explicitly.

### Analysis / figures only

```bash
ml r/4.3.0
Rscript analyze.R                      # summaries + most figures; no arguments
Rscript mode_c_interpolation.R         # the two mode_c_*_interpolation_agg PNGs
Rscript scripts/make_metric_plots.R    # the three PLOT_GENERATION.md figures
Rscript render_dendrogram.R <csv> <png>
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
