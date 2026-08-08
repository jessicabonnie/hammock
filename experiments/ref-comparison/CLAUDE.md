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

- `docs/exp_a_results.md` — **the results record; read before touching anything here**
- `docs/experiment_design.md` — full design rationale, accession tables, success criteria
- `docs/paper_outline.md` — current paper outline (Exp A only)
- `config/config.yaml` — sample metadata, nf-core output paths, sketch parameters
- `config/nfcore_samplesheet_human.csv` — nf-core input: 3 ChIP-seq samples + matched controls × 3 references
- `workflow/Snakefile` — downstream pipeline (peaks → FASTA → hammock → plots)
- `workflow/slurm_profile/config.yaml` — SLURM cluster profile for Snakemake
- `scripts/run_nfcore.sh` — launches nf-core/chipseq for GRCh38, GRCh37, CHM13
- `scripts/exp_a_{validate_plot,sweep_summary,dendrogram,metric_comparison}.R` — the four analysis passes

## Working here — standing caveats

- **There is no single lead cell; use the k ≥ 15 plateau.** "Lead cell is
  k=15, w=15" was established on the deleted `_with_ends` column. Rescored on
  `jaccard_similarity` across all 20 cells (2026-08-08,
  `estimator_ie_crossref.py`), k15_w15 ranks **fifth** at Δ 0.483 broad, behind
  k20_w30 (0.540), k20_w20 (0.533), k15_w30 (0.494) and k15_w20 (0.488). The
  2026-08-06 recomputation in `docs/exp_a_results.md` did not catch this because
  it rescored only two cells. The top five are within 0.06 of each other, so
  none is meaningfully best — quote the plateau, not an argmax. Pre-fix notes
  calling k=15/20 a "saturated low" regime are still wrong; that part holds.
- **Archived CSVs carry `jaccard_similarity_with_ends`; hammock no longer
  emits it** (v0.6.0, `CLAUDE.md` divergence #8). The stats tables in
  `docs/exp_a_results.md` were computed on it. A re-run will not reproduce
  them column-for-column — see the 2026-08-06 recomputation in that file.
  `config/config.yaml`'s `primary_sim_col` is already on `jaccard_similarity`,
  but note `workflow/Snakefile` hardcodes the column rather than reading it.
- **Parse these CSVs by column name, never by field index.** The schema
  differs between the archived and current versions.
- **`jaccard_similarity` is register-equality, not set Jaccard** — it has a
  chance-agreement floor and is not rank-faithful across pairs of differing
  set size. The calibrated `jaccard_similarity_ie` is exactly recoverable
  from the archived `containment_AB`/`containment_BA` columns via
  `J = 1/(1/C_AB + 1/C_BA − 1)`.

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
