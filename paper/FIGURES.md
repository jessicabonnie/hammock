# Rendering the manuscript figures

Each command below regenerates one manuscript figure from committed inputs.
Run from the repository root on Rockfish after loading R:

```bash
module load libjpeg/9c r/4.3.0
export PATH=/data/apps/extern/spack_on/gcc/9.3.0/r/4.3.0-5w26bclf7ogtbbwwfx337icyt6p53t3m/bin:$PATH
```

Outputs were checked on 2026-10-28 to be pixel-identical to the figures in the
manuscript repository (`hammock_paper/figures/`); the Fig 2 PDF was compared
after rasterizing at 150 dpi.

| Fig | Manuscript file | Command |
|---|---|---|
| 1 | `roadmap_blueprint_top5_marker_pairwise_comparisons.png` | `Rscript paper/database_counts/scripts/plot_roadmap_blueprint_top_markers.R` (paths are relative to the repository root) |
| 2 | `hammock_workflow_vector.pdf` | `Rscript paper/workflow_diagram/plot_workflow_diagram.R paper/figures/hammock_workflow_vector` (the argument is the output prefix; the default, `hammock_workflow`, does not match the manuscript file name) |
| 3 | `pairwise_scaling.png` | `Rscript paper/pairwise_scaling/plot_pairwise_scaling.R` (defaults read the `docs/data/*_v0120.csv` single-build inputs) |
| 4 | `interval_accuracy.png` | `Rscript paper/interval_accuracy/plot_interval_accuracy.R` |
| 5 | `cross_reference_identity.png` | `Rscript paper/cross_reference_identity/plot_cross_reference_identity.R`, then `Rscript paper/cross_reference_identity/plot_cross_reference_parameter_plateau.R`, then `Rscript paper/cross_reference_identity/assemble_cross_reference_figure.R` |
| 6, 7 | `sequence_tissue_clustering_k20_w20_8class.png`, `sequence_tissue_distance_heatmap_k20_w20_8class.png` | `bash paper/sequence_tissue_clustering/render_rehashed_tissue_figures.sh` (do not run the two plot scripts directly: their own defaults point at the k10_w30 input) |
| 8 | `parameter_objective_tradeoff_estimators_8class.png` | `Rscript paper/parameter_objective_tradeoff/plot_parameter_objective_tradeoff_estimators.R paper/figures/parameter_objective_tradeoff_estimators_8class.png 18 8` (precision 18, eight classes; the defaults give ten classes) |
| 9 | `ten_class_partition_vacillation.png` | `Rscript paper/ten_class_partition_vacillation/plot_ten_class_partition_vacillation.R` (its TSV is a hand-made layout of the seed-17, p=18, k10_w30 partitions recorded in `experiments/maurano_dhs_rehash_validation/results/summaries/rehash_run_scores.csv`) |
