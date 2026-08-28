#!/bin/bash
#SBATCH --job-name=msv2_fig3plot
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=experiments/mixed_stride_v2_figure3_reproduction/logs/finalize_%j.out
#SBATCH --error=experiments/mixed_stride_v2_figure3_reproduction/logs/finalize_%j.err

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/mixed_stride_v2_figure3_reproduction
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
mkdir -p "$EXP/logs" "$EXP/results/prepared" "$EXP/figures"

"$PYTHON" "$EXP/prepare_inputs.py" \
  --experiment-dir "$EXP" --truth docs/data/maurano_bedtools_ref.tsv

Rscript paper/pairwise_scaling/plot_pairwise_scaling.R \
  "$EXP/figures/pairwise_scaling_reproduced.png" \
  docs/data/cpp_vs_bedtools_t16_p18.csv \
  "$EXP/results/prepared/panel_b_figure3_summary.csv" \
  "$EXP/results/panel_b/bedtools.csv"

Rscript paper/pairwise_scaling/plot_pairwise_scaling.R \
  "$EXP/figures/pairwise_scaling_reproduced_expanded.png" \
  docs/data/cpp_vs_bedtools_t16_p18.csv \
  "$EXP/results/prepared/panel_b_summary.csv" \
  "$EXP/results/panel_b/bedtools.csv"
