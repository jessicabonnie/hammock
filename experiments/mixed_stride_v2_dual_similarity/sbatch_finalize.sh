#!/bin/bash
#SBATCH --job-name=msv2_replot
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=experiments/mixed_stride_v2_dual_similarity/logs/finalize_%j.out
#SBATCH --error=experiments/mixed_stride_v2_dual_similarity/logs/finalize_%j.err

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/mixed_stride_v2_dual_similarity
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
IE_EXP=experiments/mixed_stride_v2_figure3_reproduction
mkdir -p "$EXP/logs" "$EXP/results/prepared" "$EXP/figures"

"$PYTHON" "$EXP/prepare_inputs.py" \
  --register-equality "$EXP/results/register_equality.csv" \
  --inclusion-exclusion "$IE_EXP/results/prepared/panel_b_summary.csv" \
  --truth docs/data/maurano_bedtools_ref.tsv \
  --output-dir "$EXP/results/prepared"

Rscript "$EXP/plot_dual_similarity.R" \
  "$EXP/results/prepared/dual_similarity_plot.csv" \
  "$IE_EXP/results/panel_b/bedtools.csv" \
  "$EXP/figures/register_equality_vs_ie.png"
