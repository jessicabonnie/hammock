#!/bin/bash
#SBATCH --job-name=msv2_synplot
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=experiments/mixed_stride_v2_synthetic_subb_scaling/logs/finalize_%j.out
#SBATCH --error=experiments/mixed_stride_v2_synthetic_subb_scaling/logs/finalize_%j.err

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/mixed_stride_v2_synthetic_subb_scaling
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
mkdir -p "$EXP/logs" "$EXP/results/prepared" "$EXP/figures"

"$PYTHON" "$EXP/prepare_results.py" \
  --raw-dir "$EXP/results/raw" \
  --output-dir "$EXP/results/prepared"

Rscript "$EXP/plot_results.R" \
  "$EXP/results/prepared/summary.csv" \
  "$EXP/figures/synthetic_subb_scaling.png"
