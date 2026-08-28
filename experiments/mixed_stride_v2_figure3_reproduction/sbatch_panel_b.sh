#!/bin/bash
#SBATCH --job-name=msv2_fig3b
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=experiments/mixed_stride_v2_figure3_reproduction/logs/panel_b_%j.out
#SBATCH --error=experiments/mixed_stride_v2_figure3_reproduction/logs/panel_b_%j.err

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools/2.30.0
ml parallel

EXP=experiments/mixed_stride_v2_figure3_reproduction
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
BINARY=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp
mkdir -p "$EXP/logs" "$EXP/results/panel_b"

echo "node: $(hostname)  cores: $(nproc)"
echo "commit: $(git rev-parse HEAD)"
"$BINARY" --version

"$PYTHON" experiments/subB_mixed_stride/run_sweep.py \
  --binary "$BINARY" --corpus maurano --methods mixed-stride \
  --subB-values 1.0 0.3 0.1 0.03 0.01 0.008 0.005 0.003 0.001 \
  --reps 5 --precision 18 --threads 8 --metrics \
  --output "$EXP/results/panel_b/hammock.csv"

"$PYTHON" experiments/subB_mixed_stride/run_bedtools_maurano.py \
  --threads 8 --reps 5 \
  --output "$EXP/results/panel_b/bedtools.csv"
