#!/bin/bash
#SBATCH --job-name=msv2_re
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=experiments/mixed_stride_v2_dual_similarity/logs/register_equality_%j.out
#SBATCH --error=experiments/mixed_stride_v2_dual_similarity/logs/register_equality_%j.err

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml parallel

EXP=experiments/mixed_stride_v2_dual_similarity
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
BINARY=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp
mkdir -p "$EXP/logs" "$EXP/results"

echo "node: $(hostname)  cores: $(nproc)"
echo "commit: $(git rev-parse HEAD)"
"$BINARY" --version

# Omitting run_sweep.py's --metrics flag explicitly selects hammock-cpp's
# cheap --register-equality output shape.
"$PYTHON" experiments/subB_mixed_stride/run_sweep.py \
  --binary "$BINARY" --corpus maurano --methods mixed-stride \
  --subB-values 1.0 0.3 0.1 0.03 0.01 0.008 0.005 0.003 0.001 \
  --reps 5 --precision 18 --threads 8 \
  --output "$EXP/results/register_equality.csv"
