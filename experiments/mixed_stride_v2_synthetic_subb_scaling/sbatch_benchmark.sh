#!/bin/bash
#SBATCH --job-name=msv2_synsubb
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=experiments/mixed_stride_v2_synthetic_subb_scaling/logs/benchmark_%j.out
#SBATCH --error=experiments/mixed_stride_v2_synthetic_subb_scaling/logs/benchmark_%j.err

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools/2.30.0
ml parallel

EXP=experiments/mixed_stride_v2_synthetic_subb_scaling
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
BINARY=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp
mkdir -p "$EXP/logs" "$EXP/results/raw" "$EXP/figures"

echo "node: $(hostname)  cores: $(nproc)"
echo "commit: $(git rev-parse HEAD)"
"$BINARY" --version

"$PYTHON" "$EXP/run_benchmark.py" \
  --binary "$BINARY" \
  --output-dir "$EXP/results/raw" \
  --num-files 2,4,8,16,32,64,128,256,512 \
  --small-runs 20 --large-runs 3 \
  --num-intervals 10000 --precision 18 --threads 16 \
  --corpus-seed 20260810
