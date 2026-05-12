#!/bin/bash
#SBATCH --job-name=hm_subB_ms
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/subB_mixed_stride/logs/sweep_%j.out
#SBATCH --error=experiments/subB_mixed_stride/logs/sweep_%j.err

# Mode B subB sweep with --mixed-stride. Sweeps subB ∈ {1.0, 0.5, 0.25, 0.1,
# 0.05, 0.01} × 3 size classes (10k, 100k, 1M intervals) × 5 replicates;
# ground truth for accuracy is subB=1.0 on the same binary.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

mkdir -p experiments/subB_mixed_stride/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3

# Regenerate seeded BEDs only if the data dir is empty (the BEDs are
# reproducible from generate_data.py + BASE_SEED).
if [ -z "$(ls -A experiments/subB_mixed_stride/data 2>/dev/null)" ]; then
    "$PYTHON" experiments/subB_mixed_stride/generate_data.py
fi

"$PYTHON" experiments/subB_mixed_stride/run_sweep.py --threads 8

ml r/4.3.0
Rscript experiments/subB_mixed_stride/analyze.R
