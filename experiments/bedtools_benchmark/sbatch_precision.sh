#!/bin/bash
#SBATCH --job-name=hm_precision
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/precision_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/precision_%j.err

# Precision sweep: p ∈ {10,12,14,16,18}, t=8, files=64, intervals=10000,
# subB ∈ {1.0, 0.25}. Bedtools is precision/subB-independent → run once per
# data realization for reference.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

"$PYTHON" experiments/bedtools_benchmark/sweep.py \
    --axis precision \
    --precisions 10,12,14,16,18 \
    --threads 8 \
    --num-files 64 \
    --num-intervals 10000 \
    --subB-list 1.0,0.25,0.1 \
    --runs 3
