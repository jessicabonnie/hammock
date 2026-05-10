#!/bin/bash
#SBATCH --job-name=hm_intervals
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/intervals_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/intervals_%j.err

# Intervals/file sweep (the asymptotic axis): N ∈ {1k, 10k, 100k, 1M},
# p=14, t=8, files=16. num_files reduced to keep bedtools tractable at 1M.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3

"$PYTHON" experiments/bedtools_benchmark/sweep.py \
    --axis intervals \
    --intervals-list 1000,10000,100000,1000000 \
    --precision 14 \
    --threads 8 \
    --num-files 16 \
    --runs 3
