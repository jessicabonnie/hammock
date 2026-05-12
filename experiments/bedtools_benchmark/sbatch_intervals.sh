#!/bin/bash
#SBATCH --job-name=hm_intervals
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/intervals_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/intervals_%j.err

# Intervals/file sweep (the asymptotic axis): N ∈ {1k, 10k, 100k, 1M},
# p=14, t=8, files=16, subB ∈ {1.0, 0.25}. num_files reduced to keep
# bedtools tractable at 1M. Time budget bumped to 6h: at 1M intervals
# the hammock @ subB=1.0 run alone is ~5 min; the additional subB=0.25
# run roughly doubles hammock walltime per realization.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

"$PYTHON" experiments/bedtools_benchmark/sweep.py \
    --axis intervals \
    --intervals-list 1000,10000,100000,1000000 \
    --precision 14 \
    --threads 8 \
    --num-files 16 \
    --subB-list 1.0,0.25,0.1 \
    --runs 3
