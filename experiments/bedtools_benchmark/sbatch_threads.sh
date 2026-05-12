#!/bin/bash
#SBATCH --job-name=hm_threads
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/threads_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/threads_%j.err

# Thread sweep: t ∈ {1,2,4,8,16}, p=14, files=64, intervals=10000,
# subB ∈ {1.0, 0.25}. t=1 is the long pole — bedtools at 64 files = 4096
# sequential pairs.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

"$PYTHON" experiments/bedtools_benchmark/sweep.py \
    --axis threads \
    --thread-list 1,2,4,8,16 \
    --precision 14 \
    --num-files 64 \
    --num-intervals 10000 \
    --subB-list 1.0,0.25 \
    --runs 3
