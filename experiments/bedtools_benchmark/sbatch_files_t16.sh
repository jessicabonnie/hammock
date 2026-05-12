#!/bin/bash
#SBATCH --job-name=hm_files_t16
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/files_t16_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/files_t16_%j.err

# Apples-to-apples replication of the orig hammock Oct 2025 benchmark:
# t=16, p=14, 10k intervals/file, N up to 512, subB ∈ {1.0, 0.25}. The orig
# run showed hammock winning at N≥128 and 5.83× at N=512; our t=8 files run
# capped out at N=256 without reaching the wide-N regime where hammock's
# pairwise-compare cost stops dominating. This re-runs the same protocol at
# the same thread count and adds a subB=0.25 hammock line for the
# subsampling comparison.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --threads 16 \
    --precision 14 \
    --num-intervals 10000 \
    --num-files 2,4,8,16,32,64,128,256,512 \
    --subB-list 1.0,0.25 \
    --runs 3
