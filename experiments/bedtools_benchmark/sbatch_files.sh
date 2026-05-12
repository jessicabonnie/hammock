#!/bin/bash
#SBATCH --job-name=hm_files
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/files_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/files_%j.err

# num_files sweep (the scaling-with-N axis): N ∈ {2,4,8,16,32,64,128,256},
# p=14, t=8, intervals=10000, subB ∈ {1.0, 0.25}. Uses the existing
# benchmark_cpp_vs_bedtools.py (which already supports a num_files sweep).
# N=512 omitted — bedtools at 512² pairs is ~hours and would blow the time
# limit.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --threads 8 \
    --precision 14 \
    --num-intervals 10000 \
    --num-files 2,4,8,16,32,64,128,256 \
    --subB-list 1.0,0.25,0.1 \
    --runs 3
