#!/bin/bash
#SBATCH --job-name=hm_files_big
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/files_big_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/files_big_%j.err

# Extension of the files_t16 sweep to N=1024 — pushes hammock's win
# regime further into the asymptotic so the crossover plot has more
# post-crossover area. Includes N=512 as a same-machine-state reference
# point to verify continuity with the earlier sweep.

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
    --num-files 512,1024 \
    --runs 3
