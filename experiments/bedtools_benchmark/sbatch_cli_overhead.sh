#!/bin/bash
#SBATCH --job-name=cli_overhead
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=experiments/bedtools_benchmark/logs/cli_overhead_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/cli_overhead_%j.err

# Quantifies the hammock CLI vs hammock-cpp overhead gap for the
# methods/limitations text -- see measure_cli_overhead.py's docstring.
# Small and quick: order-of-magnitude characterization, not a headline figure.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
mkdir -p experiments/bedtools_benchmark/logs

echo "node: $(hostname)  cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)  cores: $(nproc)"

/home/jbonnie1/.conda/envs/hammock/bin/python3 experiments/bedtools_benchmark/measure_cli_overhead.py \
    --num-files 2,8,32,128 \
    --num-intervals 10000 \
    --precision 18 \
    --threads 16 \
    --runs 5 \
    --corpus-seed 20260811
