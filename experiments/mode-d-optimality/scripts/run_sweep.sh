#!/bin/bash
#SBATCH --job-name=mode-d-sweep
#SBATCH --partition=bigmem
#SBATCH --account=blangme2_bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/vast/blangme2/jbonnie/hammock/mode-d-optimality/results/slurm-%j.out
#SBATCH --error=/vast/blangme2/jbonnie/hammock/mode-d-optimality/results/slurm-%j.err

set -euo pipefail

EXPERIMENT_DIR="/home/jbonnie1/interval_sketch/hammock/experiments/mode-d-optimality"

ml anaconda
conda activate hammock

cd "$EXPERIMENT_DIR"

python scripts/sweep_kw.py "$@"
