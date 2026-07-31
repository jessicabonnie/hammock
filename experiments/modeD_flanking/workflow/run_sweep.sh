#!/usr/bin/env bash
#SBATCH --job-name=part2_sweep
#SBATCH --time=24:00:00
#SBATCH --partition=shared
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/modeD_flanking/logs/part2_sweep.%j.log

set -eu
EXPT=/home/jbonnie1/interval_sketch/hammock_claude/experiments/modeD_flanking
cd "$EXPT"
echo "[sweep] start $(date)"
python3 run_sweep_synthetic.py --threads 8
echo "[sweep] done  $(date)"
