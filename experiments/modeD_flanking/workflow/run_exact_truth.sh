#!/usr/bin/env bash
#SBATCH --job-name=part2_exact_kmer
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/modeD_flanking/logs/part2_exact_kmer.%j.log

set -eu
EXPT=/home/jbonnie1/interval_sketch/hammock_claude/experiments/modeD_flanking
cd "$EXPT"
echo "[exact_kmer] start $(date)"
python3 exact_kmer_jaccard.py --ks 8 10 15 20
echo "[exact_kmer] done  $(date)"
