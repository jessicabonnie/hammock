#!/bin/bash
#SBATCH --job-name=maurano_d_k5_w20_p24
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=experiments/maurano_dhs_validation/logs/slurm/d_k5_w20_p24_%j.out
#SBATCH --error=experiments/maurano_dhs_validation/logs/slurm/d_k5_w20_p24_%j.err

# One-off: run Mode D at k=5, w=20, p=24 on the 20 Maurano FASTAs.
# Companion to sbatch_k5_w5_p24.sh — tests whether w > k breaks the
# saturation that k=5, w=5 hit (all jaccards = 1.0). Renders a
# tissue-coloured dendrogram and prints the jaccard distribution.
# Does not touch the Snakemake DAG.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude/experiments/maurano_dhs_validation

HAMMOCK=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock
CSV=results/raw_d/hammock_mnmzr_p24_jaccD_k5_w20.csv
PNG=figures/mode_d_dendrogram_k5_w20_p24.png

mkdir -p results/raw_d figures logs/d logs/slurm

echo "[$(date)] hammock Mode D k=5 w=20 p=24"
"$HAMMOCK" data/maurano_fastas.txt data/maurano_fastas.txt \
    --mode D --precision 24 -k 5 -w 20 \
    --threads 8 \
    --outprefix results/raw_d/hammock \
    > logs/d/p24_k5_w20.log 2>&1

if [[ ! -s "$CSV" ]]; then
    echo "ERROR: expected output not produced: $CSV" >&2
    exit 1
fi

echo "[$(date)] off-diagonal jaccard distribution (column 9 = jaccard_similarity):"
awk -F, 'NR>1 && $1!=$2 {print $9}' "$CSV" | sort -n | uniq -c | head

echo "[$(date)] rendering dendrogram"
ml gcc/9.3.0 r/4.3.0 libjpeg/9c
Rscript render_dendrogram.R "$CSV" "$PNG"

echo "[$(date)] done"
