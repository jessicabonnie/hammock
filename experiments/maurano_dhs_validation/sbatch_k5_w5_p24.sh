#!/bin/bash
#SBATCH --job-name=maurano_d_k5_w5_p24
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=experiments/maurano_dhs_validation/logs/slurm/d_k5_w5_p24_%j.out
#SBATCH --error=experiments/maurano_dhs_validation/logs/slurm/d_k5_w5_p24_%j.err

# One-off: run Mode D at k=5, w=5, p=24 on the 20 Maurano FASTAs, then
# render a tissue-coloured dendrogram. Does not touch the Snakemake DAG —
# the resulting CSV lands in results/raw_d/ where a future snakemake run
# will see it and skip.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude/experiments/maurano_dhs_validation

HAMMOCK=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock
# "_full" tag matches --metrics below (python/hammock/outprefix.py now
# always tags Python-CLI output _ie/_re/_full; docs/seed-metrics-column-restructure.md)
# -- also needed so render_dendrogram.R's default jaccard_similarity read
# below has a column to find, and so this lands under the same _full tag
# the Snakefile's D_CSVS_BASE/D_CSVS_EXT patterns now require (see the
# "future snakemake run will see it and skip" note above).
CSV=results/raw_d/hammock_mnmzr_p24_jaccD_k5_w5_full.csv
PNG=figures/mode_d_dendrogram_k5_w5_p24.png

mkdir -p results/raw_d figures logs/d logs/slurm

echo "[$(date)] hammock Mode D k=5 w=5 p=24"
"$HAMMOCK" data/maurano_fastas.txt data/maurano_fastas.txt \
    --mode D --precision 24 -k 5 -w 5 \
    --threads 8 \
    --outprefix results/raw_d/hammock \
    --metrics \
    > logs/d/p24_k5_w5.log 2>&1

if [[ ! -s "$CSV" ]]; then
    echo "ERROR: expected output not produced: $CSV" >&2
    exit 1
fi

echo "[$(date)] rendering dendrogram"
ml gcc/9.3.0 r/4.3.0 libjpeg/9c
Rscript render_dendrogram.R "$CSV" "$PNG"

echo "[$(date)] done"
