#!/bin/bash
#SBATCH --job-name=sbr_threads
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/threads_%j.out
#SBATCH --error=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/threads_%j.err

# Thread scaling with the single 0.12.0 build.
# Reproduces experiments/bedtools_benchmark/sbatch_fig3_threading.sh, archived
# 0.7.0 run job 29670792 (c531) -> docs/data/sweep_threads_p18.csv:
#   two synthetic 64-file collections x 10,000 intervals, p = 18, subB 1.0,
#   threads/processes 1,2,4,8,16,32,48, three replicates, bedtools FIRST at
#   every setting (--bedtools-first: the current sweep.py rotates arm order,
#   which the archived run did not), hammock with --register-equality output.
# The corpus is unseeded, as in the archived run (fresh draw per replicate), so
# inputs are statistically but not byte-identical to job 29670792's.

set -euo pipefail
source /home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/common.sh
require_cpu 6248R
require_sha "$SINGLE_BINARY" "$SINGLE_SHA"
cd "$EXP"

ml bedtools/2.30.0
ml parallel
export HAMMOCK_CPP_BIN=$SINGLE_BINARY

OUT=$EXP/results/threads
FIG=$EXP/figures/threads
mkdir -p "$OUT" "$FIG"
META=$OUT/metadata_${SLURM_JOB_ID}.txt
write_metadata "$META" "$SINGLE_BINARY"

"$PY_HAMMOCK" harness/sweep.py \
    --binary "$SINGLE_BINARY" \
    --axis threads \
    --thread-list 1,2,4,8,16,32,48 \
    --precision 18 \
    --num-files 64 \
    --num-intervals 10000 \
    --subB-list 1.0 \
    --runs 3 \
    --bedtools-first \
    --output-dir "$OUT" --figures-dir "$FIG"

finish_metadata "$META"
