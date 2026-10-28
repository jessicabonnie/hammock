#!/bin/bash
#SBATCH --job-name=sbr_precision
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/precision_%j.out
#SBATCH --error=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/precision_%j.err

# Precision sweep (Figure 4B) and Maurano accuracy runs (Figure 4A p=18, and
# the in-text p=21 run) with the single 0.12.0 build.
#
# Part 1 reproduces the t=16 half of experiments/bedtools_benchmark/
# sbatch_fig3_panelB.sh, archived 0.7.0 run job 29670793 (c432) ->
# docs/data/sweep_precision_maurano_p18_t16.csv (+ results/
# sweep_precision_20260809_230812_pairs.csv):
#   20-file Maurano subset as both operands (400 timed ordered pairs),
#   p = 12..24 even, 16 threads/processes, three replicates, bedtools once per
#   replicate before the hammock precision series, hammock wall times with
#   --register-equality output plus a separate IE pass per precision.
#   The archived job also ran a t=8 series afterwards; it is not cited in the
#   manuscript and is not repeated here.
#   Known difference: the IE pass now runs the bare default (IE column only)
#   instead of --metrics (Part 9 harness change, as in confirming rerun 29769995).
#
# Part 2 (the archived p=18/p=21 estimates came from hammock 0.5.0 via the
# Python CLI, commit 6007450 in hammock_claude): one deterministic all-pairs run per precision, p = 18 and 21, with
# --metrics so the output carries both reg_eq_similarity and
# jaccard_similarity_ie in the same shape as docs/data/hammock_hll_p{18,21}_jaccB_full.csv.

set -euo pipefail
source /home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/common.sh
require_cpu 6248R
require_sha "$SINGLE_BINARY" "$SINGLE_SHA"
cd "$EXP"

ml bedtools/2.30.0
ml parallel
export HAMMOCK_CPP_BIN=$SINGLE_BINARY

OUT=$EXP/results/precision
FIG=$EXP/figures/precision
ACC=$EXP/results/accuracy
mkdir -p "$OUT" "$FIG" "$ACC"
META=$OUT/metadata_${SLURM_JOB_ID}.txt
write_metadata "$META" "$SINGLE_BINARY"

echo; echo "================ precision sweep, threads = 16 ================"
"$PY_HAMMOCK" harness/sweep.py \
    --binary "$SINGLE_BINARY" \
    --axis precision \
    --corpus maurano \
    --precisions 12,14,16,18,20,22,24 \
    --threads 16 \
    --subB-list 1.0 \
    --runs 3 \
    --output-dir "$OUT" --figures-dir "$FIG"
echo "sweep_end: $(date -Is)" | tee -a "$META"

echo; echo "================ accuracy runs, p = 18, 21 ================"
LIST=$ACC/maurano_list_${SLURM_JOB_ID}.txt
[[ ! -e "$LIST" ]] || { echo "refusing to overwrite $LIST" >&2; exit 1; }
ls /home/jbonnie1/interval_sketch/hammock_claude/experiments/subB_mixed_stride/data/maurano_sorted/*.bed > "$LIST"
[[ $(wc -l < "$LIST") -eq 20 ]] || { echo "expected 20 Maurano BEDs" >&2; exit 1; }
for P in 18 21; do
  PREFIX=$ACC/maurano_job${SLURM_JOB_ID}
  "$SINGLE_BINARY" "$LIST" "$LIST" --mode B -p "$P" --metrics -t 16 \
      -o "$PREFIX" --verbose 2>&1 | tee "$ACC/maurano_job${SLURM_JOB_ID}_p${P}.log"
  test -s "${PREFIX}_hll_p${P}_jaccB_full.csv"
done

finish_metadata "$META"
