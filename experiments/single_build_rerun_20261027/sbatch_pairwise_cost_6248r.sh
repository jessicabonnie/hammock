#!/bin/bash
#SBATCH --job-name=sbr_pcost_6248r
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/pairwise_cost_6248r_%j.out
#SBATCH --error=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/pairwise_cost_6248r_%j.err

# Per-pair estimator cost with the SINGLE 6248R build, on its own node type.
# Same protocol as sbatch_pairwise_cost_6448y.sh (64-file collections, 16
# threads, p = 12..24 even, five replicates, --register-equality vs --metrics,
# corpus seed 20261027), but on an exclusive 48-core 6248R node instead of 16
# reserved cores on a 6448Y node. This is the variant that keeps every
# manuscript benchmark on one binary; the 6448Y run is the design-matched one.

set -euo pipefail
source /home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/common.sh
require_cpu 6248R
require_sha "$SINGLE_BINARY" "$SINGLE_SHA"
cd "$EXP"
export HAMMOCK_CPP_BIN=$SINGLE_BINARY

OUT=$EXP/results/pairwise_cost_6248r
mkdir -p "$OUT"
META=$OUT/metadata_${SLURM_JOB_ID}.txt
write_metadata "$META" "$SINGLE_BINARY"

"$PY_REF" harness/pairwise_cost_by_precision.py \
    --binary "$SINGLE_BINARY" \
    --precisions 12,14,16,18,20,22,24 \
    --num-files 64 \
    --num-intervals 10000 \
    --threads 16 \
    --runs 5 \
    --corpus-seed 20261027 \
    --output-dir "$OUT"

finish_metadata "$META"
