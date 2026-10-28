#!/bin/bash
#SBATCH --job-name=sbr_pcost_6448y
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/pairwise_cost_6448y_%j.out
#SBATCH --error=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/pairwise_cost_6448y_%j.err

# Per-pair estimator cost, design-matched to the archived 0.7.0 run
# (experiments/bedtools_benchmark/sbatch_pairwise_cost.sh, job 29620391, sr09,
# Xeon Gold 6448Y, 16 reserved cores ->
# results/pairwise_cost_by_precision_20260808_183458.csv):
#   two synthetic 64-file collections x 10,000 intervals, 16 threads,
#   p = 12..24 even, five replicates, arms --register-equality ("no_metrics")
#   vs --metrics ("metrics": full block incl. inclusion-exclusion, the same
#   output shape the 0.7.0 default produced), arm order alternating by replicate.
# NOTE (2026-10-27): build job 32614873 completed the compile but exited on the
# ymm>100 gate (ymm=20, zmm=78): GCC 14 -march=native on Sapphire Rapids emits
# 512-bit AVX-512 code rather than the 256-bit code of the 6248R build. The
# binary is vectorized but at a different width, so this run is resubmitted
# without the build dependency and is NOT the same codegen as the single build.
# Binary: the db1411a/GCC 14.3 NATIVE 6448Y build from sbatch_build_6448y.sh
# (a second build of the same commit; the 6248R binary is measured on its own
# node type by sbatch_pairwise_cost_6248r.sh).
# Corpus seeded (20261027) so this run and the 6248R run sketch identical
# inputs; the archived run was unseeded.

set -euo pipefail
source /home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/common.sh
require_cpu 6448Y
cd "$EXP"

BIN=$EXP/build_sr6448y_gcc14/hammock-cpp
BUILD_SHA=$(sed -n 's/^binary_sha256: //p' "$EXP/build_metadata_sr6448y.txt")
[[ -n "$BUILD_SHA" ]] || { echo "no build metadata" >&2; exit 1; }
require_sha "$BIN" "$BUILD_SHA"
export HAMMOCK_CPP_BIN=$BIN

OUT=$EXP/results/pairwise_cost_6448y
mkdir -p "$OUT"
META=$OUT/metadata_${SLURM_JOB_ID}.txt
write_metadata "$META" "$BIN"

# Smoke run through the real kernels (as in the archived wrapper).
_probe=$(mktemp -d)
trap 'rm -rf "$_probe"' EXIT
printf 'chr1\t100\t200\nchr1\t500\t900\n' > "$_probe/a.bed"
printf 'chr1\t150\t250\nchr2\t10\t99\n'   > "$_probe/b.bed"
printf '%s\n' "$_probe/a.bed" "$_probe/b.bed" > "$_probe/list.txt"
"$BIN" "$_probe/list.txt" "$_probe/list.txt" --mode B -p 14 -t 2 -o "$_probe/out" >/dev/null

"$PY_REF" harness/pairwise_cost_by_precision.py \
    --binary "$BIN" \
    --precisions 12,14,16,18,20,22,24 \
    --num-files 64 \
    --num-intervals 10000 \
    --threads 16 \
    --runs 5 \
    --corpus-seed 20261027 \
    --output-dir "$OUT"

finish_metadata "$META"
