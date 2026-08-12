#!/bin/bash
#SBATCH --job-name=subB_nsweep
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/subB_nsweep_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/subB_nsweep_%j.err

# Supplementary figure: hammock's --subB (mixed-stride) wall-time speedup as a
# function of N (files per side), a line per subB value, all measured with
# --metrics (jaccard_similarity_ie), p=18 -- see docs/seed-subsampling-synthetic-supplement.md
# for why this job exists (the outline's earlier synthetic-corpus subB claim,
# "4.21x at p=14", never traced to a real measurement at any precision).
#
# HAMMOCK-ONLY, NO BEDTOOLS PROJECTION. A prior design considered adding a
# projected bedtools reference line (the same log-log-fit technique
# sbatch_fig3_largeN.sh / plot_largeN_supplement.R use), but a reviewer's
# argument won out: one projected bedtools curve under THREE hammock curves
# lets a reader compute three different "hammock vs bedtools" speedups from a
# single unmeasured extrapolation -- triple the mis-citation surface of the
# existing large-N supplement, for zero new bedtools information (bedtools'
# wall time does not depend on subB at all, so the fitted curve would be
# numerically identical to the one Figure S7 already publishes). This figure
# answers a narrower, hammock-internal question instead: how does subB's own
# benefit change with scale.
#
# PREDICTION TO CHECK, STATED BEFORE MEASURING (both plan-review agents
# independently derived this): subB only discounts the Theta(N) sketch-
# construction phase; the Theta(N^2) pairwise-comparison phase (which
# --metrics makes costlier) is subB-invariant and its share of wall time grows
# fast -- measured (existing archived hammock_ie_B data) at ~3% of wall at
# N=64 vs ~35% at N=1024. So the three subB lines should NARROW toward each
# other as N -> 2048, not stay parallel. If they don't, that's the actual
# finding, not measurement noise.
#
# COST. hammock_ie_B (metrics-on, subB=1.0, p=18, t=16) sums to ~866s across
# N=2..2048 for ONE replicate of ONE arm (from archived cpp_vs_bedtools_t16_p18*
# data). subB<1.0 is NOT proportionally cheaper at this scale: subsampling only
# discounts sketch construction, and the pairwise phase dominates at high N
# regardless of subB (CLAUDE.md divergence #4's Maurano single-pair table shows
# subB=0.01 only 24% faster than subB=1.0 overall). Conservative estimate: 3
# arms x 5 reps x ~866s =~ 3.6h, hence the 6h allocation (padded).
#
# CORPUS. Matches Figure 3A / sbatch_fig3_largeN.sh exactly (10k intervals/file,
# 16 threads, p=18) so this job's own subB=1.0 arm can be sanity-checked
# against the archived N=512 point (docs/data/cpp_vs_bedtools_t16_p18.csv) --
# see the R analysis script for that check. Single allocation end-to-end (no
# stitching across jobs), so the two-allocation overlap-check machinery in
# plot_largeN_supplement.R does not apply here -- this is a free correctness
# check, not a required join.
#
# --metrics-all (new flag, see benchmark_cpp_vs_bedtools.py's arms_for()):
# runs EVERY subB arm with --metrics instead of the reduced-column timed arm
# (--register-equality, was --no-metrics), unlike --metrics-arm which only
# ever adds one extra fixed arm at subB=1.0.
# Labels: hammock_ie_B (subB=1.0, byte-compatible with the existing consumers'
# exact-match filter) / hammock_ie_B_subB0.1 / hammock_ie_B_subB0.01 (both fail
# every existing R consumer's filter, so this cannot corrupt an existing
# figure's data even if a script were pointed at this CSV by mistake).
#
# BINARY PINNING (added after job 29758041 was contaminated mid-run 2026-08-11
# by a concurrent rebuild of the shared site-packages hammock-cpp -- see
# docs/seed-subsampling-synthetic-supplement.md). This job copies the live binary to
# a job-local, node-scratch path ONCE at start and points HAMMOCK_CPP_BIN
# there, so nothing outside this job can mutate the binary this run measures,
# no matter how long the sweep takes or what else touches the shared env
# meanwhile.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
LIVE_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp
PINNED_DIR="${TMPDIR:-/tmp}/subB_nsweep_${SLURM_JOB_ID:-$$}"
mkdir -p "$PINNED_DIR"
trap 'rm -rf "$PINNED_DIR"' EXIT
cp "$LIVE_BIN" "$PINNED_DIR/hammock-cpp"
chmod +x "$PINNED_DIR/hammock-cpp"
export HAMMOCK_CPP_BIN="$PINNED_DIR/hammock-cpp"

echo "node: $(hostname)  cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)  cores: $(nproc)"
echo "pinned binary: $HAMMOCK_CPP_BIN (copied from $LIVE_BIN at job start, git HEAD $(git rev-parse HEAD))"
"$HAMMOCK_CPP_BIN" --version

"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --no-bedtools \
    --threads 16 \
    --precision 18 \
    --num-intervals 10000 \
    --num-files 2,4,8,16,32,64,128,256,512,1024,2048 \
    --subB-list 1.0,0.1,0.01 \
    --metrics-all \
    --corpus-seed 20260811 \
    --runs 5 \
    "$@"
