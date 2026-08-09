#!/bin/bash
#SBATCH --job-name=fig3_panelA
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fig3_panelA_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fig3_panelA_%j.err

# Figure 3, Panel A: wall time vs N on the synthetic corpus, at p=18.
#
# WHY THIS EXISTS. Panel A has always been measured at p=14, which appears
# nowhere else in the paper: the CLI default is 18, Panel B is 18, the accuracy
# figure is 18/21/23. Runtime is a steep function of precision -- across the
# archived sweep, sketching grows 18.7x and the pairwise phase 967.7x from p=14
# to p=24 -- so the speed claim was being read off the cheapest configuration
# while the accuracy claims came from a different one. Nothing else about the
# protocol changes; this is the same harness at the default precision.
#
# It was never caught because neither archived dataset could answer it: the
# precision sweep that HAS a bedtools arm stops at p=18, and the run that
# reaches p=24 has NO bedtools arm.
#
# --exclusive --mem=0 rather than the old --partition=shared --cpus-per-task=16
# --mem=32G. The measured 24% gap between two harnesses on the same node on the
# same day is confined entirely to the sketch phase, and co-tenancy is one of
# the two surviving explanations; an exclusive node removes it as a variable.
# The workload still uses 16 threads -- the other 32 cores are guaranteed idle
# headroom, not extra parallelism.
#
# ARMS. bedtools, hammock subB=1.0 --no-metrics, hammock subB=1.0 with the
# default 9-column block (--metrics-arm). Deliberately NO subB<1 arm: Panel C is
# the subsampling panel, plot_pairwise_scaling.R:112 filters subB curves out of
# Panel A anyway, and including one re-opens a framing this figure is being
# rebuilt to avoid.
#
# All three arms rotate by replicate index, bedtools included. Until recently
# bedtools was exempt and ran first every time, immediately behind the pre-sort
# that had just walked every input into page cache.
#
# COST. ~1318 s per replicate including corpus generation (42 s of pure-Python
# generate_bed_file) and the pre-sort (11 s), so ~1.1 h for 3 replicates. The
# 6 h limit is slack, not an estimate.
#
# READING THE OUTPUT. Report within-job ratios. Do not divide an absolute from
# this job by an absolute from Panel B's -- they are separate allocations and
# may be separate nodes, and the binary is built -march=native.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

# Preflight. The binary is compiled -march=native on a login node and this
# cluster exposes no CPU-model SLURM feature, so --constraint cannot pin a
# homogeneous set and an allocation on an older microarchitecture would die of
# SIGILL. Better to find that out in the first second than 40 minutes in.
echo "node:     $(hostname)"
echo "cpu:      $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
echo "cores:    $(nproc)"
"$HAMMOCK_CPP_BIN" --version || {
    echo "hammock-cpp failed to run on this node -- likely -march=native vs a" >&2
    echo "different CPU model. Rebuild on this node type or pin -march." >&2
    exit 1
}

"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --threads 16 \
    --precision 18 \
    --num-intervals 10000 \
    --num-files 2,4,8,16,32,64,128,256,512 \
    --subB-list 1.0 \
    --metrics-arm \
    --corpus-seed 20260809 \
    --runs 3 \
    "$@"
