#!/bin/bash
#SBATCH --job-name=fig3_panelA
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
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
# ARMS. bedtools, hammock subB=1.0 --register-equality (the reduced-column
# timed arm), hammock subB=1.0 with the full metrics block (--metrics-arm).
# Deliberately NO subB<1 arm: Panel C is
# the subsampling panel, plot_pairwise_scaling.R:112 filters subB curves out of
# Panel A anyway, and including one re-opens a framing this figure is being
# rebuilt to avoid.
#
# All three arms rotate by replicate index, bedtools included. Until recently
# bedtools was exempt and ran first every time, immediately behind the pre-sort
# that had just walked every input into page cache.
#
# THE BEDTOOLS BASELINE IS THE FRAGILE PART OF THIS JOB. `bedtools jaccard` has
# no batch mode, so the workflow launches one process per pair -- N^2 of them --
# and on these nodes process creation caps near 123 exec/s and does not scale
# with cores. Job 29651772 was cancelled for exactly this: its bedtools leg ran
# at ~0.8x parallel efficiency, i.e. "t=16" meant "t<1", which would have
# inflated the reported speedup by roughly 6x.
#
# Two mitigations, both in place before this job:
#   * bedtools.sh now runs ONE process per pair (--tagstring + a single awk)
#     instead of three. Worth ~2.1x, consistently, on every config tested.
#   * every bedtools row carries mean_bedtools_parallel_eff, so the achieved
#     parallelism is recorded rather than assumed. READ IT before quoting any
#     speedup from this run.
#
# The achieved level is node-dependent and cannot be fixed from here -- measured
# 1.17x (shared/sr08), 1.68x (parallel/exclusive/c599), 2.86x
# (parallel/cpus-per-task=16/c516), all with the new dispatch. That spread is
# why the number travels in the CSV instead of being assumed to be ~1.
#
# COST: was ~4.6 h for 3 replicates on job 29651772; expect roughly half that
# now that the bedtools leg is ~2.1x faster per pair. MEASURED there -- an
# earlier version of this header said "~1318 s per replicate, so ~1.1 h", which
# was wrong by 4x.
#
# The error is worth recording because it is structural, not arithmetic: the job
# is bedtools' N^2 term at the largest N and essentially nothing else. Measured
# bedtools wall per replicate at t=16, 10k intervals/file:
#
#   N        2     4     8    16    32    64   128   256*   512*
#   wall  0.86  0.86  1.58  4.47 15.88 61.70 244.8  ~972  ~3859     (* projected)
#
# per-doubling ratio 2.83 -> 3.55 -> 3.89 -> 3.97, converging up to the
# theoretical 4.0. So N=512 alone is ~69% of the job and N in {256,512} is ~92%:
# every N below 128 is free, and the total is set almost entirely by one cell.
# Any estimate that does not start from the largest-N bedtools term is wrong.
#
# hammock is not the constraint and does not scale the same way -- it goes
# 7.33 s (N=64) -> 14.75 s (N=128), roughly linear, because at p=18 it is
# sketching-dominated (O(N) files) rather than pairwise-dominated.
#
# --time=12:00:00, not 6. The harness aggregates and writes its CSV only at the
# END, so a timeout loses the whole run rather than truncating it. 6 h did fit
# (~1.35 h margin), but the asymmetry between "a few idle hours in the
# reservation" and "lose 4.6 h of compute and start over" is not close.
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
