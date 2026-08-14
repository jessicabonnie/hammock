#!/bin/bash
#SBATCH --job-name=cli_overhead_contention
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=experiments/bedtools_benchmark/logs/cli_overhead_contention_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/cli_overhead_contention_%j.err

# Quick, targeted follow-up to job 29763124's CPU-model confound (see
# docs/seed-hammock-cpp-file-dispatch.md "Part 2 update", "Still open").
#
# job 29763124's Arm B ran --partition=parallel --exclusive and landed on
# c192 (Xeon Gold 6248R) -- a different CPU model than job 29758101's sr15
# (Xeon Gold 6448Y), so the N=32/128 ratio shift vs. that original job
# couldn't be cleanly attributed to exclusivity alone.
#
# A truly clean same-node comparison (sbatch --partition=shared --nodelist=sr15
# --exclusive) is available but queues ~3 days even though sr15 is currently
# idle -- shared's exclusive whole-node requests appear to sit far back
# regardless of target node (checked via sbatch --test-only, not assumed).
# parallel's "c*"-named 48-core 6248R pool looks like a structurally separate
# hardware generation from shared's dedicated "sr*"-named 64-core 6448Y pool,
# so there's no way to get 6448Y-class hardware quickly through parallel.
#
# This job instead holds CPU model FIXED (parallel, same 6248R-class pool
# job 29763124 already used) and varies only exclusivity: no --exclusive
# here, a 16-core slice like the ORIGINAL shared job requested, so other
# jobs can land on and contend for the rest of the node. Reproduces Arm B's
# exact config otherwise (--cpp-metrics-arm metrics, N=32/128 only -- the two
# points the falsification criterion is actually about, corpus-seed 20260811,
# --runs 6) so it is directly comparable to job 29763124's Arm B rows at the
# same N values. If contention barely moves the ratio here (same CPU model as
# 29763124), that's evidence the shared-vs-exclusive shift versus the
# ORIGINAL job leans CPU-model-driven, not contention-driven -- and vice
# versa. Doesn't resolve the cross-CPU-model question by itself, only
# isolates the contention variable on one side of it.
#
# No binary pinning: short job (~10 min), low exposure window, not worth the
# complexity that the multi-hour job 29763124 needed it for.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
mkdir -p experiments/bedtools_benchmark/logs

echo "node: $(hostname)  cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)  cores requested: 16 of $(nproc) on node"

/home/jbonnie1/.conda/envs/hammock/bin/python3 -u experiments/bedtools_benchmark/measure_cli_overhead.py \
    --num-files 32,128 \
    --num-intervals 10000 \
    --precision 18 \
    --threads 16 \
    --runs 6 \
    --corpus-seed 20260811 \
    --cpp-metrics-arm metrics \
    --output experiments/bedtools_benchmark/results/cli_overhead_contention_check_${SLURM_JOB_ID}.csv
