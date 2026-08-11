#!/bin/bash
#SBATCH --job-name=cli_overhead
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/cli_overhead_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/cli_overhead_%j.err

# ⚠️  DO NOT EDIT python/hammock/, cpp/, or bindings/ IN THIS REPO WHILE THIS
# JOB IS RUNNING (~4-6h). The hammock CLI is a pip editable install, live-
# linked to python/hammock/ in this working tree -- any edit there changes
# what the CLI measures with NO rebuild needed, and the job's own
# provenance-drift detector (see measure_cli_overhead.py) can only catch it
# after the fact, at replicate granularity, and cannot see an edit that's
# reverted between two checks. This repo has already had same-day
# concurrent-session edits collide once during this job's own review
# process -- this is a real, not hypothetical, risk window.

# Quantifies the hammock CLI vs hammock-cpp overhead gap for the
# methods/limitations text -- see measure_cli_overhead.py's docstring and
# docs/seed-hammock-cpp-file-dispatch.md Part 2.
#
# Re-measurement of job 29758101 per Part 2's "Next steps" #2, now that Part 1
# (the CLI --threads sketch-phase bug) is fixed (v0.7.1). Two arms, same node,
# same SLURM allocation, same corpus seed (arm choice doesn't feed
# derive_seed, so both arms see byte-identical corpora at a given (N, run) --
# free replication control, not a bug):
#
#   ARM B (--metrics, N=2,8,32,128) -- reruns job 29758101's config (same
#   arm, same N range, same corpus seed/precision/threads) but on --exclusive
#   instead of --partition=shared. This is the one that actually answers the
#   seed's falsification criterion: "if --exclusive brings the N=32/128
#   ratios to within ±2-4% of 1.0 (job 29628907's noise floor), the crossover
#   is shared-partition contention, not a real dispatch-strategy difference."
#   The criterion was calibrated against the --metrics-arm ratios
#   (0.864/0.807) -- testing it on a different arm would answer a different,
#   uncalibrated question.
#
#   CAVEAT, found by a second adversarial review pass: job 29758101 ran on
#   node sr15 (Xeon Gold 6448Y). --partition=parallel --exclusive does NOT
#   pin the node -- precedent (job 29652432, same partition/exclusive
#   combo) landed on c594 (Xeon Gold 6248R), a different CPU generation.
#   CLAUDE.md already flags cpu_model as load-bearing (CMakeLists bakes in
#   -march=native). So a ratio shift in Arm B could be exclusivity, could be
#   CPU model, or both -- CHECK THE cpu_model COLUMN (now stamped on every
#   row, see below) before attributing any change to exclusivity alone.
#
#   ARM A (--no-metrics, N up to 2048) -- the arm benchmark_cpp_vs_bedtools.py
#   actually uses for the headline figures (job 29758101 was --metrics only,
#   which is NOT that arm), extended to N=512/1024/2048 where the paper's real
#   throughput claims live.
#
# Still only n=1 at the job/node level -- this does not by itself promote to
# CLAUDE.md's "Re-measured and gated" bar. If Arm B's ratios come back
# meaningfully different from 0.864/0.807 (crossover survives --exclusive) or
# meaningfully closer to 1.0 (crossover was contention), that's suggestive,
# not confirmed, until repeated on a second exclusive allocation -- see
# docs/seed-benchmark-methodology.md's job-level-variance caveat.
#
# --runs 6, not 5: with exactly 2 tools and the alternate-first-mover-by-run_i
# logic, an odd --runs gives a 3:2 split, not a balanced first-mover count
# (docs/seed-benchmark-methodology.md item 7 already flagged this exact class
# of bug and prescribes "--runs a multiple of the arm count").
#
# CHECKPOINTING + BINARY PINNING (added after a second review pass, mirroring
# benchmark_cpp_vs_bedtools.py's own checkpoint/pin fix, commits
# 4f288d3/b8a082a -- see docs/seed-subsampling-synthetic-supplement.md for
# the job-29758041 incident that motivated it: a concurrent binary rebuild
# contaminated a multi-hour sweep, and cancelling it lost every N value's
# results, not just the one in progress, because nothing had been written to
# disk yet). measure_cli_overhead.py now writes its CSV atomically after
# every num_files block, so a timeout/scancel/node failure during Arm A's
# multi-hour N=512-2048 phase costs at most the in-progress N -- Arm B and
# every completed N in Arm A survive on disk regardless. This job also pins
# hammock-cpp to a job-local scratch copy so a concurrent rebuild elsewhere
# can't change what this run is timing partway through, same as the sibling
# script. NOTE: the `hammock` CLI itself CANNOT be pinned the same way --
# claude-ref-comparison is a pip editable install, and prepending a job-local
# copy to PYTHONPATH does not shadow it (checked empirically, not assumed:
# the editable finder intercepts before normal sys.path resolution). So the
# CLI's Python-level behavior stays live-linked to whatever is on disk in
# this repo's python/hammock/ for this job's whole duration -- an edit there
# needs no rebuild to change what's being measured mid-run. The Python script
# now fingerprints (git HEAD, dirty status, both binaries' mtime/size) once
# at start and after every checkpoint, and stamps a provenance_drift column +
# a loud stderr WARNING if anything changed -- a detector, not a fix, for the
# one gap pinning can't close here.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
mkdir -p experiments/bedtools_benchmark/logs

echo "node: $(hostname)  cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)  cores: $(nproc)"

PYEXE="/home/jbonnie1/.conda/envs/hammock/bin/python3"
LIVE_CPP_BIN="/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp"
PINNED_DIR="${TMPDIR:-/tmp}/cli_overhead_${SLURM_JOB_ID:-$$}"
mkdir -p "$PINNED_DIR"
trap 'rm -rf "$PINNED_DIR"' EXIT
cp "$LIVE_CPP_BIN" "$PINNED_DIR/hammock-cpp"
chmod +x "$PINNED_DIR/hammock-cpp"
export HAMMOCK_CPP_BIN="$PINNED_DIR/hammock-cpp"

echo "pinned hammock-cpp: $HAMMOCK_CPP_BIN (copied from $LIVE_CPP_BIN at job start, git HEAD $(git rev-parse HEAD))"
"$HAMMOCK_CPP_BIN" --version

echo "=== Arm B: --metrics control, N=2,8,32,128 (job 29758101's config, on --exclusive) ==="
"$PYEXE" -u experiments/bedtools_benchmark/measure_cli_overhead.py \
    --num-files 2,8,32,128 \
    --num-intervals 10000 \
    --precision 18 \
    --threads 16 \
    --runs 6 \
    --corpus-seed 20260811 \
    --cpp-metrics-arm metrics \
    --output experiments/bedtools_benchmark/results/cli_overhead_metrics_exclusive_${SLURM_JOB_ID}.csv

echo "=== Arm A: --no-metrics, N=2,8,32,128,512,1024,2048 (headline-figure arm, extended N) ==="
"$PYEXE" -u experiments/bedtools_benchmark/measure_cli_overhead.py \
    --num-files 2,8,32,128,512,1024,2048 \
    --num-intervals 10000 \
    --precision 18 \
    --threads 16 \
    --runs 6 \
    --corpus-seed 20260811 \
    --cpp-metrics-arm no-metrics \
    --output experiments/bedtools_benchmark/results/cli_overhead_nometrics_exclusive_${SLURM_JOB_ID}.csv
