#!/bin/bash
#SBATCH --job-name=hm_pairwise_cost
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/pairwise_cost_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/pairwise_cost_%j.err

# Re-measure what the similarity block costs, as a function of precision.
#
# Why this exists as a batch job when the 20260804 run was not one: that run has
# no SLURM record, so it was taken on a shared dev node where an unrelated job
# can take cores from under it. That is not hypothetical -- a 16-thread run of
# this script on devlangmead1 on 2026-08-08 was sharing 48 cores with an
# 8-thread job the whole time. --cpus-per-task=16 reserves the cores instead.
#
# Config is deliberately identical to the 20260804 run (64 files/side, 10k
# intervals, threads 16, 5 runs, p 12..24) so the two are comparable except for
# the code change and the contention.
#
# NOTE on --threads 16 vs --cpus-per-task=16: as of the pairwise-phase fix,
# hammock honors --threads for the OpenMP pairwise loop via a num_threads()
# clause, so the team is actually 16 here. Before that fix the binary spawned a
# team per core on the node regardless, which is one more reason the old numbers
# are not a clean baseline.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

mkdir -p experiments/bedtools_benchmark/logs

# The refactor env: same env the binary was built in and the extension installed
# into. (The older sbatch scripts here run the harness under the *orig* conda
# env, which works only because the harness shells out and never imports
# hammock. Not worth imitating.)
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python3

# Pin the binary explicitly rather than letting find_hammock_cpp() glob the build
# tree by mtime. check_binary_version() probes whatever path is resolved and
# refuses anything older than 0.7.0.
export HAMMOCK_CPP_BIN=/home/jbonnie1/interval_sketch/hammock_claude/build/cp310-cp310-linux_x86_64/hammock-cpp

echo "host:    $(hostname)"
echo "cpu:     $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
echo "cpus:    ${SLURM_CPUS_PER_TASK:-unset} (cgroup: $(nproc) visible)"
echo "binary:  $HAMMOCK_CPP_BIN"

# -march=native is baked in (CMakeLists.txt:42,47,68) and the binary really does
# contain AVX-512 (zmm) instructions, so a node older than the build host dies
# with SIGILL. Fail here with a clear message instead of 70 confusing runs.
if ! "$HAMMOCK_CPP_BIN" --version >/dev/null 2>&1; then
    echo "FATAL: $HAMMOCK_CPP_BIN will not run on $(hostname)." >&2
    echo "  Most likely an -march=native/AVX-512 mismatch with the build host." >&2
    echo "  Rebuild on this node type, or submit with -w <node-like-build-host>." >&2
    exit 1
fi

"$PYTHON" experiments/bedtools_benchmark/pairwise_cost_by_precision.py \
    --precisions 12,14,16,18,20,22,24 \
    --num-files 64 \
    --num-intervals 10000 \
    --threads 16 \
    --runs 5
