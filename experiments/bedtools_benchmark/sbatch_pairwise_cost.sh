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
# NOTE on --threads 16 vs --cpus-per-task=16: hammock-cpp has called
# omp_set_num_threads(args.threads) since the initial commit (hammock_cli.cpp),
# so -t 16 has always meant a 16-thread team here, sketching and pairwise alike.
# The num_threads() clause added alongside the fused pass fixed the *Python*
# path only (bindings/_core.cpp), which had ignored --threads for the pairwise
# phase. So the 20260804 baseline was NOT taken with a 48-wide team, and nothing
# about the thread count distinguishes it from this run -- the differences are
# the fused code and the reserved cores.

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
# with SIGILL. Fail here with a clear message instead of mid-sweep.
#
# --version alone is NOT a sufficient probe: it parses argv and prints a string,
# and the loader performs no ISA check, so an unsupported instruction faults only
# when the sketching/pairwise kernel actually executes it. Check the CPU flag
# directly, then force a real end-to-end run through those kernels.
if ! grep -qm1 avx512f /proc/cpuinfo; then
    echo "FATAL: $(hostname) lacks avx512f; the -march=native binary will SIGILL." >&2
    echo "  Rebuild on this node type, or target a node matching the build host." >&2
    exit 1
fi
_probe=$(mktemp -d)
trap 'rm -rf "$_probe"' EXIT
printf 'chr1\t100\t200\nchr1\t500\t900\n' > "$_probe/a.bed"
printf 'chr1\t150\t250\nchr2\t10\t99\n'   > "$_probe/b.bed"
printf '%s\n' "$_probe/a.bed" "$_probe/b.bed" > "$_probe/list.txt"
if ! "$HAMMOCK_CPP_BIN" "$_probe/list.txt" "$_probe/list.txt" \
        --mode B -p 14 -t 2 -o "$_probe/out" >/dev/null 2>&1; then
    echo "FATAL: $HAMMOCK_CPP_BIN failed a smoke run on $(hostname)." >&2
    echo "  Exercises the real sketching+pairwise kernels, so this catches an" >&2
    echo "  -march=native mismatch that --version would not." >&2
    exit 1
fi

"$PYTHON" experiments/bedtools_benchmark/pairwise_cost_by_precision.py \
    --precisions 12,14,16,18,20,22,24 \
    --num-files 64 \
    --num-intervals 10000 \
    --threads 16 \
    --runs 5
