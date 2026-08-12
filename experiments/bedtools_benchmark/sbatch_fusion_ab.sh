#!/bin/bash
#SBATCH --job-name=hm_fusion_ab
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fusion_ab_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fusion_ab_%j.err

# Settles what three separate jobs across four days could not: how much of the
# similarity block's cost the fused pass actually removed.
#
# The previous answer came from comparing a 2026-08-04 run to a 2026-08-08 run.
# That comparison is unsound, and its own control says so -- the reduced-column
# arm (--no-metrics on the pre binary, --register-equality on post -- see
# fusion_ab.py's capability-probe docstring) is byte-identical code in both
# binaries and it measured 1.27-1.53x slower on the older run (which, per
# sacct, had no SLURM allocation at all). Any absolute cross-run number
# therefore carries an unknown machine factor of that size.
#
# Here all four arms -- {pre, post} x {--metrics, reduced-column} -- run against
# one seeded corpus in one allocation on one node, permuted per replicate so arm
# is not confounded with position. Every reported quantity is a within-replicate
# ratio, so a machine factor is common-mode and cancels.
#
# The pre/post control on the reduced-column arm is the calibration: identical
# code, so it must land at 1.00, and however far it misses is the error bar on
# everything else in the job.
#
# Binaries are built from 826d7b4^ (pre) and HEAD (post) with the SAME compiler
# and cmake flags -- otherwise the comparison measures the build, not the change.
# They live in build/ab/ because /tmp is node-local and a compute node cannot see
# the login node's copy.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python3
PRE=/home/jbonnie1/interval_sketch/hammock_claude/build/ab/bin_pre
POST=/home/jbonnie1/interval_sketch/hammock_claude/build/ab/bin_post

echo "host:   $(hostname)"
echo "cpu:    $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
echo "cpus:   ${SLURM_CPUS_PER_TASK:-unset} reserved, $(nproc) visible"
echo "pre :   $PRE  md5=$(md5sum "$PRE" | cut -c1-12)"
echo "post:   $POST md5=$(md5sum "$POST" | cut -c1-12)"

# -march=native is baked in and the binaries carry AVX-512. --version does not
# execute the vector kernels and the loader performs no ISA check, so probe the
# CPU flag and then force a real run through the sketching+pairwise path.
if ! grep -qm1 avx512f /proc/cpuinfo; then
    echo "FATAL: $(hostname) lacks avx512f; these binaries will SIGILL." >&2
    exit 1
fi
_probe=$(mktemp -d)
trap 'rm -rf "$_probe"' EXIT
for i in $(seq 1 8); do
    printf 'chr1\t%d\t%d\nchr2\t%d\t%d\n' $((i*100)) $((i*100+50)) $((i*7)) $((i*7+30)) \
        > "$_probe/f$i.bed"
    echo "$_probe/f$i.bed" >> "$_probe/list.txt"
done
for b in "$PRE" "$POST"; do
    if ! "$b" "$_probe/list.txt" "$_probe/list.txt" --mode B -p 14 -t 2 \
            -o "$_probe/out" >/dev/null 2>&1; then
        echo "FATAL: $b failed a smoke run on $(hostname) (-march=native mismatch?)" >&2
        exit 1
    fi
done

"$PYTHON" experiments/bedtools_benchmark/fusion_ab.py \
    --pre-binary "$PRE" \
    --post-binary "$POST" \
    --precisions 12,14,16,18,20,22,24 \
    --num-files 64 \
    --num-intervals 10000 \
    --threads 16 \
    --runs 5 \
    --corpus-seed 20260808
