#!/bin/bash
#SBATCH --job-name=fig3_threads
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fig3_threads_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fig3_threads_%j.err

# Supplementary threading figure: how each tool uses the cores it is given.
#
# WHY. Figure 3's main panels quote a speedup at a fixed thread count, which
# silently assumes both tools converted cores into throughput. They do not, and
# the difference is structural rather than incidental:
#
#   * hammock parallelises INSIDE one process (OpenMP over a shared sketch
#     pool), reading each input file once.
#   * bedtools has no batch mode, so a pairwise workflow is N^2 independent
#     process launches under a GNU parallel wrapper, re-reading each file N
#     times.
#
# The previous version of this comment attributed the resulting slowdown to a
# "~123 exec/s process creation cap" and cited an md5sum control as
# independent confirmation. RETRACTED 2026-08-09: that control ran in the same
# polluted environment as the bug it was meant to rule out (the bedtools
# module's LD_LIBRARY_PATH carries 13 unused, GPFS-hosted directories that
# ld.so searched on every exec -- docs/bedtools-baseline-retraction.md), so it
# confirmed nothing. This run uses the now-fixed sweep.py (imports the patched
# run_bedtools(), which resolves a minimal LD_LIBRARY_PATH once) and is what
# actually measures bedtools' achieved parallelism rather than asserting a
# mechanism for it. Read the result; do not assume the old number.
#
# The threads axis is self-calibrating: t=1 is one of its own points, so
# efficiency is speedup(t)/t computed entirely within one allocation. No
# cross-run comparison is involved, which is what makes it defensible given
# that absolutes on this cluster are not portable between runs.
#
# Config matches Panel A's regime so the two are readable together: synthetic,
# 64 files/side, 10k intervals/file, p=18, subB=1.0. bedtools.sh is the fixed
# one-process-per-pair version; an older three-process-per-pair version would
# understate the baseline by ~2.1x.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

echo "node: $(hostname)  cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)  cores: $(nproc)"
"$HAMMOCK_CPP_BIN" --version

"$PYTHON" experiments/bedtools_benchmark/sweep.py \
    --axis threads \
    --thread-list 1,2,4,8,16,32,48 \
    --precision 18 \
    --num-files 64 \
    --num-intervals 10000 \
    --subB-list 1.0 \
    --runs 3 \
    "$@"
