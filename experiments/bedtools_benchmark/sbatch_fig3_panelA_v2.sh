#!/bin/bash
#SBATCH --job-name=fig3_panelA_v2
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fig3_panelA_v2_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fig3_panelA_v2_%j.err

# Panel A, second revision, needed because the low-N fix cannot be spliced
# across the earlier two jobs. job 29656140 (N=2..512, n=3) already fixed the
# LD_LIBRARY_PATH bug, but at n=3 the small-N cells were noise-dominated: N=2's
# three reps spanned 0.14-1.22s and reported hammock 2.02x faster; a 20-rep
# recheck (job 29670970, N<=32 only) found the true, tight answer is the
# OPPOSITE -- bedtools 1.53x faster, 20/20. That recheck ran on node c431,
# job 29656140 ran on c531, and their one overlapping cell (N=32) disagrees by
# 28% -- 3.35s vs 2.61s -- which is bigger than either job's own noise, so it
# is very likely a genuine node effect (consistent with this repo's own
# previously-documented 1.17-2.86x node-to-node spread on identical bedtools
# code, docs/bedtools-parallelism-caveat.md). Splicing job 29670970's N<=32
# rows into job 29656140's N>=64 rows would silently introduce that node jump
# at exactly the crossover the fix was trying to pin down.
#
# So: one job, one node, two passes at two replicate counts, matching what the
# noise actually requires. N in {2,4,8,16,32} needs many reps because the
# cells are cheap and were shown to be wrong at n=3; N in {64,...,512} does NOT
# need more reps -- those cells already had low CV at n=3 in job 29656140 (e.g.
# N=512's three reps spanned 714.4-718.5s, under 1%) and 20 reps there would
# cost about 2.4 hours for no benefit. Both passes carry --metrics-arm so the
# +IE curve stays complete across the whole N range.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
ml bedtools/2.30.0
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

echo "node: $(hostname)  cores: $(nproc)"
"$HAMMOCK_CPP_BIN" --version

# No --output flag exists (only --output-dir; the stem is
# cpp_vs_bedtools_t<threads>_<timestamp>, auto-generated). Both passes land in
# the same results dir with different timestamps, low-N first, so they are
# distinguishable afterward by mtime order -- no renaming needed here.

echo; echo "================ low N, n=20 ================"
"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --threads 16 --precision 18 \
    --num-files 2,4,8,16,32 \
    --metrics-arm \
    --runs 20 \
    --corpus-seed 20260810 \
    "$@"

echo; echo "================ high N, n=3 ================"
"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --threads 16 --precision 18 \
    --num-files 64,128,256,512 \
    --metrics-arm \
    --runs 3 \
    --corpus-seed 20260810 \
    "$@"

echo; echo "node for both passes: $(hostname)"
