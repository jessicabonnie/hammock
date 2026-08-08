#!/bin/bash
#SBATCH --job-name=hm_subB_maurano
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=experiments/subB_mixed_stride/logs/maurano_%j.out
#SBATCH --error=experiments/subB_mixed_stride/logs/maurano_%j.err

# Panel B of paper/figures/pairwise_scaling.png: the Maurano subB tradeoff and
# its bedtools baseline.
#
# Why re-run something the code change did not touch. Both legs are --no-metrics
# throughout, so the fused pairwise pass cannot have moved them. What is wrong
# with the existing data is where it was taken: sacct has no job for either leg,
# and sweep_maurano_20260512_173720.log records cpu_count=48 / 1511 GB, which is
# devlangmead1 -- a shared interactive node with no reserved cores. Panel A came
# from an sr-class allocation, so the two panels of one figure were measured
# under different contention regimes, and neither CSV recorded enough to notice.
#
# Also a pre-0.7.0 binary: it reported integer milliseconds, so the pair-phase
# column had ~1 ms granularity against a 14 ms measurement. The current binary
# reports microseconds.
#
# Threads stay at 8, matching the May run, NOT Panel A's 16. The panels answer
# different questions and are never divided by one another; changing it would
# break continuity with the subB tables in this experiment's RESULTS.md and in
# CLAUDE.md divergence #4, for no gain.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools2
ml parallel
mkdir -p experiments/subB_mixed_stride/logs

PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/interval_sketch/hammock_claude/build/cp310-cp310-linux_x86_64/hammock-cpp

echo "host:   $(hostname)"
echo "cpu:    $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
echo "cpus:   ${SLURM_CPUS_PER_TASK:-unset} reserved, $(nproc) visible"
echo "binary: $HAMMOCK_CPP_BIN"

if ! "$HAMMOCK_CPP_BIN" --version >/dev/null 2>&1; then
    echo "FATAL: $HAMMOCK_CPP_BIN will not run on $(hostname) (-march=native?)" >&2
    exit 1
fi

# hammock leg: 3 methods x 6 subB x 5 reps on the 20-sample DHS corpus.
"$PYTHON" experiments/subB_mixed_stride/run_sweep.py --corpus maurano

# Same grid again with the 9-column output (jaccard_similarity_ie plus the
# containment/cosketch block) -- the configuration CLAUDE.md actually
# recommends, which no published Panel B number has ever measured. Kept as a
# separate arm rather than a replacement so the two are directly comparable and
# switching the panel stays a decision, not a side effect: it lands in
# sweep_maurano_ie_<stamp>.csv, which fails analyze.R's ^sweep_maurano_[0-9]
# glob and so cannot be auto-adopted as the headline sweep.
"$PYTHON" experiments/subB_mixed_stride/run_sweep.py --corpus maurano --metrics

# bedtools baseline: 190 unique pairs x 5 reps, 8-way GNU parallel fan-out.
# Shared by both arms -- bedtools is indifferent to hammock's output shape.
"$PYTHON" experiments/subB_mixed_stride/run_bedtools_maurano.py --threads 8 --reps 5
