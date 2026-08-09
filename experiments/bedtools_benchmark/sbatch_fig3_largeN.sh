#!/bin/bash
#SBATCH --job-name=fig3_largeN
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fig3_largeN_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fig3_largeN_%j.err

# hammock at catalog scale: N = 512, 1024, 2048 files per side, NO bedtools arm.
#
# WHY NO BASELINE. bedtools is Theta(N^2) with a large constant -- measured 3.97x
# per doubling, 495 s per replicate at N=256 -- so one replicate costs ~2.2 h at
# N=1024 and ~8.6 h at N=2048. Three replicates at N=2048 is over a day of node
# time to strengthen a claim that is already settled: Panel A measures 14.75x at
# N=256 and the gap is still widening. The affordable move is to measure hammock
# where the exact tool cannot follow and project bedtools from its own fitted
# curve, LABELLED AS A PROJECTION.
#
# N=512 is included deliberately, even though Panel A also measures it. It is the
# overlap point: if this job's N=512 hammock number disagrees with Panel A's,
# something about the two allocations differs and the extrapolation should not be
# trusted. A join with no shared point cannot be checked at all.
#
# WHY THIS IS CHEAP. At p=18 hammock is sketching-dominated, i.e. Theta(N) files
# rather than Theta(N^2) pairs: measured 15.80 s at N=128 and 33.57 s at N=256.
# The pairwise term does grow quadratically but from a tiny base (~12.6 us/pair
# at p=18), so it only starts to matter around here -- which is itself worth
# capturing, since N=2048 is 4.2M pairs and is where the linear approximation
# should visibly break.
#
# Corpus size check, because this is where it would bite: 2 x 2048 files x 10k
# intervals is ~873 MB of BED in the node-local /tmp (881 GB free, xfs on NVMe),
# and 4096 sketches at p=18 is ~1 GB resident. Neither is close to a limit.
#
# Both arms are kept (--no-metrics and the default 9-column block) because the
# cost of the metrics block at large N has never been measured -- every number
# in CLAUDE.md's table is N=64.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

echo "node: $(hostname)  cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)  cores: $(nproc)"
echo "tmp:  $(df -hT /tmp | tail -1)"
"$HAMMOCK_CPP_BIN" --version

"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --no-bedtools \
    --threads 16 \
    --precision 18 \
    --num-intervals 10000 \
    --num-files 512,1024,2048 \
    --subB-list 1.0 \
    --metrics-arm \
    --corpus-seed 20260809 \
    --runs 3 \
    "$@"
