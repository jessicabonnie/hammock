#!/bin/bash
#SBATCH --job-name=fig3_panelB
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fig3_panelB_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fig3_panelB_%j.err

# Figure 3, Panel B: the precision frontier on the Maurano DHS corpus.
# x = Jaccard MAE vs exact bedtools, y = speedup vs bedtools, one point per p.
#
# ON MAURANO, NOT SYNTHETIC. Two reasons, the second decisive:
#
#   1. Speed and accuracy then come from the same corpus in the same job. The
#      alternative splices Maurano accuracy onto synthetic speed.
#   2. The synthetic corpus cannot carry an accuracy axis at all. Its generator
#      draws start in [0, 10 Mbp] over 24 chromosomes, so every file covers
#      ~19% of the same span independently and essentially every pair has the
#      same true J ~ 0.10 -- a delta function, with no range for a
#      MAE-vs-precision curve to be read against. Maurano spans J = 0.135-0.627.
#
# BOTH THREAD SETTINGS, IN ONE JOB. Panel A is synthetic at t=16; Panel C (the
# existing Maurano subB bars) is t=8. Running both here means Panel B can be
# drawn at t=16 to match Panel A while the t=8 series ties to Panel C, instead
# of Figure 3 carrying two thread counts across three panels with no bridge.
# Because both run inside THIS allocation, t=8 vs t=16 is a within-job
# comparison and survives the machine factor that kills cross-run absolutes.
#
# 20 files passed as both operands (the subB_mixed_stride convention) = 400
# ordered pairs, all of them timed. 20 are file-vs-itself and are excluded from
# MAE -- both tools return ~1.0 there at zero error, which is free correctness
# that dilutes MAE ~5% without measuring anything. So n_pairs reads 380.
#
# VERIFICATION GATE. At p=18 the run must reproduce
# jaccard_ie_mae_vs_bt = 1.1517e-3, which was computed independently from the
# archived docs/data/hammock_hll_p18_jaccB.csv against maurano_bedtools_ref.tsv.
# A local dry run of this exact configuration returned 1.151647e-3 (ratio
# 1.0000). A materially different value means something is wrong with the run,
# not with the gate.
#
# COST. ~208 s per replicate per thread setting -- 7 precisions x 2 passes,
# since the untimed --metrics pass re-sketches everything -- plus ~12 s of
# bedtools. About 12 min per thread setting, ~25 min total. The 4 h limit is
# slack.
#
# READING THE OUTPUT. Panel B's y-axis is speedup vs bedtools, the same
# quantity Panel A annotates, and at p=18 it will read ~1.16x while Panel A
# carries an ~8x arrow. That is not a contradiction: Maurano is N=20, well below
# the N~64 crossover. Panel B's p=18 point and Panel C's first bar are in fact
# the same measurement. Annotate it rather than leaving it to the caption.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

# Preflight -- see the note in sbatch_fig3_panelA.sh.
echo "node:     $(hostname)"
echo "cpu:      $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"
echo "cores:    $(nproc)"
"$HAMMOCK_CPP_BIN" --version || {
    echo "hammock-cpp failed to run on this node -- likely -march=native vs a" >&2
    echo "different CPU model. Rebuild on this node type or pin -march." >&2
    exit 1
}

for T in 16 8; do
    echo
    echo "================ threads = $T ================"
    "$PYTHON" experiments/bedtools_benchmark/sweep.py \
        --axis precision \
        --corpus maurano \
        --precisions 12,14,16,18,20,22,24 \
        --threads "$T" \
        --subB-list 1.0 \
        --runs 3 \
        "$@"
done
