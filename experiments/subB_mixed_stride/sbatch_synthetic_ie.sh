#!/bin/bash
#SBATCH --job-name=hm_subB_synth_ie
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=experiments/subB_mixed_stride/logs/sweep_synth_ie_%j.out
#SBATCH --error=experiments/subB_mixed_stride/logs/sweep_synth_ie_%j.err

# Synthetic-corpus counterpart to the Maurano jaccard_similarity_ie subB run
# (docs/data/maurano_subB_ie_summary.csv) that backs Figure 3B. Existing
# synthetic subB data (results/summary_synthetic_20260512_144455.csv) predates
# --metrics / jaccard_similarity_ie entirely (CLAUDE.md: metrics block is
# ~1.6x the wall cost of --no-metrics, so those timings aren't comparable to
# what Figure 3B reports). mixed-stride only, subB in {1.0, 0.1, 0.01} to
# match the Figure 3B bars exactly. All three size classes (10k/100k/1M) kept
# so we can see whether the speedup is stable within the synthetic corpus
# itself, same as the original table did.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

mkdir -p experiments/subB_mixed_stride/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3

"$PYTHON" experiments/subB_mixed_stride/run_sweep.py \
    --corpus synthetic \
    --methods mixed-stride \
    --subB-values 1.0 0.1 0.01 \
    --precision 18 \
    --threads 8 \
    --reps 5 \
    --metrics \
    --binary /home/jbonnie1/interval_sketch/hammock_claude/build/cp310-cp310-linux_x86_64/hammock-cpp
