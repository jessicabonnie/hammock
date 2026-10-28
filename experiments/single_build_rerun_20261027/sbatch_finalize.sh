#!/bin/bash
#SBATCH --job-name=sbr_finalize
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/finalize_%j.out
#SBATCH --error=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/finalize_%j.err

# Side-by-side summaries (archived 0.7.0 / 0.5.0 vs 0.12.0) into
# summary/finalize_<jobid>/. Submitted with afterany so it runs even if one
# benchmark fails; missing benchmarks are recorded in NOTES.txt.

set -euo pipefail
source /home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/common.sh
cd "$EXP"
"$PY_REF" finalize.py
