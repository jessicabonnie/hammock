#!/bin/bash
#SBATCH --job-name=maurano_fill_highk_w
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=experiments/maurano_dhs_validation/logs/fill_highk_w_%j.out
#SBATCH --error=experiments/maurano_dhs_validation/logs/fill_highk_w_%j.err

# Scoped fill-in: Mode D sweep for the missing high-k x large-w cells.
# Grid: k in {15,20,25} x w in {300,500} x p in {20,22,23,24} = 24 invocations
# (all satisfy w >= k). These complete the raw_d set so analyze.R can be
# regenerated over a uniform grid. See experiments/maurano_dhs_validation/README.md.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP_DIR=experiments/maurano_dhs_validation
mkdir -p "$EXP_DIR/logs"

# Refactored hammock package is installed in claude-ref-comparison
# (see memory/reference_new_hammock_env.md). Use that env's python3 and hammock.
CONDA_ENV=/home/jbonnie1/.conda/envs/claude-ref-comparison
PYTHON="$CONDA_ENV/bin/python3"
HAMMOCK="$CONDA_ENV/bin/hammock"

# 1. Run the 24 missing Mode D cells.
"$PYTHON" "$EXP_DIR/run_sweep_d.py" \
    --ks 15 20 25 \
    --ws 300 500 \
    --ps 20 22 23 24 \
    --threads 8 \
    --hammock "$HAMMOCK"

# 2. R analysis over the now-complete raw_d set (needs r/4.3.0; Cairo for PNG).
ml r/4.3.0
Rscript "$EXP_DIR/analyze.R"
