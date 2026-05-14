#!/bin/bash
#SBATCH --job-name=maurano_dhs
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=experiments/maurano_dhs_validation/logs/sweep_%j.out
#SBATCH --error=experiments/maurano_dhs_validation/logs/sweep_%j.err

# Full A/B/C/D sweep on the Maurano fetal-tissue DHS corpus against bedtools
# pairwise Jaccard. See experiments/maurano_dhs_validation/README.md.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP_DIR=experiments/maurano_dhs_validation
mkdir -p "$EXP_DIR/logs"

# The refactored hammock package is installed in claude-ref-comparison
# (see memory/reference_new_hammock_env.md). Use that env's python so
# the driver scripts find hammock via DEFAULT_HAMMOCK.
CONDA_ENV=/home/jbonnie1/.conda/envs/claude-ref-comparison
PYTHON="$CONDA_ENV/bin/python3"
HAMMOCK="$CONDA_ENV/bin/hammock"

# 1. Prepare data (idempotent: symlinks BEDs/FASTAs, builds bedtools ref).
bash "$EXP_DIR/prepare_data.sh"

# 2. A/B/C sweep (HLL, 20 BEDs, p in {18,21,23}; Mode C expA + subB at p=21).
"$PYTHON" "$EXP_DIR/run_sweep_abc.py" --threads 8 --hammock "$HAMMOCK"

# 3. Mode D (k,w,p) sweep over 20 FASTAs. The full grid has 5 x 7 x 4 = 140
#    nominal configs (minus w<k); each run ~5-10 min at p=24.
"$PYTHON" "$EXP_DIR/run_sweep_d.py"   --threads 8 --hammock "$HAMMOCK"

# 4. R analysis (needs r/4.3.0; uses Cairo because system R has no X11).
ml r/4.3.0
Rscript "$EXP_DIR/analyze.R"
