#!/bin/bash
set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/mixed_stride_v2_figure3_reproduction
mkdir -p "$EXP/logs" "$EXP/results/panel_b" "$EXP/figures"

job_b=$(sbatch --parsable "$EXP/sbatch_panel_b.sh")
job_plot=$(sbatch --parsable --dependency="afterok:${job_b}" \
  "$EXP/sbatch_finalize.sh")

echo "panel B:  $job_b"
echo "finalize: $job_plot (afterok:$job_b)"
