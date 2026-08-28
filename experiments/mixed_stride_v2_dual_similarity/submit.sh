#!/bin/bash
set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/mixed_stride_v2_dual_similarity
mkdir -p "$EXP/logs" "$EXP/results" "$EXP/figures"

job_re=$(sbatch --parsable "$EXP/sbatch_register_equality.sh")
job_plot=$(sbatch --parsable --dependency="afterok:${job_re}" \
  "$EXP/sbatch_finalize.sh")

echo "register-equality: $job_re"
echo "finalize:          $job_plot (afterok:$job_re)"
