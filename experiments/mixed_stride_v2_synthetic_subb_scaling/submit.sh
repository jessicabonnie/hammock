#!/bin/bash
set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/mixed_stride_v2_synthetic_subb_scaling
mkdir -p "$EXP/logs" "$EXP/results/raw" "$EXP/figures"

job_benchmark=$(sbatch --parsable "$EXP/sbatch_benchmark.sh")
job_plot=$(sbatch --parsable --dependency="afterok:${job_benchmark}" \
  "$EXP/sbatch_finalize.sh")

echo "benchmark: $job_benchmark"
echo "finalize:  $job_plot (afterok:$job_benchmark)"
