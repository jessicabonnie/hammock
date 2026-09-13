#!/bin/bash
#SBATCH --job-name=hm_ms_vs_hash
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=experiments/mixed_stride_hash_threshold_speed/logs/run_%j.out
#SBATCH --error=experiments/mixed_stride_hash_threshold_speed/logs/run_%j.err

set -euo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-/home/jbonnie1/interval_sketch/hammock_claude}"
cd "$REPO_ROOT"

EXPERIMENT=experiments/mixed_stride_hash_threshold_speed
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
BINARY="$REPO_ROOT/build/cp310-cp310-linux_x86_64/hammock-cpp"
RUN_DIR="$EXPERIMENT/results/job_${SLURM_JOB_ID}"

mkdir -p "$EXPERIMENT/logs" "$EXPERIMENT/results"
export OMP_NUM_THREADS=8
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "host: $(hostname)"
echo "affinity cores: $(nproc)"
echo "commit: $(git rev-parse HEAD)"
"$BINARY" --version

"$PYTHON" "$EXPERIMENT/run_benchmark.py" \
  --binary "$BINARY" --output-dir "$RUN_DIR" \
  --rounds 10 --warmup-rounds 1 --precision 18 --threads 8 \
  --seed 42 --gate-seed 31337

"$PYTHON" "$EXPERIMENT/summarize.py" "$RUN_DIR"
