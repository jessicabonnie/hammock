#!/bin/bash
#SBATCH --job-name=fig3b_vec
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=experiments/vectorized_rebuild_20261027/logs/fig3b_%j.out
#SBATCH --error=experiments/vectorized_rebuild_20261027/logs/fig3b_%j.err

# Figure 3B rerun with the db1411a build compiled on a 6248R node.
# Protocol identical to experiments/mixed_stride_v2_figure3_reproduction.

set -euo pipefail
cpu=$(lscpu | sed -n 's/^Model name:[[:space:]]*//p')
[[ "$cpu" == *"6248R"* ]] || { echo "unexpected CPU: $cpu" >&2; exit 1; }
cd /home/jbonnie1/interval_sketch/hammock_claude

ml bedtools/2.30.0
ml parallel

EXP=experiments/vectorized_rebuild_20261027/fig3b
PYTHON=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python
BINARY=/vast/blangme2/jbonnie/hammock-db1411a-builds/c6248r-gcc14/hammock-cpp
EXPECTED_SHA=$(sed -n "s/^binary_sha256: //p" experiments/vectorized_rebuild_20261027/build_metadata.txt)

actual_sha=$(sha256sum "$BINARY" | cut -d' ' -f1)
if [[ "$actual_sha" != "$EXPECTED_SHA" ]]; then
  echo "binary sha256 mismatch: $actual_sha" >&2; exit 1
fi

{
  echo "job_id: ${SLURM_JOB_ID}"
  echo "node: $(hostname)"
  echo "cores: $(nproc)"
  echo "cpu: $(lscpu | sed -n 's/^Model name:[[:space:]]*//p')"
  echo "benchmark_commit: $(git rev-parse HEAD)"
  echo "binary: $BINARY"
  echo "binary_sha256: $actual_sha"
  echo "binary_version: $("$BINARY" --version)"
  echo "bedtools_version: $(bedtools --version)"
  echo "parallel_version: $(parallel --version | head -1)"
  echo "start: $(date -Is)"
} | tee "$EXP/results/metadata.txt"

"$PYTHON" experiments/subB_mixed_stride/run_sweep.py \
  --binary "$BINARY" --corpus maurano --methods mixed-stride \
  --subB-values 1.0 0.3 0.1 0.03 0.01 0.008 0.005 0.003 0.001 \
  --reps 5 --precision 18 --threads 8 --metrics \
  --output "$EXP/results/panel_b/hammock.csv"

"$PYTHON" experiments/subB_mixed_stride/run_bedtools_maurano.py \
  --threads 8 --reps 5 \
  --output "$EXP/results/panel_b/bedtools.csv"

echo "end: $(date -Is)" | tee -a "$EXP/results/metadata.txt"
