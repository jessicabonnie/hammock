#!/bin/bash
set -euo pipefail
: "${HAMMOCK_REHASH_RUN_ROOT:?set HAMMOCK_REHASH_RUN_ROOT}"
: "${HAMMOCK_REHASH_SOURCE_COMMIT:?set HAMMOCK_REHASH_SOURCE_COMMIT}"
: "${SLURM_ARRAY_TASK_ID:?submit as an array}"
phase="${HAMMOCK_REHASH_PHASE:-primary}"
repo="${HAMMOCK_REHASH_RUN_ROOT}/hammock"
experiment="${repo}/experiments/maurano_dhs_rehash_validation"
export PYTHONPATH="${HAMMOCK_REHASH_RUN_ROOT}/runtime/site-packages"
export LD_PRELOAD="/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/libstdc++.so.6"
mkdir -p "${experiment}/logs/slurm"
cd "${repo}"
exec /usr/bin/time -v -o "${experiment}/logs/slurm/rehash_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.time" \
  /home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python \
  "${experiment}/scripts/run_rehash_sweep.py" --config "${experiment}/config.yaml" \
  --hammock "${experiment}/workflow/hammock_runtime.sh" \
  --source-commit "${HAMMOCK_REHASH_SOURCE_COMMIT}" --phase "${phase}" \
  --job-index "${SLURM_ARRAY_TASK_ID}"
