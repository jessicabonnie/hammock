#!/bin/bash
set -euo pipefail
: "${HAMMOCK_REHASH_RUN_ROOT:?set HAMMOCK_REHASH_RUN_ROOT to the disposable shared run directory}"
export PYTHONPATH="${HAMMOCK_REHASH_RUN_ROOT}/runtime/site-packages"
export LD_PRELOAD="/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/libstdc++.so.6"
exec /home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python -m hammock "$@"
