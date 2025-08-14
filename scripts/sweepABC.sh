#!/usr/bin/env bash
set -euo pipefail

# Requires: GNU parallel (install on RHEL/CentOS: sudo dnf install -y parallel)
ml parallel
expA=(0.1 0.2 0.5 0.8 1 1.2 1.5 1.8 2)
precision=($(seq 20 1 24))
bedfiles=$1

if [ -z "$bedfiles" ]; then
  echo "Usage: $0 <bedfiles>"
  exit 1
fi

hammock $bedfiles $bedfiles --mode A --precision 22
hammock $bedfiles $bedfiles --mode B --precision 22

parallel --will-cite --jobs "$(nproc)" --halt now,fail=1 --eta \
  "echo 'Running combo: expA={1} precision={2}' && \
  hammock $bedfiles $bedfiles --mode C --expA {1} --precision {2}" \
::: "${expA[@]}" ::: "${precision[@]}"
