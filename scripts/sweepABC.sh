#!/usr/bin/env bash
set -euo pipefail

# Requires: GNU parallel (install on RHEL/CentOS: sudo dnf install -y parallel)
ml parallel
expA=(2.5 3 3.5 4)
# expA=(0.1 0.2 0.5 0.8 1 1.2 1.5 1.8 2 2.5 3 3.5 4)
subB=(0.01 0.05 0.1 0.2 0.4 0.6 0.8 1)
precision=($(seq 24 1 26))
bedfiles=$1

if [ -z "$bedfiles" ]; then
  echo "Usage: $0 <bedfiles>"
  exit 1
fi

for p in "${precision[@]}"; do
  hammock $bedfiles $bedfiles --mode A --precision $p
  hammock $bedfiles $bedfiles --mode B --precision $p
done

#"$(nproc)"

# parallel --will-cite --jobs 28 --halt now,fail=1 --eta \
#   "echo 'Running combo: expA={1} precision={2}' && \
#   hammock $bedfiles $bedfiles --mode C --expA {1} --precision {2}" \
# ::: "${expA[@]}" ::: "${precision[@]}"


parallel --will-cite --jobs 28 --halt now,fail=1 --eta \
  "echo 'Running combo: subB={1} precision={2}' && \
  hammock $bedfiles $bedfiles --mode C --subB {1} --precision {2}" \
::: "${subB[@]}" ::: "${precision[@]}"
