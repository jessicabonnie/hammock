#!/usr/bin/env bash
set -euo pipefail

# Requires: GNU parallel (install on RHEL/CentOS: sudo dnf install -y parallel)
ml parallel
klen=(8 10 15 20 25 50 100)
window=(8 10 20 30 50 100 200 500)
precision=($(seq 21 1 24))

fastafiles=$1

if [ -z "$fastafiles" ]; then
  echo "Usage: $0 <fastafiles>"
  exit 1
fi

#"$(nproc)"
for p in "${precision[@]}"; do
  echo "Running precision $p..."
  parallel --will-cite --jobs 28  --halt now,fail=1 --eta \
    "echo 'Running combo: klen={1} window={2} precision=$p' && \
    hammock $fastafiles $fastafiles --mode D --kmer_size {1} --window_size {2} --precision $p" \
  ::: "${klen[@]}" ::: "${window[@]}"
done
