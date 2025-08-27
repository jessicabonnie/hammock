#!/usr/bin/env bash
set -euo pipefail

# Requires: GNU parallel (install on RHEL/CentOS: sudo dnf install -y parallel)
ml parallel

# Memory considerations:
# - Precision 12: 2^12 = 4,096 registers = 4KB
# - Precision 14: 2^14 = 16,384 registers = 16KB  
# - Precision 16: 2^16 = 65,536 registers = 64KB
# - Precision 18: 2^18 = 262,144 registers = 256KB
# - Precision 20: 2^20 = 1,048,576 registers = 1MB
# - Precision 24: 2^24 = 16,777,216 registers = 16MB (included, manageable)
# - Precision 26: 2^26 = 67,108,864 registers = 64MB (excluded, causes memory errors)
expA=(0.5 1 1.5 2 2.5 3 3.5 4)
# expA=(0.1 0.2 0.5 0.8 1 1.2 1.5 1.8 2 2.5 3 3.5 4)
subB=(0.01 0.05 0.1 0.2 0.4 0.6 0.8 1)
# Precision range: includes 24 but avoids 26 to prevent memory errors
# precision=($(seq 24 1 26))  # Original: 26 causes memory errors (64MB per sketch)
precision=($(seq 20 2 24))     # New: includes 24, avoids 26, prevents memory issues
bedfiles=$1

if [ -z "$bedfiles" ]; then
  echo "Usage: $0 <bedfiles>"
  exit 1
fi

#for p in "${precision[@]}"; do
  hammock $bedfiles $bedfiles --mode A --precision 24
 hammock $bedfiles $bedfiles --mode B --precision 24
#done

#"$(nproc)"
parallel --will-cite --jobs 28 --halt now,fail=1 --eta \
  "echo 'Running combo: expA={1} precision=24' && \
  hammock $bedfiles $bedfiles --mode C --expA {1} --precision 24" \
::: "${expA[@]}" 


# parallel --will-cite --jobs 28 --halt now,fail=1 --eta \
#   "echo 'Running combo: expA={1} precision={2}' && \
#   hammock $bedfiles $bedfiles --mode C --expA {1} --precision {2}" \
# ::: "${expA[@]}" ::: "${precision[@]}"


# Using full path to hammock to ensure we use the development version, not system-installed version
# Reduced jobs from 28 to 16 to prevent memory issues with high precision sketches
#parallel --will-cite --jobs 16 --halt now,fail=1 --eta \
#  "echo 'Running combo: subB={1} precision={2}' && \
#  hammock $bedfiles $bedfiles --mode C --subB {1} --precision {2}" \
#::: "${subB[@]}" ::: "${precision[@]}"

parallel --will-cite --jobs 16 --halt now,fail=1 --eta \
 "echo 'Running combo: subB={1} precision=24' && \
 hammock $bedfiles $bedfiles --mode C --subB {1} --precision 24" \
::: "${subB[@]}" 