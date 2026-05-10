#!/usr/bin/env bash
# Pairwise jaccard between two lists of BED files using bedtools.
#
# NOTE: BED files must already be sorted by chr and position. The benchmark
# driver sorts them once before timing starts so sorting time isn't charged
# to bedtools.
#
# Usage: ./bedtools.sh <file1_list> <file2_list> [num_threads]

# Module names on Rockfish (versioned). Fallbacks for other systems.
ml bedtools2 2>/dev/null || ml bedtools 2>/dev/null || true
ml parallel 2>/dev/null || true

set -euo pipefail

file1=$1
file2=$2
NUM_THREADS=${3:-1}

mapfile -t file1_lines < "$file1"
mapfile -t file2_lines < "$file2"

run_jaccard() {
    local jac
    jac=$(bedtools jaccard -a "$1" -b "$2" | tail -n 1 | cut -f 3)
    printf '%s\t%s\t%s\n' "$1" "$2" "$jac"
}
export -f run_jaccard

if command -v parallel &>/dev/null && [ "$NUM_THREADS" -gt 1 ]; then
    parallel --jobs "$NUM_THREADS" run_jaccard ::: "${file1_lines[@]}" ::: "${file2_lines[@]}"
else
    if [ "$NUM_THREADS" -gt 1 ]; then
        echo "Warning: GNU parallel not found, running sequentially." >&2
    fi
    for a in "${file1_lines[@]}"; do
        for b in "${file2_lines[@]}"; do
            run_jaccard "$a" "$b"
        done
    done
fi
