#!/usr/bin/env bash
# Pairwise jaccard between two lists of BED files using bedtools.
#
# NOTE: BED files must already be sorted by chr and position. The benchmark
# driver sorts them once before timing starts so sorting time isn't charged
# to bedtools.
#
# Usage: ./bedtools.sh <file1_list> <file2_list> [num_threads]

# Pinned, and loudly. This used to read
#     ml bedtools2 2>/dev/null || ml bedtools 2>/dev/null || true
# which resolved to bedtools2/2.27.1 -- a 2017 build -- while
# estimator_compare.py hardcoded 2.30.0, so two parts of the same experiment
# were quietly comparing against different bedtools. Worse, the trailing
# `|| true` meant that if no module loaded at all the script continued and
# failed later, or silently used whatever `bedtools` happened to be on PATH.
#
# 2.30.0 is safe to pin: on three Maurano pairs it returns output byte-identical
# to 2.27.1 on all four columns (intersection, union, jaccard, n_intersections),
# so no archived jaccard changes and the accuracy gates are unaffected.
#
# HAMMOCK_BEDTOOLS_MODULE overrides the pin for a machine with different module
# names; an empty value skips the load and uses PATH.
BEDTOOLS_MODULE="${HAMMOCK_BEDTOOLS_MODULE-bedtools/2.30.0}"
if [ -n "$BEDTOOLS_MODULE" ]; then
    ml "$BEDTOOLS_MODULE" 2>/dev/null || {
        echo "bedtools.sh: could not load module '$BEDTOOLS_MODULE'." >&2
        echo "  Set HAMMOCK_BEDTOOLS_MODULE to a module that exists, or to the" >&2
        echo "  empty string to use whatever bedtools is on PATH." >&2
        exit 1
    }
fi
ml parallel 2>/dev/null || true

command -v bedtools >/dev/null 2>&1 || {
    echo "bedtools.sh: no bedtools on PATH after module load" >&2
    exit 1
}
# To stderr, so it lands in the job log without touching the stdout the caller
# parses per-pair jaccards out of. Which binary produced a timing is not
# recoverable after the fact otherwise.
echo "bedtools.sh: $(bedtools --version) at $(command -v bedtools)" >&2

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
