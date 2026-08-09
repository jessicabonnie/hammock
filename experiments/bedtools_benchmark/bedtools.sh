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

# Buffer the result so the pair-count check below can run before anything
# reaches stdout -- a partial stream that the caller has already started
# parsing is exactly the failure this check exists to prevent.
TMP_OUT=$(mktemp)
trap 'rm -f "$TMP_OUT"' EXIT

EXPECTED=$(( ${#file1_lines[@]} * ${#file2_lines[@]} ))

# ONE process per pair, not three.
#
# This used to run an exported bash function doing
#     jac=$(bedtools jaccard -a "$1" -b "$2" | tail -n 1 | cut -f 3)
# which costs a shell + bedtools + tail + cut for every pair. That matters far
# more than it looks, because a pairwise bedtools workflow launches a process
# per pair -- N^2 of them -- and on Rockfish compute nodes process creation is
# capped near 123 exec/s and does NOT scale with cores. Under that cap, exec
# count is the throughput, so trimming 3 execs to 1 is close to a 3x change in
# what the baseline appears to cost.
#
# Measured, 1024 pairs at --jobs 16 on a `parallel` node, against 12.95 ms/pair
# of serial bedtools work:
#
#   fn + bedtools|tail|cut (old)     16.41 s   parallel efficiency 0.81x
#   fn + bedtools|awk                12.23 s                       1.08x
#   --tagstring, bedtools only       7.89 s                        1.68x
#
# 1.68x out of a possible 16x is still poor -- that is the process-creation cap,
# not something this script can fix -- but it is the fairest available shot at
# the baseline, and understating bedtools would overstate our own speedup.
#
# --tagstring prefixes each of bedtools' two output lines with the pair, so the
# per-pair text handling disappears and ONE awk parses the whole stream at the
# end. Fields are: 1=file_a 2=file_b 3=intersection 4=union 5=jaccard
# 6=n_intersections. The `$3+0==$3` test keeps only rows whose intersection
# column is numeric, which drops bedtools' repeated header line.
parse_stream() {
    awk -F'\t' '$3 != "" && $3+0 == $3 { print $1 "\t" $2 "\t" $5 }'
}

if command -v parallel &>/dev/null && [ "$NUM_THREADS" -gt 1 ]; then
    parallel --jobs "$NUM_THREADS" --tagstring '{1}\t{2}' \
        bedtools jaccard -a {1} -b {2} \
        :::: "$file1" :::: "$file2" | parse_stream > "$TMP_OUT"
else
    if [ "$NUM_THREADS" -gt 1 ]; then
        echo "bedtools.sh: GNU parallel not found -- running SEQUENTIALLY." >&2
        echo "  Timings from this run are single-threaded regardless of the" >&2
        echo "  requested thread count; do not report them as parallel." >&2
    fi
    for a in "${file1_lines[@]}"; do
        for b in "${file2_lines[@]}"; do
            bedtools jaccard -a "$a" -b "$b" \
                | awk -F'\t' -v q="$a" -v r="$b" 'NR==2 {print q "\t" r "\t" $3}'
        done
    done > "$TMP_OUT"
fi

# A dropped pair used to be invisible: the caller would just see fewer rows and
# average over whatever survived. Failing loudly here is the difference between
# a wrong MAE and a stopped run.
GOT=$(wc -l < "$TMP_OUT")
if [ "$GOT" -ne "$EXPECTED" ]; then
    echo "bedtools.sh: expected $EXPECTED pairs, got $GOT." >&2
    echo "  Some bedtools invocations produced no parseable jaccard." >&2
    exit 1
fi
cat "$TMP_OUT"
