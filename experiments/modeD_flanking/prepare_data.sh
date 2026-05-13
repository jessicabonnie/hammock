#!/usr/bin/env bash
# Stage data for modeD_flanking. Idempotent.
#
# What it does:
#   1. Creates /vast symlinks for data/, results/, logs/.
#   2. Symlinks the Maurano DHS experiment dir under data/maurano_link/
#      so Part 1 can read raw Mode D CSVs + FASTAs + bedtools ref directly.
#   3. Creates the synthetic/ subdir for Part 2 (filled in by generate_synthetic.py).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

VAST_ROOT="/vast/blangme2/jbonnie/hammock_claude_experiments/modeD_flanking"
MAURANO_DIR="/home/jbonnie1/interval_sketch/hammock_claude/experiments/maurano_dhs_validation"

mkdir -p "$VAST_ROOT"/{data,results,logs}
for sub in data results logs; do
  if [ ! -e "$sub" ]; then
    ln -sfn "$VAST_ROOT/$sub" "$sub"
  fi
done

# Part 1: link maurano_dhs_validation as a read-only source.
if [ ! -e data/maurano_link ]; then
  ln -sfn "$MAURANO_DIR" data/maurano_link
fi

# Part 2: synthetic FASTA pairs land here.
mkdir -p data/synthetic

echo "Staged:"
echo "  data/maurano_link -> $MAURANO_DIR"
echo "  data/synthetic    -> (empty; run generate_synthetic.py for Part 2)"
