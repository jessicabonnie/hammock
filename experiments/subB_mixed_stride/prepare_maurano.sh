#!/usr/bin/env bash
# Symlink the 20 Maurano fetal-tissue DHS BED files into data/maurano/.
# Source dir is the same one maurano_dhs_validation uses.
#
# Idempotent: re-running just re-creates the symlinks.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

OLD_BEDS_DIR="/home/jbonnie1/interval_sketch/hammock/experiments/dnase1-hypersensitivity/data/maurano.dnaseI"

if [ ! -d "$OLD_BEDS_DIR" ]; then
    echo "ERROR: source BED dir not found: $OLD_BEDS_DIR" >&2
    exit 1
fi

mkdir -p data/maurano

n=0
for bed in "$OLD_BEDS_DIR"/*.bed; do
    [ -e "$bed" ] || continue
    ln -sfn "$bed" "data/maurano/$(basename "$bed")"
    n=$((n + 1))
done

echo "linked $n BED files into data/maurano/"
ls -lh data/maurano/ | head -5
echo "..."
echo "total: $(ls data/maurano/*.bed | wc -l) files"
