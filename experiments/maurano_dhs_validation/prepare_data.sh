#!/usr/bin/env bash
# Prepare the Maurano DHS corpus for the maurano_dhs_validation experiment.
#
# What it does:
#   1. Ensures data/, results/, logs/ exist as /vast symlinks.
#   2. Symlinks the 20 BED files and 20 hg19 FASTA files from the old
#      hammock experiment dir into data/maurano.dnaseI/ and data/fastas/.
#   3. Writes data/maurano_files.txt, data/maurano_fastas.txt (absolute paths).
#   4. Copies data/maurano_filenames_key.tsv (tissue metadata).
#   5. (Re)builds data/maurano_bedtools_ref.tsv — pairwise bedtools Jaccard
#      ground truth — only if missing.
#
# Idempotent. Safe to re-run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

VAST_ROOT="/vast/blangme2/jbonnie/hammock_claude_experiments/maurano_dhs_validation"
OLD_DATA="/home/jbonnie1/interval_sketch/hammock/experiments/dnase1-hypersensitivity/data"
OLD_BEDS_DIR="$OLD_DATA/maurano.dnaseI"
OLD_FAS_DIR="$OLD_DATA/fastas"
OLD_KEY="$OLD_DATA/maurano_filenames_key.tsv"

mkdir -p "$VAST_ROOT"/{data,results,logs}
for sub in data results logs; do
  if [ ! -e "$sub" ]; then
    ln -sfn "$VAST_ROOT/$sub" "$sub"
  fi
done

mkdir -p data/maurano.dnaseI data/fastas

# Symlink BEDs.
if [ ! -d "$OLD_BEDS_DIR" ]; then
  echo "ERROR: old BED dir not found: $OLD_BEDS_DIR" >&2
  echo "Re-download per the original tutorial (see README) or fix the path." >&2
  exit 1
fi
for bed in "$OLD_BEDS_DIR"/*.bed; do
  base=$(basename "$bed")
  ln -sfn "$bed" "data/maurano.dnaseI/$base"
done

# Symlink FASTAs (already produced by bedtools getfasta against hg19).
if [ -d "$OLD_FAS_DIR" ]; then
  for fa in "$OLD_FAS_DIR"/*.fa; do
    base=$(basename "$fa")
    ln -sfn "$fa" "data/fastas/$base"
  done
else
  echo "WARN: no FASTA dir at $OLD_FAS_DIR; Mode D sweep will fail until you run bedtools getfasta." >&2
fi

# File lists with absolute paths so hammock's filepaths_file works from anywhere.
ls -1 data/maurano.dnaseI/*.bed | xargs -I {} realpath {} > data/maurano_files.txt
if compgen -G "data/fastas/*.fa" > /dev/null; then
  ls -1 data/fastas/*.fa | xargs -I {} realpath {} > data/maurano_fastas.txt
fi

# Tissue metadata.
if [ -f "$OLD_KEY" ]; then
  cp -f "$OLD_KEY" data/maurano_filenames_key.tsv
fi

# Bedtools pairwise reference (ground truth).
REF=data/maurano_bedtools_ref.tsv
if [ ! -s "$REF" ]; then
  if ! command -v bedtools >/dev/null 2>&1; then
    echo "WARN: bedtools not on PATH; skipping reference build. Symlinking old reference if present." >&2
    if [ -f "$OLD_DATA/maurano_bedtools_ref.tsv" ]; then
      ln -sfn "$OLD_DATA/maurano_bedtools_ref.tsv" "$REF"
    fi
  else
    echo "Computing bedtools pairwise Jaccard reference -> $REF ..."
    {
      printf "file1\tfile2\tintersection\tunion\tjaccard\tn_intersections\n"
      while read -r a; do
        while read -r b; do
          ja=$(bedtools jaccard -a "$a" -b "$b" 2>/dev/null | tail -n1)
          inter=$(awk '{print $1}' <<<"$ja")
          un=$(awk '{print $2}' <<<"$ja")
          j=$(awk '{print $3}' <<<"$ja")
          ni=$(awk '{print $4}' <<<"$ja")
          printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$(basename "$a")" "$(basename "$b")" "$inter" "$un" "$j" "$ni"
        done < data/maurano_files.txt
      done < data/maurano_files.txt
    } > "$REF"
  fi
fi

echo "Prepared $(wc -l < data/maurano_files.txt) BEDs and $(wc -l < data/maurano_fastas.txt 2>/dev/null || echo 0) FASTAs."
echo "Bedtools reference: $REF"
