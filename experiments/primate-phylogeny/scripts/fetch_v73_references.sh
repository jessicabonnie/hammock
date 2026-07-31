#!/usr/bin/env bash
# =============================================================================
# fetch_v73_references.sh — pull Ensembl release-73 (Sept 2013) toplevel DNA
# FASTAs for the species needed by primate-phylogeny Phase 1+.
#
# Villar 2015 peaks (E-MTAB-2633) were called against these exact assemblies,
# so chromosome names and coordinates line up without any liftOver.
#
# Output: one decompressed .fa per species under $DEST_DIR. .fai is generated
# by bedtools getfasta on first use.
# =============================================================================
set -euo pipefail

DEST_DIR="${DEST_DIR:-/vast/blangme2/jbonnie/hammock/primate-phylogeny/references}"
BASE="http://ftp.ensembl.org/pub/release-73/fasta"

mkdir -p "$DEST_DIR"

# Format: ucsc_code  species_lc  AssemblyName
declare -a refs=(
  "rheMac  macaca_mulatta        MMUL_1"
  "calJac  callithrix_jacchus    C_jacchus3.2.1"
  "canFam  canis_familiaris      CanFam3.1"
  "bosTau  bos_taurus            UMD3.1"
  "monDom  monodelphis_domestica BROADO5"
)

fetch_one() {
  local code="$1" species_lc="$2" asm="$3"
  local species_cap="${species_lc^}"
  local fname="${species_cap}.${asm}.73.dna.toplevel.fa.gz"
  local url="${BASE}/${species_lc}/dna/${fname}"
  local out_gz="${DEST_DIR}/${code}.fa.gz"
  local out_fa="${DEST_DIR}/${code}.fa"

  if [[ -s "$out_fa" ]]; then
    echo "[$code] already present ($(du -h "$out_fa" | cut -f1)), skipping"
    return 0
  fi

  echo "[$code] $url"
  curl -fsSL "$url" -o "$out_gz"
  gunzip -f "$out_gz"
  echo "[$code] done — $(du -h "$out_fa" | cut -f1)"
}

export -f fetch_one
export DEST_DIR BASE

# Run all 5 downloads in parallel (background + wait).
pids=()
for spec in "${refs[@]}"; do
  read -r code species_lc asm <<<"$spec"
  fetch_one "$code" "$species_lc" "$asm" &
  pids+=($!)
done

fail=0
for pid in "${pids[@]}"; do
  wait "$pid" || fail=1
done

if (( fail )); then
  echo "FAIL: at least one download failed" >&2
  exit 1
fi
echo "All references staged under: $DEST_DIR"
ls -lh "$DEST_DIR"
