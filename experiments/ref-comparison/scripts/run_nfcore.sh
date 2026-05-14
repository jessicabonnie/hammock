#!/usr/bin/env bash
# =============================================================================
# run_nfcore.sh — Launch nf-core/chipseq for claude-ref-comparison
# =============================================================================
#
# Runs nf-core/chipseq for each genome × peak type (8 runs total):
#   GRCh38, GRCh37, CHM13 (human) × {broad, narrow} + mm10 (mouse) × {broad, narrow}
#
# Alignment steps are shared via -resume caching — only MACS3 onward differs
# between broad and narrow runs for the same genome.
#
# After all runs complete, run the downstream Snakemake pipeline:
#   snakemake --profile workflow/slurm_profile/
#
# Prerequisites:
#   - nextflow installed and on PATH
#   - Singularity available on cluster nodes
#   - nf-core/chipseq pulled: nextflow pull nf-core/chipseq -r 2.1.0
#
# =============================================================================

set -euo pipefail

# ── Argument parsing ───────────────────────────────────────────────────────────
# Usage: bash run_nfcore.sh [REF ...]
#   REF can be: GRCh38, GRCh37, CHM13, mm10, or "all" (default if no args given)
# Examples:
#   bash run_nfcore.sh CHM13          # only CHM13 broad + narrow
#   bash run_nfcore.sh GRCh38 mm10    # GRCh38 and mm10, broad + narrow each
#   bash run_nfcore.sh all            # all 8 runs (same as no arguments)

REQUESTED_REFS=("${@:-all}")

should_run() {
  local ref="$1"
  for r in "${REQUESTED_REFS[@]}"; do
    if [[ "${r}" == "all" || "${r}" == "${ref}" ]]; then
      return 0
    fi
  done
  return 1
}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXPERIMENT_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="${EXPERIMENT_DIR}/config"
NFCORE_OUT="/vast/blangme2/jbonnie/hammock/claude-ref-comparison/nfcore"
export NXF_SINGULARITY_CACHEDIR="/vast/blangme2/jbonnie/hammock/singularity_cache"

NF_CONFIG="${EXPERIMENT_DIR}/nextflow.config"

# Local references (no-alt analysis sets)
GRCH38_FASTA="/data/blangme2/jessica/mus_homo/references/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
GRCH37_FASTA="/data/blangme2/jessica/mus_homo/references/hg19/hg19.fa"
CHM13_FASTA="/data/blangme2/fasta/chm13_v2.0/chm13v2.0.fasta"
MM10_FASTA="/data/blangme2/jessica/mus_homo/references/mm10/mm10_no_alt_analysis_set_ENCODE.fasta"

# GTF annotations
GRCH38_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
GRCH37_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz"
CHM13_GTF="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz"
MM10_GTF="/data/blangme2/jessica/mus_homo/references/mm10/mm10.2021-04-23.ncbiRefSeq.gtf"

# Helper: run nf-core/chipseq for a given genome and peak type
run_chipseq() {
  local label="$1" samplesheet="$2" outdir="$3" fasta="$4" gtf="$5" gsize="$6" narrow="$7"

  local peak_args=""
  if [[ "${narrow}" == "true" ]]; then
    peak_args="--narrow_peak"
  fi

  echo "[${label}]"
  nextflow run nf-core/chipseq -r 2.1.0 \
    -c "${NF_CONFIG}" \
    --input "${samplesheet}" \
    --outdir "${outdir}" \
    --fasta "${fasta}" \
    --gtf "${gtf}" \
    --macs_gsize "${gsize}" \
    ${peak_args} \
    -resume
  echo ""
}

echo "=== nf-core/chipseq runs for claude-ref-comparison ==="
echo "Output root: ${NFCORE_OUT}"
echo "Requested references: ${REQUESTED_REFS[*]}"
echo "Alignment steps are cached — only MACS3 onward re-runs for the second peak type."
echo ""

# ── GRCh38 ──────────────────────────────────────────────────────────────────
if should_run GRCh38; then
  run_chipseq "Human GRCh38 broad" \
    "${CONFIG_DIR}/nfcore_samplesheet_human.csv" \
    "${NFCORE_OUT}/GRCh38" \
    "${GRCH38_FASTA}" "${GRCH38_GTF}" 2700000000 false

  run_chipseq "Human GRCh38 narrow" \
    "${CONFIG_DIR}/nfcore_samplesheet_human.csv" \
    "${NFCORE_OUT}/GRCh38_narrow" \
    "${GRCH38_FASTA}" "${GRCH38_GTF}" 2700000000 true
fi

# ── GRCh37 ──────────────────────────────────────────────────────────────────
if should_run GRCh37; then
  run_chipseq "Human GRCh37 broad" \
    "${CONFIG_DIR}/nfcore_samplesheet_human.csv" \
    "${NFCORE_OUT}/GRCh37" \
    "${GRCH37_FASTA}" "${GRCH37_GTF}" 2700000000 false

  run_chipseq "Human GRCh37 narrow" \
    "${CONFIG_DIR}/nfcore_samplesheet_human.csv" \
    "${NFCORE_OUT}/GRCh37_narrow" \
    "${GRCH37_FASTA}" "${GRCH37_GTF}" 2700000000 true
fi

# ── CHM13 (T2T v2.0) ───────────────────────────────────────────────────────
if should_run CHM13; then
  run_chipseq "Human CHM13 broad" \
    "${CONFIG_DIR}/nfcore_samplesheet_human.csv" \
    "${NFCORE_OUT}/CHM13" \
    "${CHM13_FASTA}" "${CHM13_GTF}" 3055000000 false

  run_chipseq "Human CHM13 narrow" \
    "${CONFIG_DIR}/nfcore_samplesheet_human.csv" \
    "${NFCORE_OUT}/CHM13_narrow" \
    "${CHM13_FASTA}" "${CHM13_GTF}" 3055000000 true
fi

# ── mm10 ────────────────────────────────────────────────────────────────────
if should_run mm10; then
  run_chipseq "Mouse mm10 broad" \
    "${CONFIG_DIR}/nfcore_samplesheet_mouse.csv" \
    "${NFCORE_OUT}/mm10" \
    "${MM10_FASTA}" "${MM10_GTF}" 1870000000 false

  run_chipseq "Mouse mm10 narrow" \
    "${CONFIG_DIR}/nfcore_samplesheet_mouse.csv" \
    "${NFCORE_OUT}/mm10_narrow" \
    "${MM10_FASTA}" "${MM10_GTF}" 1870000000 true
fi

echo "=== All nf-core runs complete ==="
echo "Peak files:"
echo "  Broad:  ${NFCORE_OUT}/<ref>/bwa/merged_library/macs3/broad_peak/*_peaks.broadPeak"
echo "  Narrow: ${NFCORE_OUT}/<ref>_narrow/bwa/merged_library/macs3/narrow_peak/*_peaks.narrowPeak"
echo ""
echo "Next: run the downstream Snakemake pipeline:"
echo "  cd ${EXPERIMENT_DIR}"
echo "  snakemake --profile workflow/slurm_profile/"
