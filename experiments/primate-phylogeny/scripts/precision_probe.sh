#!/usr/bin/env bash
# =============================================================================
# precision_probe.sh — sweep HLL precision at the two best (k, w) cells from
# Phase 1 to see whether the cross-species similarity dynamic range widens.
#
# Phase 1 finding: at p=24, cross-species Jaccard sits in 0.64–0.73 — only
# 0.09 of spread to encode ~180 Myr of mammalian divergence. This probe tests
# whether p=26 or p=28 produces a wider, more topology-resolving spread.
#
# Cells probed: k∈{8, 10}, w=10. Precisions: 24 (baseline, reused from
# Phase 1), 26, 28. Outputs go under OUTDIR/precision_probe/.
# =============================================================================
set -euo pipefail

ROOT="/home/jbonnie1/interval_sketch/hammock_claude/experiments/primate-phylogeny"
OUTDIR="/vast/blangme2/jbonnie/hammock/primate-phylogeny/results/H3K4me3/precision_probe"
FASTA_LIST="/vast/blangme2/jbonnie/hammock/primate-phylogeny/results/H3K4me3/fasta_lists/all.txt"
HAMMOCK="hammock"

mkdir -p "$OUTDIR"

declare -a runs=(
  "8  10 26"
  "8  10 28"
  "10 10 26"
  "10 10 28"
)

# Submit each as its own sbatch — they parallelize on the cluster.
jobids=()
for spec in "${runs[@]}"; do
  read -r k w p <<<"$spec"
  cell_dir="$OUTDIR/k${k}_w${w}_p${p}"
  mkdir -p "$cell_dir"
  log="$cell_dir/hammock.log"
  # Memory scales ~4× per precision step (2^p registers per sketch × 7 samples).
  # p=26 → ~512 MB, p=28 → ~2 GB. Give 16 GB headroom.
  jid=$(sbatch --parsable \
    --partition=parallel \
    --account=blangme2_scal \
    --mem=16000M \
    --time=120 \
    --cpus-per-task=8 \
    --job-name="precprobe_k${k}_w${w}_p${p}" \
    --output="$cell_dir/slurm.out" \
    --error="$cell_dir/slurm.err" \
    --wrap="set -euo pipefail; ${HAMMOCK} ${FASTA_LIST} ${FASTA_LIST} --mode D --outprefix ${cell_dir}/probe -k ${k} -w ${w} --precision ${p} --full-paths 2>${log}")
  echo "submitted k=${k} w=${w} p=${p} → SLURM ${jid}"
  jobids+=("$jid")
done

echo
echo "Waiting on jobs: ${jobids[*]}"
# Poll squeue every 30 s until all 4 jobs leave the queue.
while true; do
  remaining=$(squeue -h -j "$(IFS=,; echo "${jobids[*]}")" 2>/dev/null | wc -l)
  echo "  $(date +%H:%M:%S) — $remaining of ${#jobids[@]} still in queue"
  (( remaining == 0 )) && break
  sleep 30
done

echo
echo "=== cross-species similarity spread per cell ==="
printf "%-22s %8s %8s %8s %8s %8s\n" "cell" "min" "max" "spread" "hsa_rh" "Mmus_mD"
for spec in "${runs[@]}"; do
  read -r k w p <<<"$spec"
  cell_dir="$OUTDIR/k${k}_w${w}_p${p}"
  csv=$(ls "${cell_dir}"/probe*.csv 2>/dev/null | head -1 || true)
  if [[ -z "$csv" || ! -s "$csv" ]]; then
    printf "%-22s  (no csv)\n" "k${k}_w${w}_p${p}"
    continue
  fi
  awk -F, -v label="k${k}_w${w}_p${p}" '
    NR==1 {
      for (i=1;i<=NF;i++) h[$i]=i
      jcol=h["jaccard_similarity"]
      next
    }
    {
      a=$h["file1"]; b=$h["file2"]
      sub(/.*\//,"",a); sub(/.*\//,"",b)
      sub(/\.fa$/,"",a); sub(/\.fa$/,"",b)
      if (a==b) next
      v=$jcol+0
      if (v < mn || mn=="") mn=v
      if (v > mx) mx=v
      key=a"|"b
      val[key]=v
    }
    END {
      hsa_rh = val["hsa|rheMac"]; if (hsa_rh=="") hsa_rh = val["rheMac|hsa"]
      mm_md  = val["Mmus|monDom"]; if (mm_md=="") mm_md = val["monDom|Mmus"]
      printf "%-22s %8.4f %8.4f %8.4f %8.4f %8.4f\n", label, mn, mx, mx-mn, hsa_rh, mm_md
    }' "$csv"
done

echo
echo "Baseline (p=24, from Phase 1):"
for kw in k8_w10 k10_w10; do
  csv="/vast/blangme2/jbonnie/hammock/primate-phylogeny/results/H3K4me3/${kw}/phylo_mnmzr_p24_jaccD_${kw}.csv"
  awk -F, -v label="${kw}_p24" '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; jcol=h["jaccard_similarity"]; next }
    {
      a=$h["file1"]; b=$h["file2"]
      sub(/.*\//,"",a); sub(/.*\//,"",b)
      sub(/\.fa$/,"",a); sub(/\.fa$/,"",b)
      if (a==b) next
      v=$jcol+0
      if (v < mn || mn=="") mn=v
      if (v > mx) mx=v
      key=a"|"b
      val[key]=v
    }
    END {
      hsa_rh = val["hsa|rheMac"]; if (hsa_rh=="") hsa_rh = val["rheMac|hsa"]
      mm_md  = val["Mmus|monDom"]; if (mm_md=="") mm_md = val["monDom|Mmus"]
      printf "%-22s %8.4f %8.4f %8.4f %8.4f %8.4f\n", label, mn, mx, mx-mn, hsa_rh, mm_md
    }' "$csv"
done
