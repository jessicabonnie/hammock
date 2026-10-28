#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
data_dir="${script_dir}/data"
figure_dir="${repo_root}/paper/figures"
key="${data_dir}/maurano_filenames_key.tsv"
work_dir="$(mktemp -d)"
trap 'rm -rf -- "${work_dir}"' EXIT

verify_sha256() {
  local path="$1" expected="$2" observed
  observed="$(sha256sum "${path}" | cut -d' ' -f1)"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "Rehash input checksum mismatch: ${path}" >&2
    echo "expected ${expected}" >&2
    echo "observed ${observed}" >&2
    exit 1
  fi
}

verify_sha256 "${data_dir}/p12_seed00042_k10_w30.csv" \
  5c7014c64295b0dec60387e7f4d78f9d0a9b9e390a1406d7a7c59800fcc50e73
verify_sha256 "${data_dir}/p18_seed00042_k10_w30.csv" \
  de146ace640cfa153e1dc63c4bcb17e48f10597f0b641951d193e14fd6a9e179
verify_sha256 "${data_dir}/p24_seed00042_k10_w30.csv" \
  bd0b42953d0bd44207f6b3f5668dba6cfabd2eca75f24c13193084a134b6ffe2
verify_sha256 "${data_dir}/p18_seed00042_k20_w20.csv" \
  1ff426bf2a092b4cd14de24ccc1bc53f729b111ccea9c5962e3f903e6bcf3434
verify_sha256 "${data_dir}/p18_seed00012_k20_w20.csv" \
  6489cceb0e207f8c6d651b6fd358bdefb2103422afbfd6eb86efc47841f854bc

if [[ "${1:-}" == "--verify-only" ]]; then
  echo "All tissue-figure inputs match the validated rehash checksums."
  exit 0
fi

render_dendrogram() {
  local input="$1" output="$2" axis_break="$3" clusters="$4"
  local pdf="${work_dir}/$(basename "${output}" .png).pdf"
  local png_base="${work_dir}/$(basename "${output}" .png)"
  Rscript "${script_dir}/plot_sequence_tissue_clustering.R" \
    "${input}" "${key}" "${pdf}" jaccard_similarity_ie A \
    "${axis_break}" 6 true 18 "${clusters}" "" none tissue
  pdftoppm -png -r 300 -singlefile "${pdf}" "${png_base}"
  cp "${png_base}.png" "${output}"
}

render_heatmap() {
  local input="$1" output="$2" clusters="$3"
  local pdf="${work_dir}/$(basename "${output}" .png).pdf"
  local png_base="${work_dir}/$(basename "${output}" .png)"
  Rscript "${script_dir}/plot_distance_heatmap.R" \
    "${input}" "${key}" "${pdf}" jaccard_similarity_ie B false true 18 false \
    "${clusters}" tissue
  pdftoppm -png -r 300 -singlefile "${pdf}" "${png_base}"
  cp "${png_base}.png" "${output}"
}

k10_input="${data_dir}/p18_seed00042_k10_w30.csv"
k20_input="${data_dir}/p18_seed00042_k20_w20.csv"

render_dendrogram "${k10_input}" "${figure_dir}/sequence_tissue_clustering.png" 0.05 10
render_heatmap "${k10_input}" "${figure_dir}/sequence_tissue_distance_heatmap.png" 10
render_dendrogram "${k10_input}" "${figure_dir}/sequence_tissue_clustering_8class.png" 0.05 8
render_heatmap "${k10_input}" "${figure_dir}/sequence_tissue_distance_heatmap_8class.png" 8

render_dendrogram "${k20_input}" "${figure_dir}/sequence_tissue_clustering_k20_w20.png" 0.30 10
render_heatmap "${k20_input}" "${figure_dir}/sequence_tissue_distance_heatmap_k20_w20.png" 10
render_dendrogram "${k20_input}" "${figure_dir}/sequence_tissue_clustering_k20_w20_8class.png" 0.30 8
render_heatmap "${k20_input}" "${figure_dir}/sequence_tissue_distance_heatmap_k20_w20_8class.png" 8

Rscript "${script_dir}/combine_precision_dendrograms.R"
Rscript "${script_dir}/combine_precision_dendrograms.R" \
  "${figure_dir}/sequence_tissue_clustering_precision_8class.png" 8
Rscript "${script_dir}/combine_precision_heatmaps.R"
Rscript "${script_dir}/combine_precision_heatmaps.R" \
  "${figure_dir}/sequence_tissue_distance_heatmap_precision_8class.png" 8
Rscript "${script_dir}/combine_k20_w20_cluster_cuts.R"

echo "Rendered rehashed tissue dendrograms and heatmaps in ${figure_dir}"
