#!/usr/bin/env Rscript
# =============================================================================
# OBSOLETE as of hammock v0.6.0 -- DOES NOT RUN on current output.
#
# This compares the two Mode D Jaccard columns, but `jaccard_similarity_with_ends`
# no longer exists (CLAUDE.md divergence #8, docs/mode-d-ends-removal.md), so
# cluster_one() below will fail on any CSV written by v0.6.0+. It still works on
# the archived CSVs under results/, which retain the old 20-column schema.
# Kept as a record of the question; it was never run (no column_comparison.tsv
# on disk). Delete once the archived comparison is no longer of interest.
# =============================================================================
# Recompute clustering ARI for both Mode D similarity columns
# (jaccard_similarity vs jaccard_similarity_with_ends) across all (k, w) cells.
#
# Mirrors scripts/cluster_plot.R exactly: complete-linkage hclust on
# 1 - similarity, cut to k = 5 clusters, ARI vs tissue and species labels.
#
# Output: TSV at results/column_comparison.tsv
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

results_dir <- "/vast/blangme2/jbonnie/hammock/mus-homo/results"
meta_file <- "/home/jbonnie1/interval_sketch/hammock_claude/experiments/mus-homo/config/samples.tsv"
out_tsv <- file.path(results_dir, "column_comparison.tsv")

meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  select(sample_id, tissue, species)

# pair-counting ARI (no extra package dependency)
adjusted_rand_index <- function(a, b) {
  stopifnot(length(a) == length(b))
  ct <- table(a, b)
  n <- sum(ct)
  if (n < 2) return(0)
  sum_comb <- sum(choose(ct, 2))
  sum_a <- sum(choose(rowSums(ct), 2))
  sum_b <- sum(choose(colSums(ct), 2))
  total <- choose(n, 2)
  expected <- sum_a * sum_b / total
  max_idx <- 0.5 * (sum_a + sum_b)
  if (max_idx == expected) return(0)
  (sum_comb - expected) / (max_idx - expected)
}

cluster_one <- function(csv_path, sim_col) {
  mat_long <- read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(sample_a = str_remove(basename(file1), "\\.fa$"),
           sample_b = str_remove(basename(file2), "\\.fa$"))
  samples <- unique(c(mat_long$sample_a, mat_long$sample_b))
  sim_wide <- mat_long %>%
    select(sample_a, sample_b, similarity = !!sym(sim_col)) %>%
    pivot_wider(names_from = sample_b, values_from = similarity) %>%
    column_to_rownames("sample_a")
  sim_wide <- sim_wide[samples, samples]
  diag(sim_wide) <- 1.0
  dist_mat <- as.dist(1 - as.matrix(sim_wide))
  hc <- hclust(dist_mat, method = "complete")
  k <- length(unique(meta$tissue))
  clusters <- cutree(hc, k = k)
  cl_df <- tibble(sample_id = names(clusters), cluster = clusters) %>%
    left_join(meta, by = "sample_id")
  list(
    ari_tissue = adjusted_rand_index(cl_df$tissue, cl_df$cluster),
    ari_species = adjusted_rand_index(cl_df$species, cl_df$cluster),
    clusters = cl_df
  )
}

kw_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
kw_dirs <- kw_dirs[str_detect(kw_dirs, "^k\\d+_w\\d+$")]

rows <- list()
for (kw in sort(kw_dirs)) {
  parts <- str_match(kw, "^k(\\d+)_w(\\d+)$")
  k <- as.integer(parts[2]); w <- as.integer(parts[3])
  csv_files <- list.files(file.path(results_dir, kw), pattern = "\\.csv$",
                          full.names = TRUE)
  stopifnot(length(csv_files) == 1)
  res_no <- cluster_one(csv_files[1], "jaccard_similarity")
  res_we <- cluster_one(csv_files[1], "jaccard_similarity_with_ends")
  rows[[length(rows) + 1]] <- tibble(
    k = k, w = w,
    ari_tissue_no_ends     = res_no$ari_tissue,
    ari_species_no_ends    = res_no$ari_species,
    tissue_dominant_no_ends = as.integer(res_no$ari_tissue > res_no$ari_species),
    ari_tissue_with_ends   = res_we$ari_tissue,
    ari_species_with_ends  = res_we$ari_species,
    tissue_dominant_with_ends = as.integer(res_we$ari_tissue > res_we$ari_species),
    delta_ari_tissue       = res_no$ari_tissue - res_we$ari_tissue
  )
}
out <- bind_rows(rows) %>% arrange(k, w) %>%
  mutate(across(where(is.numeric) & !c(k, w, starts_with("tissue_dominant")),
                ~ round(.x, 4)))
write_tsv(out, out_tsv)
cat(sprintf("Wrote %s\n\n", out_tsv))

# ── Summary ───────────────────────────────────────────────────────────────
summarize_col <- function(df, suffix) {
  tdc <- df[[paste0("tissue_dominant_", suffix)]]
  at  <- df[[paste0("ari_tissue_",  suffix)]]
  as_ <- df[[paste0("ari_species_", suffix)]]
  list(td = sum(tdc),
       max_t = max(at), max_s = max(as_),
       above_03 = sum(at > 0.3))
}
cat("=== SUMMARY ===\n")
cat(sprintf("%-32s %-7s %-10s %-10s %s\n",
            "column", "td/20", "maxARI_t", "maxARI_s", "ARI_t>0.3"))
for (suf in c("no_ends" = "jaccard_similarity",
              "with_ends" = "jaccard_similarity_with_ends")) {
  k_ <- names(which(c("no_ends" = "jaccard_similarity",
                      "with_ends" = "jaccard_similarity_with_ends") == suf))
  s <- summarize_col(out, k_)
  cat(sprintf("%-32s %4d/20 %+10.4f %+10.4f  %3d/20\n",
              suf, s$td, s$max_t, s$max_s, s$above_03))
}

cat("\n=== Tissue-dominant cells (ARI_t > ARI_s), jaccard_similarity ===\n")
td_no <- out %>% filter(tissue_dominant_no_ends == 1L) %>%
  arrange(desc(ari_tissue_no_ends))
print(td_no %>% select(k, w, ari_tissue_no_ends, ari_species_no_ends),
      n = Inf)

cat("\n=== Tissue-dominant cells, jaccard_similarity_with_ends ===\n")
td_we <- out %>% filter(tissue_dominant_with_ends == 1L) %>%
  arrange(desc(ari_tissue_with_ends))
print(td_we %>% select(k, w, ari_tissue_with_ends, ari_species_with_ends),
      n = Inf)

cat("\n=== Cells where jaccard_similarity beats _with_ends on tissue ARI by >0.05 ===\n")
beats <- out %>% filter(delta_ari_tissue > 0.05) %>%
  arrange(desc(delta_ari_tissue))
print(beats %>% select(k, w, ari_tissue_no_ends, ari_tissue_with_ends,
                       delta_ari_tissue, ari_species_no_ends,
                       ari_species_with_ends),
      n = Inf)

cat(sprintf("\nVerification (with_ends): species ARI range = %.4f .. %.4f\n",
            min(out$ari_species_with_ends), max(out$ari_species_with_ends)))
