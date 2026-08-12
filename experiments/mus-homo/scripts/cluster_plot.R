#!/usr/bin/env Rscript
# =============================================================================
# mus-homo/scripts/cluster_plot.R
# Tissue-over-species clustering on DNase-seq peak FASTAs (human + mouse)
#
# Input:  hammock all-vs-all similarity CSV
#         samples.tsv (with at least sample_id, tissue, species columns)
# Output: dendrogram.png, pca.png, cluster_assignments.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(Cairo)
  library(ggdendro)
  library(cowplot)
})

cairo_png_in <- function(filename, width, height, ...) {
  Cairo::CairoPNG(filename, width = width, height = height,
                  units = "in", res = 300, ...)
}

# ── Snakemake I/O ─────────────────────────────────────────────────────────
mat_file          <- snakemake@input[["matrix"]]
meta_file         <- snakemake@input[["metadata"]]
out_dendro        <- snakemake@output[["dendrogram"]]
out_pca           <- snakemake@output[["pca"]]
out_clusters      <- snakemake@output[["clusters"]]
# workflow/Snakefile:155 passes config.yaml's `primary_sim_col`
# ("jaccard_similarity") as this param -- an inert legacy default now, kept
# untouched deliberately (a static Snakemake `params:` literal has no
# fallback capability of its own). The live column name is resolved below
# from mat_long's actual columns instead: prefer `reg_eq_similarity`, fall
# back to this Snakemake-passed literal for archived pre-rename CSVs.
snakemake_sim_col <- snakemake@params[["sim_col"]]

mat_long <- read_csv(mat_file, show_col_types = FALSE)
meta     <- read_tsv(meta_file, show_col_types = FALSE) %>%
            select(sample_id, tissue, species, ref)

if ("reg_eq_similarity" %in% names(mat_long)) {
  sim_col <- "reg_eq_similarity"
} else {
  message("cluster_plot.R: 'reg_eq_similarity' column not found, falling back to ",
          "Snakemake-passed sim_col '", snakemake_sim_col, "' (archived pre-rename CSV?)")
  sim_col <- snakemake_sim_col
}

mat_long <- mat_long %>%
  mutate(
    sample_a = str_remove(basename(file1), "\\.fa$"),
    sample_b = str_remove(basename(file2), "\\.fa$"),
  )

samples <- unique(c(mat_long$sample_a, mat_long$sample_b))
sim_wide <- mat_long %>%
  select(sample_a, sample_b, similarity = !!sym(sim_col)) %>%
  pivot_wider(names_from = sample_b, values_from = similarity) %>%
  column_to_rownames("sample_a")
sim_wide <- sim_wide[samples, samples]
diag(sim_wide) <- 1.0

dist_mat <- as.dist(1 - as.matrix(sim_wide))
hc <- hclust(dist_mat, method = "complete")

tissues <- unique(meta$tissue)
k       <- length(tissues)
clusters <- cutree(hc, k = k)
cluster_df <- tibble(sample_id = names(clusters), cluster = clusters) %>%
  left_join(meta, by = "sample_id")
write_tsv(cluster_df, out_clusters)

# ── Dendrogram ────────────────────────────────────────────────────────────
dendro_data <- dendro_data(as.dendrogram(hc))
label_with_meta <- dendro_data$labels %>%
  left_join(meta, by = c("label" = "sample_id"))

species_colors <- c(human = "#2166ac", mouse = "#d73027")

p_dendro <- ggplot() +
  geom_segment(data = dendro_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5) +
  geom_text(data = label_with_meta,
            aes(x = x, y = y - 0.01,
                label = paste0(tissue, " (", str_sub(species, 1, 2), ")"),
                color = species),
            hjust = 1, size = 3.5) +
  scale_color_manual(values = species_colors) +
  coord_flip() + scale_y_reverse() +
  labs(title    = "mus-homo: DNase-seq peak-FASTA clustering",
       subtitle = "Minimizer sketch similarity, complete linkage",
       x = NULL, y = "Distance (1 - Jaccard similarity)", color = "Species") +
  theme_minimal(12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(out_dendro, p_dendro, width = 10, height = 7, dpi = 300, device = cairo_png_in)

# ── PCA ───────────────────────────────────────────────────────────────────
pca_res <- prcomp(as.matrix(sim_wide), scale. = FALSE)
pca_df <- as_tibble(pca_res$x[, 1:min(3, ncol(pca_res$x))]) %>%
  mutate(sample_id = rownames(pca_res$x)) %>%
  left_join(meta, by = "sample_id")
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = species, shape = tissue)) +
  geom_point(size = 4) +
  geom_text(aes(label = tissue), nudge_y = 0.02, size = 3, show.legend = FALSE) +
  scale_color_manual(values = species_colors) +
  labs(title    = "mus-homo: PCA of minimizer sketch similarity",
       subtitle = "Points colored by species, shaped by tissue",
       x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
       color = "Species", shape = "Tissue") +
  theme_cowplot(12)

ggsave(out_pca, p_pca, width = 9, height = 7, dpi = 300, device = cairo_png_in)

message("Done. Outputs: ", out_dendro, ", ", out_pca, ", ", out_clusters)
