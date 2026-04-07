#!/usr/bin/env Rscript
# =============================================================================
# scripts/exp_b_cluster_plot.R
# Experiment B: Tissue-over-Species Clustering
#
# Reproduces the clustering hierarchy from:
#   Yue et al. (2014) Nature 515:355-364  (Mouse ENCODE)
#   Lin et al. (2014) PNAS 111:17224-17229
#   Roadmap Epigenomics Consortium (2015) Nature 518:317-330
#
# Input:  similarity_matrix.tsv  (all-vs-all pairwise Jaccard similarities)
#         exp_b_metadata.tsv     (sample_id, tissue, species, mark)
# Output: dendrogram.png
#         pca.png
#         cluster_assignments.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggdendro)
  library(cowplot)
  library(factoextra)
  library(RColorBrewer)
})

# ── Snakemake I/O ─────────────────────────────────────────────────────────
mat_file      <- snakemake@input[["matrix"]]
meta_file     <- snakemake@input[["metadata"]]
out_dendro    <- snakemake@output[["dendrogram"]]
out_pca       <- snakemake@output[["pca"]]
out_clusters  <- snakemake@output[["clusters"]]

# ── Load data ─────────────────────────────────────────────────────────────
mat  <- read_tsv(mat_file,  show_col_types = FALSE) %>%
  column_to_rownames("sample_a")
meta <- read_tsv(meta_file, show_col_types = FALSE)

# Convert similarity to distance
dist_mat <- as.dist(1 - as.matrix(mat))

# ── Hierarchical clustering (complete linkage — as in Roadmap 2015) ────────
hc <- hclust(dist_mat, method = "complete")

# Extract cluster assignments at k = number of tissues
tissues <- unique(meta$tissue)
k       <- length(tissues)
clusters <- cutree(hc, k = k)

cluster_df <- tibble(
  sample_id = names(clusters),
  cluster   = clusters
) %>%
  left_join(meta, by = "sample_id")

write_tsv(cluster_df, out_clusters)

# ── Dendrogram ────────────────────────────────────────────────────────────
# Color by species; label by tissue_species
label_df <- meta %>%
  mutate(label = paste0(tissue, " (", str_sub(species, 1, 2), ")"))

dendro_data <- dendro_data(as.dendrogram(hc))

# Add color metadata to labels
label_with_meta <- dendro_data$labels %>%
  left_join(meta, by = c("label" = "sample_id"))

species_colors <- c(human = "#2166ac", mouse = "#d73027")

p_dendro <- ggplot() +
  geom_segment(
    data = dendro_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    size = 0.5
  ) +
  geom_text(
    data = label_with_meta,
    aes(x = x, y = y - 0.01, label = paste0(tissue, "\n(", species, ")"),
        color = species),
    hjust = 1, size = 3, angle = 0
  ) +
  scale_color_manual(values = species_colors) +
  coord_flip() +
  scale_y_reverse() +
  labs(
    title    = "Experiment B: Tissue-over-Species Clustering",
    subtitle = paste0(
      "Minimizer sketch similarity, complete linkage hierarchical clustering\n",
      "Reproduces Yue et al. (2014) and Lin et al. (2014)"
    ),
    x     = NULL,
    y     = "Distance (1 − Jaccard similarity)",
    color = "Species"
  ) +
  theme_minimal(12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(out_dendro, p_dendro, width = 10, height = 7, dpi = 300)

# ── PCA ───────────────────────────────────────────────────────────────────
# PCA on the similarity matrix rows (each sample as a point)
pca_res <- prcomp(as.matrix(mat), scale. = FALSE)

pca_df <- as_tibble(pca_res$x[, 1:3]) %>%
  mutate(sample_id = rownames(pca_res$x)) %>%
  left_join(meta, by = "sample_id")

var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

tissue_shapes <- setNames(
  seq_along(unique(pca_df$tissue)),
  unique(pca_df$tissue)
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                              color = species, shape = tissue, label = tissue)) +
  geom_point(size = 4) +
  geom_text(nudge_y = 0.01, size = 3, show.legend = FALSE) +
  scale_color_manual(values = species_colors) +
  scale_shape_manual(values = tissue_shapes) +
  labs(
    title    = "PCA of Minimizer Sketch Similarity",
    subtitle = "Points colored by species, shaped by tissue",
    x        = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
    y        = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
    color    = "Species",
    shape    = "Tissue"
  ) +
  theme_cowplot(12)

ggsave(out_pca, p_pca, width = 9, height = 7, dpi = 300)

message("Done. Outputs: ", out_dendro, ", ", out_pca, ", ", out_clusters)

