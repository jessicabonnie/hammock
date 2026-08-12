#!/usr/bin/env Rscript
# =============================================================================
# primate-phylogeny/scripts/build_phylogeny.R
# Neighbor-joining trees from hammock peak-FASTA similarity matrix.
#
# Mode D emits 7 similarity metrics on the minimizer HLL (jaccard,
# jaccard_ie, containment_AB, containment_BA, cosketch_geom, cosketch_arith,
# cosketch_max). Containment / cosketch isolate `|A ∩ B| / |A|` from
# union-size penalty, so they typically widen the cross-species margin
# between phylogenetically close pairs and long-branch artifact pairs.
#
# This script builds one NJ tree per metric for a small curated set:
# jaccard + cosketch_max + cosketch_geom = 3 trees per (k, w) cell.
# Containment_AB/BA are asymmetric and not used directly for NJ; their
# information lives in cosketch_max (the symmetric upper bound).
#
# Until hammock v0.6.0 this also built three `_with_ends` twins; that column
# family was removed (CLAUDE.md divergence #8), so `*_we.*` outputs on /vast
# are stale artifacts of earlier runs.
#
# Input:  hammock all-vs-all similarity CSV (Mode D output, 17 cols)
#         samples.tsv (with species_code, common_name, scientific_name, clade)
# Output: nj_tree_{jacc,cosmax,cosgeom}_mz.{newick,png}
#         dist_matrix_{jacc,cosmax,cosgeom}_mz.tsv
#         metric_spreads.tsv  (per-column spreads + key-pair margins)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(Cairo)
  library(ape)
})

cairo_png_in <- function(filename, width, height, ...) {
  Cairo::CairoPNG(filename, width = width, height = height,
                  units = "in", res = 300, ...)
}

# ── Snakemake I/O ─────────────────────────────────────────────────────────
mat_file       <- snakemake@input[["matrix"]]
meta_file      <- snakemake@input[["metadata"]]
out_spread     <- snakemake@output[["spread_tsv"]]
out_newick_jmz <- snakemake@output[["newick_jacc_mz"]]
out_newick_cmz <- snakemake@output[["newick_cosmax_mz"]]
out_newick_gmz <- snakemake@output[["newick_cosgeom_mz"]]
out_png_jmz    <- snakemake@output[["tree_png_jacc_mz"]]
out_png_cmz    <- snakemake@output[["tree_png_cosmax_mz"]]
out_png_gmz    <- snakemake@output[["tree_png_cosgeom_mz"]]
out_dist_jmz   <- snakemake@output[["dist_jacc_mz"]]
out_dist_cmz   <- snakemake@output[["dist_cosmax_mz"]]
out_dist_gmz   <- snakemake@output[["dist_cosgeom_mz"]]

mat_long <- read_csv(mat_file, show_col_types = FALSE)
meta     <- read_tsv(meta_file, show_col_types = FALSE) %>%
            filter(!is.na(species_code) & species_code != "") %>%
            select(species_code, common_name, scientific_name, clade)

mat_long <- mat_long %>%
  mutate(
    sample_a = str_remove(basename(file1), "\\.fa$"),
    sample_b = str_remove(basename(file2), "\\.fa$"),
  )

samples <- unique(c(mat_long$sample_a, mat_long$sample_b))

clade_colors <- c(
  primate        = "#1f77b4",
  scandentia     = "#17becf",
  rodent         = "#ff7f0e",
  lagomorph      = "#bcbd22",
  carnivore      = "#2ca02c",
  cetartiodactyl = "#d62728",
  marsupial      = "#9467bd"
)

# Register-equality column: prefer the new name, fall back to the legacy
# `jaccard_similarity` name for archived pre-rename CSVs. Resolved once
# against the actually-loaded data (not hardcoded at any call site below),
# and logged at most once per run even though it feeds two consumers
# (the metric-spread table and build_and_plot()'s tree call).
.reg_eq_fallback_logged <- FALSE
resolve_reg_eq_col <- function(df) {
  if ("reg_eq_similarity" %in% names(df)) return("reg_eq_similarity")
  if (!.reg_eq_fallback_logged) {
    message("build_phylogeny.R: 'reg_eq_similarity' column not found, ",
            "falling back to legacy 'jaccard_similarity' (archived pre-rename CSV?)")
    .reg_eq_fallback_logged <<- TRUE
  }
  "jaccard_similarity"
}
reg_eq_col <- resolve_reg_eq_col(mat_long)

# ── Metric-spread summary (all hammock similarity columns) ────────────────
all_metrics <- c(
  reg_eq_col, "jaccard_similarity_ie",
  "containment_AB", "containment_BA",
  "cosketch_geom", "cosketch_arith", "cosketch_max"
)
present <- intersect(all_metrics, names(mat_long))

cross_only <- mat_long %>% filter(sample_a != sample_b)
get_pair_val <- function(col, a, b) {
  v <- cross_only %>% filter(sample_a == a, sample_b == b) %>% pull(!!sym(col))
  if (length(v) == 0) NA_real_ else v[1]
}

spread_rows <- lapply(present, function(col) {
  vals <- cross_only[[col]]
  # Canonical label regardless of which name was actually read, so
  # metric_spreads.tsv reads `reg_eq_similarity` in both the prefers and
  # fallback branches.
  metric_label <- if (col == reg_eq_col) "reg_eq_similarity" else col
  tibble(
    metric            = metric_label,
    min_cross         = min(vals, na.rm = TRUE),
    max_cross         = max(vals, na.rm = TRUE),
    spread            = max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE),
    hsa_rheMac        = get_pair_val(col, "hsa", "rheMac"),
    Mmus_monDom       = get_pair_val(col, "Mmus", "monDom"),
    primate_vs_LBA    = get_pair_val(col, "hsa", "rheMac") -
                         get_pair_val(col, "Mmus", "monDom")
  )
}) %>% bind_rows()
write_tsv(spread_rows, out_spread)

# ── Tree-builder helper ───────────────────────────────────────────────────
build_and_plot <- function(sim_col, out_newick, out_png, out_dist, label_tag) {
  if (!(sim_col %in% names(mat_long))) {
    # Column missing — emit empty outputs to keep snakemake happy
    file.create(out_newick); file.create(out_dist)
    cairo_png_in(out_png, width = 4, height = 2)
    plot.new(); text(0.5, 0.5, paste("missing column:", sim_col))
    dev.off()
    return(invisible(NULL))
  }
  sim_wide <- mat_long %>%
    select(sample_a, sample_b, similarity = !!sym(sim_col)) %>%
    pivot_wider(names_from = sample_b, values_from = similarity) %>%
    column_to_rownames("sample_a")
  sim_wide <- sim_wide[samples, samples]
  diag(sim_wide) <- 1.0

  sim_mat  <- (as.matrix(sim_wide) + t(as.matrix(sim_wide))) / 2
  dist_mat <- 1 - sim_mat
  diag(dist_mat) <- 0

  write_tsv(as_tibble(dist_mat, rownames = "species_code"), out_dist)

  if (length(samples) < 3) {
    tree <- ape::read.tree(text = sprintf("(%s:%.4f, %s:%.4f);",
                                          samples[1], dist_mat[1, 2] / 2,
                                          samples[2], dist_mat[1, 2] / 2))
  } else {
    tree <- ape::nj(as.dist(dist_mat))
    tree <- ape::ladderize(tree)
  }
  ape::write.tree(tree, file = out_newick)

  tip_meta <- tibble(species_code = tree$tip.label) %>%
    left_join(meta, by = "species_code") %>%
    mutate(color = clade_colors[clade])

  cairo_png_in(out_png, width = 9, height = 8)
  op <- par(mar = c(2, 1, 3, 8))
  plot(tree,
       tip.color = tip_meta$color[match(tree$tip.label, tip_meta$species_code)],
       label.offset = 0.005, cex = 1.0, no.margin = FALSE,
       main = sprintf("NJ tree — %s", label_tag))
  if (length(samples) >= 3) ape::add.scale.bar()
  legend("bottomleft", legend = names(clade_colors), col = clade_colors,
         pch = 19, bty = "n", cex = 0.8)
  par(op)
  dev.off()
}

# ── Build the three chosen trees ──────────────────────────────────────────
# reg_eq_col was resolved above, against the actually-loaded data -- not a
# hardcoded literal here, so a missing-column silent-empty-output failure
# (see build_and_plot's guard above) can only occur if neither the new nor
# the legacy column name is present, not merely because the new name is.
build_and_plot(reg_eq_col,                      out_newick_jmz, out_png_jmz, out_dist_jmz, "minimizer Jaccard")
build_and_plot("cosketch_max",                  out_newick_cmz, out_png_cmz, out_dist_cmz, "minimizer cosketch_max")
build_and_plot("cosketch_geom",                 out_newick_gmz, out_png_gmz, out_dist_gmz, "minimizer cosketch_geom")

message("Done. Spreads: ", out_spread)
