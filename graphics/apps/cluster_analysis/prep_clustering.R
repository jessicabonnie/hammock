#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: prep_clustering.R <input_csv_or_tsv> <accession_key.tsv> [out.rds]")
}
input <- args[[1]]
acc_key <- args[[2]]
out <- if (length(args) >= 3) args[[3]] else file.path(dirname(input), paste0(basename(input), ".prep.rds"))

# Try to source utils from common locations
base <- getwd()
cands <- unique(c(
  file.path(base, "scripts", "utils.R"),
  file.path(base, "..", "..", "scripts", "utils.R"),
  file.path(base, "..", "scripts", "utils.R")
))
ok <- FALSE
for (p in cands) { if (file.exists(p)) { try({ source(p); ok <- TRUE }, silent = TRUE); if (ok) break } }
if (!ok) stop("Unable to locate scripts/utils.R")

fmt <- detect_file_format(input)
sim_df <- if (fmt == 'hammock') parse_hammock_format(input) else parse_bedtools_format(input)
meta <- extract_params_from_filename(input)
obj <- list(similarity = as.matrix(sim_df), meta = meta, source_file = basename(input))
saveRDS(obj, out)
cat(sprintf("Wrote cache: %s\n", out))


