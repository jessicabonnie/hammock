#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: crosscheck_r_clustering.R INPUT.csv METRIC OUTPUT.csv")
}

input <- read.csv(args[[1]], check.names = FALSE)
metric <- args[[2]]
output <- args[[3]]
required <- c("file1", "file2", metric)
if (!all(required %in% names(input))) {
  stop(sprintf("input lacks named columns: %s", paste(setdiff(required, names(input)), collapse = ", ")))
}

strip_known <- function(value) {
  name <- basename(value)
  sub("\\.(fasta|fna|fa|bed)(\\.gz)?$", "", name, ignore.case = TRUE)
}
left <- vapply(input$file1, strip_known, character(1))
right <- vapply(input$file2, strip_known, character(1))
samples <- unique(left)
if (!setequal(samples, unique(right)) || length(samples) != 20) {
  stop("expected identical 20-sample axes")
}
matrix <- matrix(NA_real_, nrow = length(samples), ncol = length(samples),
                 dimnames = list(samples, samples))
for (index in seq_len(nrow(input))) {
  matrix[left[[index]], right[[index]]] <- as.numeric(input[[metric]][[index]])
}
if (any(!is.finite(matrix)) || !isTRUE(all.equal(matrix, t(matrix), tolerance = 1e-9)) ||
    !isTRUE(all.equal(unname(diag(matrix)), rep(1, length(samples)), tolerance = 1e-12))) {
  stop("matrix invariants failed")
}

tree <- hclust(as.dist(1 - matrix), method = "average")
partition <- cutree(tree, k = 10)
rows <- expand.grid(file1 = samples, file2 = samples, stringsAsFactors = FALSE)
rows$cocluster <- as.integer(partition[rows$file1] == partition[rows$file2])
rows$cluster <- unname(partition[rows$file1])
write.csv(rows, output, row.names = FALSE, quote = FALSE)
