#!/usr/bin/env Rscript
# R code to generate PCA plots with separate calculations for each organism
require(tidyr)
require(dplyr)
require(tibble)
require(readr)
require(stringr)
require(ggplot2)
require(gridExtra)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript encode_PCA.R <encode.report> <hammock.results> [file_prefix]")
}

# Set default prefix if not provided
file_prefix <- if (length(args) == 3) args[3] else "encode"

# Functions

prepare_report <- function(filepath){
  # First read the data to see what columns are available
  report_raw <- read.table(
    filepath,
    skip = 1, sep = "\t", header=TRUE, stringsAsFactors = FALSE)
  

  
  # Check for common variations of column names
  # R converts spaces to periods in column names, so "Target of assay" becomes "Target.of.assay"
  target_col <- NULL
  if("Target.of.assay" %in% colnames(report_raw)) {
    target_col <- "Target.of.assay"
  } else if("Target" %in% colnames(report_raw)) {
    target_col <- "Target"
  } else if("Assay" %in% colnames(report_raw)) {
    target_col <- "Assay"
  } else if("Experiment.target" %in% colnames(report_raw)) {
    target_col <- "Experiment.target"
  } else {
    # Find any column that contains "target" and "assay" or similar patterns
    target_candidates <- grep("target.*assay|assay.*target|target.of|experiment.*target", colnames(report_raw), ignore.case = TRUE, value = TRUE)
    
    if(length(target_candidates) == 0) {
      # Broaden search to any column with "target" or "assay"
      target_candidates <- grep("target|assay", colnames(report_raw), ignore.case = TRUE, value = TRUE)
    }
    
    if(length(target_candidates) > 0) {
      # Try to pick the most likely one based on common ENCODE patterns
      # Priority order: most specific to least specific
      priority_patterns <- c(
        "^Target\\.of\\.assay$",      # Exact match for Target.of.assay
        "Target\\.of\\.",             # Target.of.* (like Target.of.assay)
        "Experiment\\.target",        # Experiment.target
        "Assay\\.title",             # Assay.title
        "Target\\.label",            # Target.label
        "^Target$"                   # Just "Target"
      )
      
      for(pattern in priority_patterns) {
        matches <- target_candidates[grep(pattern, target_candidates, ignore.case = TRUE)]
        if(length(matches) > 0) {
          target_col <- matches[1]
          break
        }
      }
      
      # If no priority patterns match, use the first one
      if(is.null(target_col)) {
        target_col <- target_candidates[1]
      }
    } else {
      stop("Could not find a target/assay column. Available columns: ", paste(colnames(report_raw), collapse = ", "))
    }
  }
  
  # Select columns with the identified target column
  report <- report_raw %>%
    select(Accession, Biosample.term.name, Organism, all_of(target_col), Files) %>%
    rename(Target.of.assay = all_of(target_col)) %>%
    mutate(ID=paste(Accession, Biosample.term.name, Organism, Target.of.assay, sep=":"))
  
  return(report)
}

create_lookup <- function(report){
  # First, convert the Files column to character in case it's a factor
  report$Files <- as.character(report$Files)
  
  # Create a lookup table by extracting ENCF IDs from the Files column
  lookup_table <- report %>%
    # Split the Files string into individual files
    mutate(file_ids = strsplit(Files, ",")) %>%
    # Unnest to get one row per file ID
    unnest(file_ids) %>%
    # Extract just the ENCF ID from the path using regex
    mutate(file_id = str_extract(file_ids, "ENCF[F0-9A-Z]+")) %>%
    # Remove rows where no ENCF ID was found
    filter(!is.na(file_id)) %>%
    # Keep only the relevant columns
    select(file_id, ID, Biosample.term.name, Organism, Target.of.assay)
  

  
  return(lookup_table)
}

create_heatmap_matrix <- function(results) {
  # Create the jaccard matrix for the heatmap
  jaccard_subtable <- results  %>%
    select(file1, file2, jaccard) %>%
    arrange(file1, file2) %>%
    pivot_wider(names_from = file2, 
                values_from = jaccard) %>%
    as.data.frame() 
  
  # Extract row names and prepare matrix
  row_names <- jaccard_subtable$file1
  row.names(jaccard_subtable) <- row_names
  jaccard_matrix <- as.matrix(jaccard_subtable[, -1])
  
  return(list(matrix = jaccard_matrix, row_names = row_names, col_names = colnames(jaccard_matrix)))
}

create_mapping_tables <- function(results) {
  # Create mappings between file IDs and their attributes
  file1_biosample_map <- results %>%
    select(file1, biosample1) %>%
    distinct() %>%
    deframe()
    
  file2_biosample_map <- results %>%
    select(file2, biosample2) %>%
    distinct() %>%
    deframe()
    
  file1_organism_map <- results %>%
    select(file1, organism1) %>%
    distinct() %>%
    deframe()
    
  file2_organism_map <- results %>%
    select(file2, organism2) %>%
    distinct() %>%
    deframe()

  file1_target_map <- results %>%
    select(file1, target1) %>%
    distinct() %>%
    deframe()
    
  file2_target_map <- results %>%
    select(file2, target2) %>%
    distinct() %>%
    deframe()
    
  return(list(
    file1_biosample = file1_biosample_map,
    file2_biosample = file2_biosample_map,
    file1_organism = file1_organism_map,
    file2_organism = file2_organism_map,
    file1_target = file1_target_map,
    file2_target = file2_target_map
  ))
}

prepare_results <- function(filepath, lookup, simcol="jaccard_similarity_with_ends"){
  # Read and prepare the results
  results <- readr::read_delim(filepath, show_col_types = FALSE) %>% 
    rename_with(~"jaccard", all_of(simcol)) %>%
    select(file1, file2, jaccard) %>%
    # Extract just the ENCF ID from the filename (remove all extensions and path)
    # Handle various file formats: .fa, .bed, .bed.gz, .narrowPeak, etc.
    mutate(
      file1_id = str_extract(basename(file1), "ENCF[F0-9A-Z]+"),
      file2_id = str_extract(basename(file2), "ENCF[F0-9A-Z]+")
    ) %>%
    # Debug: show what we extracted
    mutate(
      file1_orig = file1,
      file2_orig = file2
    )
  

  
  # Perform the joins with the clean IDs
  joined_results <- results %>%
    left_join(lookup, by = c("file1_id" = "file_id")) %>%
    rename(
      biosample1 = Biosample.term.name,
      organism1 = Organism,
      target1 = Target.of.assay,
      ID_file1 = ID
    ) %>%
    left_join(lookup, by = c("file2_id" = "file_id")) %>%
    rename(
      biosample2 = Biosample.term.name,
      organism2 = Organism,
      target2 = Target.of.assay,
      ID_file2 = ID
    )
  
  return(joined_results)
}

create_organism_specific_pca <- function(results, organism_name) {
  # Filter results to only include pairs where both files are from the same organism
  organism_results <- results %>%
    filter(organism1 == organism_name & organism2 == organism_name)
  
  cat("  Found", nrow(organism_results), "pairs for", organism_name, "\n")
  
  if(nrow(organism_results) == 0) {
    cat("  No data found for", organism_name, "\n")
    return(NULL)
  }
  
  # Create the jaccard matrix for this organism
  jaccard_subtable <- organism_results %>%
    select(file1, file2, jaccard) %>%
    arrange(file1, file2) %>%
    pivot_wider(names_from = file2, 
                values_from = jaccard) %>%
    as.data.frame() 
  
  # Extract row names and prepare matrix
  row_names <- jaccard_subtable$file1
  row.names(jaccard_subtable) <- row_names
  jaccard_matrix <- as.matrix(jaccard_subtable[, -1])
  
  # Replace NA values with 0 for PCA
  jaccard_matrix[is.na(jaccard_matrix)] <- 0
  
  # Perform PCA on the organism-specific similarity matrix
  pca_result <- prcomp(jaccard_matrix, center = TRUE, scale. = TRUE)
  
  # Create a data frame with PCA results and metadata
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    file_name = rownames(pca_result$x)
  )
  
  # Add metadata
  biosample_map <- organism_results %>%
    select(file1, biosample1) %>%
    distinct() %>%
    deframe()
  
  # Handle cases where mapping might be empty
  if(length(biosample_map) > 0) {
    pca_df$biosample <- biosample_map[pca_df$file_name]
  } else {
    pca_df$biosample <- NA
  }
  
  pca_df$organism <- organism_name
  
  # Get variance explained for axis labels
  var_explained <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
  
  return(list(
    pca_df = pca_df,
    var_explained = var_explained,
    organism = organism_name
  ))
}

create_organism_pca_plot <- function(pca_data) {
  if(is.null(pca_data)) {
    return(NULL)
  }
  
  pca_df <- pca_data$pca_df
  var_explained <- pca_data$var_explained
  organism <- pca_data$organism
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = biosample)) +
    geom_point(size = 5, alpha = 0.8) +
    labs(
      title = paste0("PCA of Jaccard Similarity Matrix - ", organism),
      subtitle = paste0("PC1 (", var_explained[1], "% variance) vs PC2 (", var_explained[2], "% variance)"),
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Biosample"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  return(p)
}

create_global_pca_plot <- function(results, matrix_data, mappings) {
  # Prepare data for PCA
  # Use the jaccard similarity matrix, but we need to handle missing values
  pca_matrix <- matrix_data$matrix
  
  # Replace NA values with 0 for PCA
  pca_matrix[is.na(pca_matrix)] <- 0
  
  # Perform PCA on the similarity matrix
  pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
  
  # Create a data frame with PCA results and metadata
  # Note: Each row in the similarity matrix represents one sample's similarity profile
  # to all other samples, so PC coordinates represent each sample's position in similarity space
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    file_name = rownames(pca_result$x)
  )
  
  # Add metadata from mappings (using file1/row sample metadata since PCA rows = samples)
  pca_df$biosample <- mappings$file1_biosample[pca_df$file_name]
  pca_df$organism <- mappings$file1_organism[pca_df$file_name]
  
  # Create a shape mapping for organisms
  shape_values <- c(16, 17, 15, 18, 19, 8, 11, 12, 13, 14, 20, 21, 22, 23, 24, 25, 3, 4, 5, 6)
  
  # Get variance explained for axis labels
  var_explained <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
  
  # Create the plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = biosample, shape = organism)) +
    geom_point(size = 5, alpha = 0.8) +
    scale_shape_manual(values = shape_values[1:length(unique(pca_df$organism))]) +
    labs(
      title = "PCA of Jaccard Similarity Matrix (All Organisms)",
      subtitle = paste0("PC1 (", var_explained[1], "% variance) vs PC2 (", var_explained[2], "% variance)"),
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Biosample",
      shape = "Organism"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 16)),
      shape = guide_legend(override.aes = list(color = "black"))
    )
  
  return(p)
}

# Main execution
encode_report <- args[1]
hammock_results <- args[2]
report <- prepare_report(encode_report)
lookup <- create_lookup(report)
hammock_data <- prepare_results(hammock_results, lookup)

# Prepare shared data structures for global PCA
matrix_data <- create_heatmap_matrix(hammock_data)
mappings <- create_mapping_tables(hammock_data)

# Generate global PCA plot (all organisms together)
cat("Generating global PCA plot...\n")
global_pca_filename <- paste0(file_prefix, "_pca.pdf")

tryCatch({
  global_pca_plot <- create_global_pca_plot(hammock_data, matrix_data, mappings)
  pdf(global_pca_filename, width = 16, height = 12)
  print(global_pca_plot)
  dev.off()
  cat("Global PCA plot saved to:", global_pca_filename, "\n")
}, error = function(e) {
  print(paste("Error generating global PCA plot:", e$message))
  if (dev.cur() != 1) {
    dev.off()
  }
})

# Get unique organisms
organisms <- unique(c(hammock_data$organism1, hammock_data$organism2))
organisms <- organisms[!is.na(organisms)]

# Generate PCA for each organism
pca_plots <- list()
pca_data_list <- list()

for(org in organisms) {
  cat("Processing PCA for organism:", org, "\n")
  pca_data <- create_organism_specific_pca(hammock_data, org)
  if(!is.null(pca_data)) {
    pca_data_list[[org]] <- pca_data
    pca_plots[[org]] <- create_organism_pca_plot(pca_data)
  }
}

# Create PDF with one plot per page
pca_filename <- paste0(file_prefix, "_organism_specific_pca.pdf")

tryCatch({
  pdf(pca_filename, width = 16, height = 12)
  
  for(org in names(pca_plots)) {
    if(!is.null(pca_plots[[org]])) {
      print(pca_plots[[org]])
      cat("Generated plot for:", org, "\n")
    }
  }
  
  dev.off()
  cat("Organism-specific PCA plots saved to:", pca_filename, "\n")
  
}, error = function(e) {
  print(paste("Error generating organism-specific PCA plots:", e$message))
  # Close the PDF device if it's still open
  if (dev.cur() != 1) {
    dev.off()
  }
})

# Print summary statistics
cat("\nPCA Summary:\n")
for(org in names(pca_data_list)) {
  if(!is.null(pca_data_list[[org]])) {
    n_samples <- nrow(pca_data_list[[org]]$pca_df)
    var_explained <- pca_data_list[[org]]$var_explained
    cat(sprintf("%s: %d samples, PC1=%.1f%%, PC2=%.1f%%\n", 
                org, n_samples, var_explained[1], var_explained[2]))
  }
} 