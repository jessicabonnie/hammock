#!/usr/bin/env Rscript
# R code to generate and save a heatmap of the data
require(tidyr)
require(dplyr)
require(tibble)
require(gplots)
require(RColorBrewer)
require(readr)
require(stringr)
require(ggplot2)
require(gridExtra)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript encode_heatmap.R <encode.report> <hammock.results> [file_prefix]")
}

# Set default prefix if not provided
file_prefix <- if (length(args) == 3) args[3] else "encode"
print(paste("Using file prefix:", file_prefix))

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
  
  # Debug: print results after join
  # print("Sample after join:")
  # print(head(select(joined_results, file1_id, file2_id, ID_file1, ID_file2)))
  
  return(joined_results)
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

graph_biosample_heatmap <- function(results, matrix_data, mappings) {
  # Get unique biosamples and assign colors
  all_biosamples <- unique(c(results$biosample1, results$biosample2))
  all_biosamples <- all_biosamples[!is.na(all_biosamples)]  # Remove NA values
  biosample_colors <- rainbow(length(all_biosamples))
  names(biosample_colors) <- all_biosamples
  
  # Get colors for rows and columns
  row_colors <- biosample_colors[mappings$file1_biosample[matrix_data$row_names]]
  col_colors <- biosample_colors[mappings$file2_biosample[matrix_data$col_names]]
  
  # Replace NAs with gray
  row_colors[is.na(row_colors)] <- "gray"
  col_colors[is.na(col_colors)] <- "gray"
  
  # Create the heatmap
  heatmap.2(matrix_data$matrix,
           col = brewer.pal(9, "Blues"),
           density.info = "none",
           trace = "none", 
           cexRow = 0.7,
           cexCol = 0.7,
           RowSideColors = row_colors,
           ColSideColors = col_colors,
           margins = c(15, 15),
           main = "Jaccard Similarity by Biosample",
           lhei = c(1, 4),  # Make side colors more prominent
           lwid = c(1, 4))  # Make side colors more prominent
  
  # Add legend
  legend("topright", 
         legend = all_biosamples,
         fill = biosample_colors,
         title = "Biosamples",
         cex = 0.8)
}

graph_organism_heatmap <- function(results, matrix_data, mappings) {
  # Get unique organisms and assign colors
  all_organisms <- unique(c(results$organism1, results$organism2))
  all_organisms <- all_organisms[!is.na(all_organisms)]  # Remove NA values
  organism_colors <- brewer.pal(min(9, max(3, length(all_organisms))), "Set1")
  if(length(all_organisms) > 9) {
    organism_colors <- colorRampPalette(organism_colors)(length(all_organisms))
  }
  names(organism_colors) <- all_organisms
  
  # Get colors for rows and columns
  row_colors <- organism_colors[mappings$file1_organism[matrix_data$row_names]]
  col_colors <- organism_colors[mappings$file2_organism[matrix_data$col_names]]
  
  # Replace NAs with gray
  row_colors[is.na(row_colors)] <- "gray"
  col_colors[is.na(col_colors)] <- "gray"
  
  # Create the heatmap
  heatmap.2(matrix_data$matrix,
           col = brewer.pal(9, "Blues"),
           density.info = "none",
           trace = "none", 
           cexRow = 0.7,
           cexCol = 0.7,
           RowSideColors = row_colors,
           ColSideColors = col_colors,
           margins = c(15, 15),
           main = "Jaccard Similarity by Organism",
           lhei = c(1, 4),  # Make side colors more prominent
           lwid = c(1, 4))  # Make side colors more prominent
  
  # Add legend
  legend("topright", 
         legend = all_organisms,
         fill = organism_colors,
         title = "Organisms",
         cex = 0.8)
}

graph_target_heatmap <- function(results, matrix_data, mappings) {
  # Get unique targets and assign colors
  all_targets <- unique(c(results$target1, results$target2))
  all_targets <- all_targets[!is.na(all_targets)]  # Remove NA values
  target_colors <- brewer.pal(min(9, max(3, length(all_targets))), "Set2")
  if(length(all_targets) > 9) {
    target_colors <- colorRampPalette(target_colors)(length(all_targets))
  }
  names(target_colors) <- all_targets
  
  # Get colors for rows and columns
  row_colors <- target_colors[mappings$file1_target[matrix_data$row_names]]
  col_colors <- target_colors[mappings$file2_target[matrix_data$col_names]]
  
  # Replace NAs with gray
  row_colors[is.na(row_colors)] <- "gray"
  col_colors[is.na(col_colors)] <- "gray"
  
  # Create the heatmap
  heatmap.2(matrix_data$matrix,
           col = brewer.pal(9, "Blues"),
           density.info = "none",
           trace = "none", 
           cexRow = 0.7,
           cexCol = 0.7,
           RowSideColors = row_colors,
           ColSideColors = col_colors,
           margins = c(15, 15),
           main = "Jaccard Similarity by Target of Assay",
           lhei = c(1, 4),  # Make side colors more prominent
           lwid = c(1, 4))  # Make side colors more prominent
  
  # Add legend
  legend("topright", 
         legend = all_targets,
         fill = target_colors,
         title = "Target of Assay",
         cex = 0.8)
}

create_pca_plot <- function(results, matrix_data, mappings) {
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
  pca_df$target <- mappings$file1_target[pca_df$file_name]
  pca_df$biosample <- mappings$file1_biosample[pca_df$file_name]
  pca_df$organism <- mappings$file1_organism[pca_df$file_name]
  
  # Debug: Print samples with PC2 > 2
  outlier_samples <- pca_df[pca_df$PC2 > 2, ]
  if(nrow(outlier_samples) > 0) {
    print("Samples with PC2 > 2:")
    print(outlier_samples)
  }
  
  # Handle cases where there are too many organisms for distinct shapes
  # ggplot2 has about 25 distinct shapes, organisms are typically fewer so this should work well
  organism_counts <- table(pca_df$organism)
  
  # Create a shape mapping for organisms
  shape_values <- c(16, 17, 15, 18, 19, 8, 11, 12, 13, 14, 20, 21, 22, 23, 24, 25, 3, 4, 5, 6)
  
  # Get variance explained for axis labels
  var_explained <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
  
  # Create the plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = biosample, shape = organism)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_shape_manual(values = shape_values[1:length(unique(pca_df$organism))]) +
    facet_wrap(~ target, scales = "free") +
    labs(
      title = "PCA of Jaccard Similarity Matrix by Target of Assay",
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Biosample",
      shape = "Organism"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      strip.text = element_text(size = 9, face = "bold")
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 16)),
      shape = guide_legend(override.aes = list(color = "black"))
    )
  
  return(p)
}

create_pca_plot_by_organism <- function(results, matrix_data, mappings) {
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
  pca_df$target <- mappings$file1_target[pca_df$file_name]
  pca_df$biosample <- mappings$file1_biosample[pca_df$file_name]
  pca_df$organism <- mappings$file1_organism[pca_df$file_name]
  
  # Debug: Print samples with PC2 > 2
  outlier_samples <- pca_df[pca_df$PC2 > 2, ]
  if(nrow(outlier_samples) > 0) {
    print("Samples with PC2 > 2 (from organism-faceted plot):")
    print(outlier_samples)
  }
  
  # Handle cases where there are too many targets for distinct shapes
  # ggplot2 has about 25 distinct shapes, so we'll limit to the most common targets if needed
  target_counts <- table(pca_df$target)
  
  # Create a shape mapping for targets
  shape_values <- c(16, 17, 15, 18, 19, 8, 11, 12, 13, 14, 20, 21, 22, 23, 24, 25, 3, 4, 5, 6)
  
  # Get variance explained for axis labels
  var_explained <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
  
  # Create the plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = biosample, shape = target)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_shape_manual(values = shape_values[1:length(unique(pca_df$target))]) +
    facet_wrap(~ organism, scales = "free") +
    labs(
      title = "PCA of Jaccard Similarity Matrix by Organism",
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Biosample",
      shape = "Target of Assay"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      strip.text = element_text(size = 9, face = "bold")
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 16)),
      shape = guide_legend(override.aes = list(color = "black"))
    )
  
  return(p)
}


# Create a matrix for the heatmap
encode_report <- args[1]
hammock_results <- args[2]
report <- prepare_report(encode_report)
lookup <- create_lookup(report)
hammock_data <- prepare_results(hammock_results, lookup)

# Prepare shared data structures
matrix_data <- create_heatmap_matrix(hammock_data)
mappings <- create_mapping_tables(hammock_data)

# Generate three separate heatmaps and PCA plot
biosample_filename <- paste0(file_prefix, "_biosample.pdf")
organism_filename <- paste0(file_prefix, "_organism.pdf")
target_filename <- paste0(file_prefix, "_target.pdf")
pca_filename <- paste0(file_prefix, "_pca.pdf")

print(paste("Generating biosample heatmap:", biosample_filename))
pdf(biosample_filename, width=12, height=12)
graph_biosample_heatmap(hammock_data, matrix_data, mappings)
dev.off()

print(paste("Generating organism heatmap:", organism_filename))
pdf(organism_filename, width=12, height=12)
graph_organism_heatmap(hammock_data, matrix_data, mappings)
dev.off()

print(paste("Generating target heatmap:", target_filename))
pdf(target_filename, width=12, height=12)
graph_target_heatmap(hammock_data, matrix_data, mappings)
dev.off()

print(paste("Generating combined PCA plots:", pca_filename))
pca_plot1 <- create_pca_plot(hammock_data, matrix_data, mappings)
pca_plot2 <- create_pca_plot_by_organism(hammock_data, matrix_data, mappings)

# Combine both plots into a single PDF
pdf(pca_filename, width = 16, height = 20)
grid.arrange(pca_plot1, pca_plot2, ncol = 1)
dev.off()

print("All heatmaps and PCA plot generated successfully!")

