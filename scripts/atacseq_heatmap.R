# R code to generate and save a heatmap of the data
require (tidyverse)
require(gplots)
require(RColorBrewer)
require(readr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Please provide two file paths: encode.report and hammock.results")
}

encode_report <- args[1]
hammock_results <- args[2]

# Functions


graph_me <- function(results) {
  # First, create mappings between file IDs and their biosamples and organisms
  # This preserves the relationship before pivoting
  file1_biosample_map <- results %>%
    select(file1, biosample1) %>%
    distinct() %>%
    deframe() # Converts to a named vector
    
  file2_biosample_map <- results %>%
    select(file2, biosample2) %>%
    distinct() %>%
    deframe()
    
  # Add mappings for organisms
  file1_organism_map <- results %>%
    select(file1, organism1) %>%
    distinct() %>%
    deframe()
    
  file2_organism_map <- results %>%
    select(file2, organism2) %>%
    distinct() %>%
    deframe()
    
  # Get unique biosamples and assign colors
  all_biosamples <- unique(c(results$biosample1, results$biosample2))
  biosample_colors <- rainbow(length(all_biosamples))
  names(biosample_colors) <- all_biosamples
  
  # Get unique organisms and assign colors
  all_organisms <- unique(c(results$organism1, results$organism2))
  organism_colors <- brewer.pal(min(9, max(3, length(all_organisms))), "Set1")
  if(length(all_organisms) > 9) {
    organism_colors <- colorRampPalette(organism_colors)(length(all_organisms))
  }
  names(organism_colors) <- all_organisms
  
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
  
  # Column names are in colnames(jaccard_matrix)
  col_names <- colnames(jaccard_matrix)
  
  # Now use the mapping to get colors for rows and columns
  row_biosample_colors <- biosample_colors[file1_biosample_map[row_names]]
  col_biosample_colors <- biosample_colors[file2_biosample_map[col_names]]
  
  row_organism_colors <- organism_colors[file1_organism_map[row_names]]
  col_organism_colors <- organism_colors[file2_organism_map[col_names]]
  
  # Debug if needed
  print("Row biosamples match check:")
  print(head(data.frame(
    file = row_names,
    biosample = file1_biosample_map[row_names],
    organism = file1_organism_map[row_names]
  )))
  
  # Since manual drawing doesn't work, let's try a simpler approach:
  # Create two separate heatmaps side by side or use a different strategy
  
  # Strategy: Create organism-specific symbols or use different approach
  # For now, let's make a pattern where organism info is encoded differently
  
  # Create patterned/modified colors that encode organism information
  # We'll modify the biosample colors based on organism
  
  # Create modified colors: darker for one organism, lighter for another
  create_organism_modified_color <- function(base_color, organism, all_organisms) {
    if(length(all_organisms) == 2) {
      if(organism == all_organisms[1]) {
        # Make darker for first organism
        return(adjustcolor(base_color, alpha.f = 1.0, red.f = 0.7, green.f = 0.7, blue.f = 0.7))
      } else {
        # Keep original for second organism  
        return(base_color)
      }
    } else {
      # For more than 2 organisms, use alpha variations
      org_index <- which(all_organisms == organism)
      alpha_val <- 0.4 + 0.6 * (org_index / length(all_organisms))
      return(adjustcolor(base_color, alpha.f = alpha_val))
    }
  }
  
  # Apply organism modifications to biosample colors
  row_modified_colors <- mapply(create_organism_modified_color, 
                               row_biosample_colors, 
                               file1_organism_map[row_names],
                               MoreArgs = list(all_organisms = all_organisms))
  
  col_modified_colors <- mapply(create_organism_modified_color, 
                               col_biosample_colors, 
                               file2_organism_map[col_names],
                               MoreArgs = list(all_organisms = all_organisms))
  
  print("Creating heatmap with organism-modified biosample colors...")
  print(paste("Organisms:", paste(all_organisms, collapse = ", ")))
  
  # Create the heatmap with modified colors that encode both biosample and organism
  heatmap.2(jaccard_matrix,
           col = brewer.pal(9, "Blues"),
           density.info = "none",
           trace = "none", 
           cexRow = 0.7,
           cexCol = 0.7,
           RowSideColors = row_modified_colors,
           ColSideColors = col_modified_colors,
           margins = c(12, 12))
  
  # Add comprehensive legends explaining the encoding
  legend("topright", 
         legend = all_biosamples,
         fill = biosample_colors,
         title = "Biosamples (Base Colors)",
         cex = 0.6)
         
  # Create legend for organism encoding
  if(length(all_organisms) == 2) {
    org_legend_colors <- c(
      adjustcolor(biosample_colors[1], alpha.f = 1.0, red.f = 0.7, green.f = 0.7, blue.f = 0.7),
      biosample_colors[1]
    )
    org_legend_labels <- paste(all_organisms, c("(darker)", "(original)"))
  } else {
    org_legend_colors <- sapply(1:length(all_organisms), function(i) {
      alpha_val <- 0.4 + 0.6 * (i / length(all_organisms))
      adjustcolor(biosample_colors[1], alpha.f = alpha_val)
    })
    org_legend_labels <- paste(all_organisms, "(varying intensity)")
  }
  
  legend("bottomright", 
         legend = org_legend_labels,
         fill = org_legend_colors,
         title = "Organism Encoding",
         cex = 0.6)
         
  # Add explanatory text
  mtext("Color intensity indicates organism", side = 1, line = -2, cex = 0.8)
}





prepare_report <- function(filepath){
  report <- read.table(
  filepath,
  skip = 1, sep = "\t", header=TRUE, stringsAsFactors = FALSE) %>%
  select(Accession, Biosample.term.name, Organism, Files) %>%
  mutate(ID=paste(Accession, Biosample.term.name, Organism, sep=":"))
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
    # Extract just the ENCF ID from the path
    mutate(file_id = gsub("^.*/([^/]+)/$", "\\1", file_ids)) %>%
    # Keep only the relevant columns
    select(file_id, ID, Biosample.term.name, Organism)
  
  # Debug: print a few rows to verify extraction worked
  print(head(lookup_table))
  
  return(lookup_table)
}

prepare_results <- function(filepath, lookup, simcol="jaccard_similarity"){
  # Read and prepare the results
  results <- readr::read_delim(filepath, show_col_types = FALSE) %>% 
    rename_with(~"jaccard", all_of(simcol)) %>%
    select(file1, file2, jaccard) %>%
    # Remove file extensions
    mutate(
      file1_clean = sub("\\.(bed|bed\\.gz)$", "", file1),
      file2_clean = sub("\\.(bed|bed\\.gz)$", "", file2)
    ) %>%
    # Extract just the base filename (ENCFXXXX portion)
    mutate(
      file1_id = basename(file1_clean),
      file2_id = basename(file2_clean)
    )
  
  # Debug: print samples from results and lookup before joining
  print("Sample of file IDs in results before join:")
  print(head(select(results, file1_id, file2_id)))
  
  print("Sample of file IDs in lookup table:")
  print(head(select(lookup, file_id)))
  
  # Perform the joins with the clean IDs
  joined_results <- results %>%
    left_join(lookup, by = c("file1_id" = "file_id")) %>%
    rename(
      biosample1 = Biosample.term.name,
      organism1 = Organism,
      ID_file1 = ID
    ) %>%
    left_join(lookup, by = c("file2_id" = "file_id")) %>%
    rename(
      biosample2 = Biosample.term.name,
      organism2 = Organism,
      ID_file2 = ID
    )
  
  # Debug: print results after join
  print("Sample after join:")
  print(head(select(joined_results, file1_id, file2_id, ID_file1, ID_file2)))
  
  return(joined_results)
}

# Create a matrix for the heatmap
# filepaths
report_path <- paste0(atac.seq,"report_2025_0205/experiment_report_2025_2_5_20h_24m.tsv")
results_path <-paste0(atac.seq,"report_2025_0205/data/mouseonly_hll_p20_jaccC.csv")
# Generate heatmap using base R
png("heatmap.png", width=800, height=600)
graph_me(hammock_data)
dev.off()

