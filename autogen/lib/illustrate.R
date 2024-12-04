require(tidyr)
require(dplyr)
require(ggplot2)
require(data.table)
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide the directory path as an argument")
}

pattern <- ""
if (length(args) >= 2) {
  pattern <- args[2]
}

# Set working directory to provided path
setwd(args[1])
# setwd("interval_sketch//hammock/cmode/cmode_gen/remap_part/part1k/")

# Get list of all csv files
files_all <- list.files(path=".", pattern="jacc[ABC]\\.csv$")
# When pattern is empty string, grepl() matches everything since empty string is in all strings
files <- files_all[grepl(pattern, files_all)]


dumb_function <- function(filelist){
  # Split into exact and non-exact lists
  exact_files <- filelist[grepl("exact", filelist)]
  sketch_files <- filelist[!grepl("exact", filelist)]
  hll_files <- sketch_files[grepl("hll", sketch_files)]
  balance_files <- sketch_files[grepl("-", sketch_files)]
  mh_files <- sketch_files[grepl("_mh_", sketch_files)]
  hll_balance_files <- sketch_files[grepl("hll", sketch_files) & grepl("-", sketch_files)]
  hll_unbal_files <- sketch_files[grepl("hll", sketch_files) & !grepl("-", sketch_files)]
  mh_balance_files <- sketch_files[grepl("_mh_", sketch_files) & grepl("-", sketch_files)]
  mh_unbal_files <- sketch_files[grepl("_mh_", sketch_files) & !grepl("-", sketch_files)]
  exact_balance_files <- filelist[grepl("exact", filelist) & grepl("-", filelist)]
  exact_unbal_files <- filelist[grepl("exact", filelist) & !grepl("-", filelist)]
 
    combine_versions <- function(list_files, balval=FALSE){
      data.table::rbindlist(lapply(list_files, function(x) {
    dt <- data.table::fread(x)
    if (! "balance" %in% names(dt)) dt[, balance := balval]
    
    if ("intersect_inexact" %in% names(dt)) data.table::setnames(dt, "intersect_inexact", "intersect")
    dt}),
     use.names = TRUE)
  }

  # Read and combine all files using data.table
  dt_balance <- combine_versions(balance_files, balval = TRUE) 
  dt_hll_bal <- combine_versions(hll_balance_files, balval = TRUE)
  dt_hll_unbal <- combine_versions(hll_unbal_files, balval = FALSE)
  dt_hll <- rbind(dt_hll_bal, dt_hll_unbal)

  dt_mh_bal <- combine_versions(mh_balance_files, balval = TRUE)
  dt_mh_unbal <- combine_versions(mh_unbal_files, balval = FALSE)
  dt_mh <- rbind(dt_mh_bal, dt_mh_unbal)

  dt_exact_bal <- combine_versions(exact_balance_files, balval = TRUE)
  dt_exact_unbal <- combine_versions(exact_unbal_files, balval = FALSE)
dt_exact <- rbind(dt_exact_bal, dt_exact_unbal)

  dt <- rbind(dt_hll, dt_exact, dt_mh)
  dt_balance <- rbind(dt_hll_bal, dt_exact_bal, dt_mh_bal)
  dt_unbal <- rbind(dt_hll_unbal, dt_exact_unbal, dt_mh_unbal)

  dt_unbal_wider <- dt_unbal %>%
    pivot_wider(
      id_cols = c(bed1, bed2, mode, subsample, balance),
      names_from = sketch_type,
      values_from = jaccardfunc,
      names_sep = "_"
  ) 

  dt_bal_wider <- dt_balance %>%
    pivot_wider(
      id_cols = c(bed1, bed2, mode, subsample, balance),
      names_from = sketch_type,
      values_from = jaccardfunc,
      names_sep = "_"
  ) 

  dt_wider <- rbind(dt_bal_wider, dt_unbal_wider)

  add_gen_vals <- function(dtin){
    dtin %>% 
    mutate(genmode = case_when(
        grepl("modeA", bed2) ~ "A",
        grepl("modeB", bed2) ~ "B",
        grepl("no_overlap.bed", bed2) ~ "AB",
        TRUE ~ "NA"  # default case
    )) %>%
    mutate(genfrac = case_when(
        # Look for pattern like "_0.1.bed" or "_0.33.bed" etc
        grepl("_0\\.[0-9]+\\.bed", bed2) ~ {
            extracted <- gsub(".*_(0\\.[0-9]+)\\.bed.*", "\\1", bed2)
            as.numeric(extracted)},
        TRUE ~ 0
    )) %>% # Replace any NA with 0
    mutate(genfrac = ifelse(is.na(genfrac), 0, genfrac))  %>% 
    mutate(genpad = if_else(
        condition = grepl("nopad", bed2),
        true = "nopad",
        false = "pad"
    ))
  }
  dt_genx <- add_gen_vals(dt_wider)
  dt_exact_genx <- add_gen_vals(dt_exact)
  return(list(dt_genx, dt_exact_genx))
}
dts <- dumb_function(files)
dt_gen <- dts[[1]]
dt_exact_gen <- dts[[2]]

dt_exact_other <- dt_exact_gen %>%
  pivot_wider(
    id_cols = c(bed1, bed2, subsample, genmode, genfrac,
                genpad, actfrac, balance),
    names_from = mode,
    values_from = c(union, intersect, jaccardfunc),
    names_sep = "_"
  ) %>%
  mutate(A_multiplier = ifelse(balance,1 - abs(subsample), 1)) %>%
  mutate(fakeCUnion = union_B * abs(subsample) + union_A * A_multiplier,
    fakeCIntersect = intersect_B * abs(subsample) + intersect_A * A_multiplier,
    fakeCJaccard = fakeCIntersect / fakeCUnion
  )

ggplot(dt_gen) + geom_point(aes(x=exact, y=hyperloglog))

ggplot(filter(dt_gen, mode != "C")) + geom_point(aes(x=exact, y=hyperloglog, color=balance))
