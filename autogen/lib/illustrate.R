require(tidyr)
require(dplyr)
require(ggplot2)
require(data.table)
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide the directory path as an argument")
}

# Set working directory to provided path
setwd(args[1])
# setwd("interval_sketch//hammock/cmode/cmode_gen/remap_part/part1k/")

# Get list of all csv files
files <- list.files(path=".", pattern="jacc[ABC]\\.csv$")

# Split into exact and non-exact lists
exact_files <- files[grepl("exact", files)]
sketch_files <- files[!grepl("exact", files)]
hll_files <- files[grepl("hll", files)]
balance_files <- files[grepl("-", files)]
mh_files <- files[grepl("_mh_", files)]
hll_balance_files <- files[grepl("hll", files) & grepl("-", files)]
hll_unbal_files <- files[grepl("hll", files) & !grepl("-", files)]
mh_balance_files <- files[grepl("_mh_", files) & grepl("-", files)]
mh_unbal_files <- files[grepl("_mh_", files) & !grepl("-", files)]


exact_balance_files <- files[grepl("exact", files) & grepl("-", files)]
exact_unbal_files <- files[grepl("exact", files) & !grepl("-", files)]

combine_versions <- function(list_files, balval=FALSE){
  data.table::rbindlist(lapply(list_files, function(x) {
  dt <- data.table::fread(x)
  if (! "balance" %in% names(dt)) dt[, balance := balval]
  
  if ("intersect_inexact" %in% names(dt)) data.table::setnames(dt, "intersect_inexact", "intersect")
  dt
}), use.names = TRUE)
}

# Read and combine all files using data.table
# dt <- data.table::rbindlist(lapply(files, data.table::fread))
dt_balance <- combine_versions(balance_files, balval = TRUE) 
#data.table::rbindlist(lapply(balance_files, data.table::fread))
# dt_hll <- data.table::rbindlist(lapply(hll_files, data.table::fread))
dt_hll_bal <- combine_versions(hll_balance_files, balval = TRUE)
dt_hll_unbal <- combine_versions(hll_unbal_files, balval = FALSE)
dt_hll <- rbind(dt_hll_bal, dt_hll_unbal)
dt_mh_bal <- combine_versions(mh_balance_files, balval = TRUE)
dt_mh_unbal <- combine_versions(mh_unbal_files, balval = FALSE)
dt_mh <- rbind(dt_mh_bal, dt_mh_unbal)
# dt_exact <- data.table::rbindlist(lapply(exact_files, data.table::fread))
dt_exact_bal <- combine_versions(exact_balance_files, balval = TRUE)
dt_exact_unbal <- combine_versions(exact_unbal_files, balval = FALSE)
dt_exact <- rbind(dt_exact_bal, dt_exact_unbal)

dt <- rbind(dt_hll, dt_exact, dt_mh)
dt_balance <- rbind(dt_hll_bal, dt_exact_bal, dt_mh_bal)
dt_unbal <- rbind(dt_hll_unbal, dt_exact_unbal, dt_mh_unbal)
# Widen the data table keeping bed1, bed2, mode, and subsample as static columns
# dt_exact <- dt %>% filter(sketch_type == "exact") %>%
#   pivot_wider(
#     id_cols = c(bed1, bed2, mode, subsample),
#     names_from = sketch_type,
#     values_from = c(jaccardfunc, union, intersect_inexact),
#     names_sep = "_"
#   )

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
        # Look for pattern like "_2." or "_3." etc
        grepl("_\\d\\.", bed2) ~ as.numeric(gsub(".*_(\\d)\\.*.*", "\\1", bed2)),
        TRUE ~ 0
    )) %>%
    mutate(genpad = if_else(
        condition = grepl("nopad", bed2),
        true = "nopad",
        false = "pad"
    )) %>%
    mutate(actfrac = case_when( 
      genfrac == 0 ~ 0,
      genpad == "nopad" ~ 1/genfrac,
      genpad == "pad" ~ 1/((2*genfrac)-1)))
}

dt_gen <- add_gen_vals(dt_wider)
dt_exact_gen <- add_gen_vals(dt_exact)

dt_exact_other <- dt_exact_gen %>%
   pivot_wider(
    id_cols = c(bed1, bed2, subsample, genmode, genfrac, genpad, actfrac, balance),
    names_from = mode,
    values_from =c(union, intersect, jaccardfunc),
    names_sep = "_"
  ) %>%
  mutate(A_multiplier = ifelse(balance,1-abs(subsample), 1)) %>%
  mutate(fakeCUnion = union_B * abs(subsample) + union_A * A_multiplier,
  fakeCIntersect = intersect_B * abs(subsample) + intersect_A * A_multiplier,
  fakeCJaccard = fakeCIntersect / fakeCUnion
  ) 

   ggplot(dt_gen) + geom_point(aes(x=exact, y=hyperloglog))
