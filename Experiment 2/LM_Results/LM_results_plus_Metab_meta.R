# ============================================================
# Join Downsampling Results with Metabolite Metadata
#
# Reads the per-metabolite, per-iteration linear model results
# produced by DownSamplingByPoolSizeLM_capped_resampling.R
# (LM_population_removal_all_iterations.csv) and joins them against
# a metabolite metadata key (name, metabolism category, HMDB ID) to
# produce a single combined file used downstream for functional
# module analysis (Function.R).
#
# A handful of columns not needed for that downstream analysis
# (Diet_abs_coef, Diet_sig, Removed_STD, Removed_HSD, Diet_p, Column,
# Ion.Mode, Dataset) are dropped to keep the combined file smaller.
#
# Output: LM_downsample_results_META.csv and .rds, both written to
# this folder.
# ============================================================

# set working directory
# This script expects to be run with the working directory set to
# Experiment 2/LM_Results/ within the repository. It reads
# LM_population_removal_all_iterations.csv from this folder and
# ../Key/Metab_meta_key.csv, and writes LM_downsample_results_META.csv
# and .rds here.

# =========================
# Load packages
# =========================
library(dplyr)

# =========================
# Clear environment
# =========================
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)


# read files
results_df <- read.csv("LM_population_removal_all_iterations.csv", header = TRUE, stringsAsFactors = FALSE)
meta_key <- read.csv("../Key/Metab_meta_key.csv", header = TRUE, stringsAsFactors = FALSE)

# Add metabolite meta data to results
results_with_meta <- results_df %>%
  left_join(meta_key, by = "Metabolite")

# remove excess data

results_small <- results_with_meta %>%
  select(-Diet_abs_coef, 
         -Diet_sig, 
         -Removed_STD, 
         -Removed_HSD, 
         -Diet_p,
         -Column,
         -Ion.Mode,
         -Dataset)

library(data.table) # use this to write big files faster

fwrite(
  results_small,
  "LM_downsample_results_META.csv"
)

# use this filetype for analysis later
saveRDS(results_small, "LM_downsample_results_META.rds")

