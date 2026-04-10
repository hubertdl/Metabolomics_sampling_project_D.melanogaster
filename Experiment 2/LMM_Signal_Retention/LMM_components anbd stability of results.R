# ============================================================
# Modeling Determinants of Metabolite Signal Retention
#
# This script evaluates the factors that determine whether 
# individual metabolites remain detectable across replicate 
# downsampling and different pool-size conditions.
#
# Workflow:
# 1. Import combined linear model results from the downsampling 
#    analysis across metabolite panels.
# 2. Format pool size, replicate number, and dataset variables 
#    for downstream modeling.
# 3. Define a reference condition (PoolSize = 100, full replication) 
#    and calculate the reference absolute diet effect size for each 
#    metabolite.
# 4. Quantify metabolite signal retention as the frequency with 
#    which each metabolite remains significant (FDR < 0.05) across 
#    downsampling iterations.
# 5. Quantify variability in metabolite effect estimates as the 
#    standard deviation of diet coefficients across iterations.
# 6. Combine reference effect size, detection frequency, and 
#    coefficient variability into a single modeling dataset.
# 7. Fit a linear mixed-effects model testing how reference effect 
#    size, coefficient variability, pool size, and replicate number 
#    influence metabolite signal retention, while accounting for 
#    differences among datasets and metabolites.
#
# This analysis identifies the major drivers of metabolite signal 
# retention under reduced sampling and quantifies how biological 
# effect size and estimate stability shape the robustness of 
# metabolomic inference.
# ============================================================

# set working directory
setwd("")

rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

library(dplyr)

# =========================
# Load data
# =========================
results_df <- readRDS("./LM_downsample_results_META.rds") # large output file from downsampling analysis

# =========================
# Basic preprocessing
# =========================
results_df <- results_df %>%
  mutate(
    PoolSize = as.character(PoolSize),
    N_remaining = as.numeric(as.character(N_remaining))
  )

results_df <- results_df %>%
  rename(Dataset = Panel)

# define reference

full_n <- max(results_df$N_remaining, na.rm = TRUE)

reference_effects <- results_df %>%
  filter(
    PoolSize == "100",
    N_remaining == full_n
  ) %>%
  group_by(Dataset, Metabolite_meta) %>%
  summarise(
    ref_abs_effect = mean(abs(Diet_coef), na.rm = TRUE),
    .groups = "drop"
  )

# build retention (response variable)

stability_df <- results_df %>%
  mutate(sig = Diet_FDR < 0.05) %>%
  group_by(Dataset, PoolSize, N_remaining, Metabolite_meta) %>%
  summarise(
    sig_freq = mean(sig, na.rm = TRUE),
    .groups = "drop"
  )

# quantify varibility

variance_df <- results_df %>%
  group_by(Dataset, PoolSize, N_remaining, Metabolite_meta) %>%
  summarise(
    coef_sd = sd(Diet_coef, na.rm = TRUE),
    .groups = "drop"
  )

# combine into one modeling dataset

model_df <- stability_df %>%
  left_join(reference_effects, by = c("Dataset", "Metabolite_meta")) %>%
  left_join(variance_df, by = c("Dataset", "PoolSize", "N_remaining", "Metabolite_meta"))

summary(model_df)

# turn coef_sd nas into 0
model_df <- model_df %>%
  mutate(coef_sd = ifelse(is.na(coef_sd), 0, coef_sd))

summary(model_df)

# convert varibles into correct types
model_df <- model_df %>%
  mutate(
    Dataset = factor(Dataset),
    PoolSize = factor(PoolSize),
    Metabolite_meta = factor(Metabolite_meta)
  )

# MOdel it

library(lme4)

model <- lmer(
  sig_freq ~ ref_abs_effect + coef_sd + PoolSize + N_remaining +
    (1 | Dataset) + (1 | Metabolite_meta),
  data = model_df
)

summary(model)
