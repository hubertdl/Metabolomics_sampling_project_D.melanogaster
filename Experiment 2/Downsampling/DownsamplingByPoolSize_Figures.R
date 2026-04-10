# ============================================================
# Analysis and Visualization of Downsampling Results:
# Signal Detection, Effect Size Concordance, and Retention
#
# This script analyzes results from replicate downsampling and 
# generates figures quantifying how pool size and replication 
# influence metabolite signal detection and effect size stability.
#
# Workflow:
# 1. Import per-metabolite and summary results from the 
#    downsampling linear model analysis.
# 2. Format datasets and define factor levels for consistent 
#    comparisons across LC-MS panels and pool sizes.
# 3. Fit linear mixed-effects models to test the effects of 
#    replicate number and pool size on the number of detected 
#    diet-associated metabolites.
# 4. Define a reference set of significant metabolites at 
#    PoolSize = 100 with full replication.
# 5. Classify detected metabolites as true positives or false 
#    positives relative to the reference set, and quantify 
#    detection as a percentage of reference metabolites.
# 6. Generate line plots showing true and false positive 
#    detection across replicate downsampling and pool sizes.
# 7. Assess effect size concordance by comparing diet effect 
#    estimates at smaller pool sizes (5, 50) against the 
#    reference (100), including correlation analysis.
# 8. Bin metabolites into high, medium, and low effect-size 
#    categories based on the reference condition.
# 9. Quantify and visualize signal retention across 
#    downsampling conditions for each effect-size bin.
# 10. Export summary tables and publication-quality figures.
#
# This analysis demonstrates how sampling design influences 
# statistical power, signal retention, and the stability of 
# metabolite effect size estimates.
# ============================================================

# set working directory
setwd("")

# =========================
# Load packages
# =========================
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)

# =========================
# Clear environment
# =========================
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

# =========================
# Read in results files
# =========================
files <- list.files(
  "DownSamplingXPool_Results/",
  pattern = "\\.csv$",
  full.names = TRUE
)

data_list <- lapply(files, read.csv, header = TRUE)

# use filenames as list names first
names(data_list) <- tools::file_path_sans_ext(basename(files))

# rename to cleaner object names
names(data_list) <- gsub(
  "LM_population_removal_all_iterations",
  "all_results_df",
  names(data_list)
)

names(data_list) <- gsub(
  "LM_population_removal_iteration_summary",
  "summary_results_df",
  names(data_list)
)

# bring into environment
list2env(data_list, envir = .GlobalEnv)

# =========================
# Basic formatting
# =========================
dataset_levels <- c("C18_pos", "C18_neg", "HILIC_pos", "HILIC_neg")
pool_levels <- c("5", "50", "100")
n_remaining_levels <- sort(unique(all_results_df$N_remaining), decreasing = TRUE)

all_results_df$Dataset <- factor(all_results_df$Dataset, levels = dataset_levels)
all_results_df$PoolSize <- factor(all_results_df$PoolSize, levels = pool_levels)
all_results_df$N_remaining <- factor(all_results_df$N_remaining, levels = n_remaining_levels)

summary_results_df$Dataset <- factor(summary_results_df$Dataset, levels = dataset_levels)
summary_results_df$PoolSize <- factor(summary_results_df$PoolSize, levels = pool_levels)
summary_results_df$N_remaining <- factor(summary_results_df$N_remaining, levels = n_remaining_levels)

# numeric version for model
summary_results_df$N_remaining_num <- as.numeric(as.character(summary_results_df$N_remaining))

# full replicate number
full_n <- max(as.numeric(as.character(all_results_df$N_remaining)), na.rm = TRUE)

# =========================
# Test whether replicate number and pool size affect diet-signal detection
# =========================
model_power <- lmer(
  Diet_sig_n ~ N_remaining_num * PoolSize + (1 | Dataset),
  data = summary_results_df
)

summary(model_power)
anova(model_power)

# =========================
# Figure 1: True / false positives across replicate downsampling
# As percent of reference significant metabolites
# Reference = PoolSize 100 at full replicates
# =========================

# reference metabolites
reference_hits <- all_results_df %>%
  filter(
    PoolSize == "100",
    as.numeric(as.character(N_remaining)) == full_n,
    Diet_sig == TRUE
  ) %>%
  distinct(Dataset, Metabolite) %>%
  mutate(InReference = TRUE)

# number of reference metabolites per dataset
reference_counts <- reference_hits %>%
  group_by(Dataset) %>%
  summarise(
    Reference_N = n(),
    .groups = "drop"
  )

# classify significant metabolites
sig_results <- all_results_df %>%
  filter(Diet_sig == TRUE) %>%
  left_join(reference_hits, by = c("Dataset", "Metabolite")) %>%
  mutate(
    InReference = ifelse(is.na(InReference), FALSE, InReference),
    Class = ifelse(InReference, "True Positive", "False Positive")
  )

# count TP / FP per iteration
tp_fp_iter <- sig_results %>%
  group_by(Dataset, PoolSize, N_remaining, Iteration, Class) %>%
  summarise(N = n(), .groups = "drop") %>%
  left_join(reference_counts, by = "Dataset") %>%
  mutate(
    Percent = 100 * N / Reference_N
  )

# summarise across iterations
tp_fp_summary <- tp_fp_iter %>%
  group_by(Dataset, PoolSize, N_remaining, Class) %>%
  summarise(
    mean_pct = mean(Percent, na.rm = TRUE),
    q25_pct = quantile(Percent, 0.25, na.rm = TRUE),
    q75_pct = quantile(Percent, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

tp_fp_summary$Dataset <- factor(tp_fp_summary$Dataset, levels = dataset_levels)
tp_fp_summary$PoolSize <- factor(tp_fp_summary$PoolSize, levels = pool_levels)
tp_fp_summary$N_remaining <- factor(tp_fp_summary$N_remaining, levels = n_remaining_levels)
tp_fp_summary$Class <- factor(
  tp_fp_summary$Class,
  levels = c("True Positive", "False Positive")
)

fig_tp_fp_lines <- ggplot(
  tp_fp_summary,
  aes(
    x = N_remaining,
    y = mean_pct,
    group = interaction(PoolSize, Class)
  )
) +
  geom_ribbon(
    aes(
      ymin = q25_pct,
      ymax = q75_pct,
      fill = PoolSize
    ),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(
    aes(
      color = PoolSize,
      linetype = Class
    ),
    linewidth = 1
  ) +
  geom_point(
    aes(
      color = PoolSize,
      shape = Class
    ),
    size = 2.2
  ) +
  facet_wrap(~ Dataset, scales = "fixed") +
  scale_linetype_manual(
    values = c(
      "True Positive" = "solid",
      "False Positive" = "dashed"
    )
  ) +
  scale_shape_manual(
    values = c(
      "True Positive" = 16,
      "False Positive" = 1
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  theme_classic() +
  labs(
    x = "Replicates remaining per diet",
    y = "Percent of reference significant metabolites",
    color = "Pool size",
    linetype = "",
    shape = "",
    fill = "Pool size",
    title = "True and false positives across replicate downsampling",
    subtitle = "Solid = true positives, dashed = false positives; values scaled to PoolSize 100 at full replicates (= 100%)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

print(fig_tp_fp_lines)

# =========================
# Figure 3: Effect-size concordance across pool sizes
# Full replicates only
# =========================

# full replicate results only
full_results <- all_results_df %>%
  filter(as.numeric(as.character(N_remaining)) == full_n)

# significant reference metabolites at PoolSize 100 full replicates
reference_sig <- full_results %>%
  filter(
    PoolSize == "100",
    Diet_sig == TRUE
  ) %>%
  distinct(Dataset, Metabolite) %>%
  mutate(IsReferenceSig = TRUE)

# build pairwise comparison tables
ref_100 <- full_results %>%
  filter(PoolSize == "100") %>%
  select(
    Dataset, Metabolite,
    Effect_100 = Diet_coef
  )

comp_5 <- full_results %>%
  filter(PoolSize == "5") %>%
  select(
    Dataset, Metabolite,
    Effect_small = Diet_coef
  ) %>%
  left_join(ref_100, by = c("Dataset", "Metabolite")) %>%
  mutate(Comparison = "n = 5 vs n = 100")

comp_50 <- full_results %>%
  filter(PoolSize == "50") %>%
  select(
    Dataset, Metabolite,
    Effect_small = Diet_coef
  ) %>%
  left_join(ref_100, by = c("Dataset", "Metabolite")) %>%
  mutate(Comparison = "n = 50 vs n = 100")

effect_compare_df <- bind_rows(comp_5, comp_50) %>%
  left_join(reference_sig, by = c("Dataset", "Metabolite")) %>%
  mutate(
    IsReferenceSig = ifelse(is.na(IsReferenceSig), FALSE, IsReferenceSig)
  )

effect_compare_df$Dataset <- factor(effect_compare_df$Dataset, levels = dataset_levels)
effect_compare_df$Comparison <- factor(
  effect_compare_df$Comparison,
  levels = c("n = 5 vs n = 100", "n = 50 vs n = 100")
)

# correlation labels
cor_labels <- effect_compare_df %>%
  group_by(Dataset, Comparison) %>%
  summarise(
    r_all = cor(Effect_100, Effect_small, use = "complete.obs"),
    r_sig = cor(
      Effect_100[IsReferenceSig],
      Effect_small[IsReferenceSig],
      use = "complete.obs"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "r = ", round(r_all, 2), " (all)\n",
      "r = ", round(r_sig, 2), " (sig)"
    )
  )

# shared axis limits
axis_lim <- max(abs(c(effect_compare_df$Effect_100, effect_compare_df$Effect_small)), na.rm = TRUE)
axis_lim <- ceiling(axis_lim * 10) / 10

cor_labels <- cor_labels %>%
  mutate(
    x = -0.95 * axis_lim,
    y =  0.95 * axis_lim
  )

fig_effect_scatter <- ggplot(effect_compare_df, aes(x = Effect_100, y = Effect_small)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "gray40", linewidth = 0.6
  ) +
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "gray85", linewidth = 0.4) +
  geom_point(
    color = "gray80",
    size = 2,
    alpha = 0.8
  ) +
  geom_point(
    data = effect_compare_df %>% filter(IsReferenceSig),
    aes(color = Dataset),
    size = 2.4,
    alpha = 0.9
  ) +
  geom_text(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  facet_grid(Dataset ~ Comparison) +
  coord_fixed(
    xlim = c(-axis_lim, axis_lim),
    ylim = c(-axis_lim, axis_lim)
  ) +
  theme_classic() +
  labs(
    x = "Diet effect size at n = 100",
    y = "Diet effect size at smaller pool size",
    title = "Effect-size concordance across pool sizes",
    subtitle = "Full replicates only; colored points = metabolites significant at PoolSize 100"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )

print(fig_effect_scatter)

# =========================
# Save outputs
# =========================
if (!dir.exists("DownSamplingXPool_Figures")) {
  dir.create("DownSamplingXPool_Figures")
}

write.csv(
  reference_effect_df,
  "DownSamplingXPool_Figures/Reference_effect_size_retention_summary.csv",
  row.names = FALSE
)

write.csv(
  effect_compare_df,
  "DownSamplingXPool_Figures/EffectSize_Concordance_byPool_points.csv",
  row.names = FALSE
)

write.csv(
  cor_labels,
  "DownSamplingXPool_Figures/EffectSize_Concordance_byPool_correlations.csv",
  row.names = FALSE
)

ggsave(
  "DownSamplingXPool_Figures/DietSignal_TPtoFP_Line_byPool.tif",
  fig_tp_fp_lines,
  width = 14,
  height = 10,
  dpi = 300,
  compression = "lzw"
)


##################################
# effect size (binned) retention #
##################################


# =========================
# 1. Define reference condition
# =========================
full_n <- max(all_results_df$N_remaining, na.rm = TRUE)

reference_df <- all_results_df %>%
  filter(
    PoolSize == "100",
    N_remaining == full_n
  ) %>%
  group_by(Dataset, Metabolite) %>%
  summarise(
    Ref_Effect = mean(Diet_coef, na.rm = TRUE),
    Ref_Abs_Effect = mean(Diet_abs_coef, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# 2. Bin metabolites into high / medium / low within each dataset
# =========================
reference_bins <- reference_df %>%
  group_by(Dataset) %>%
  mutate(
    Effect_Bin = case_when(
      Ref_Abs_Effect <= quantile(Ref_Abs_Effect, 1/3, na.rm = TRUE) ~ "Low",
      Ref_Abs_Effect <= quantile(Ref_Abs_Effect, 2/3, na.rm = TRUE) ~ "Medium",
      TRUE ~ "High"
    )
  ) %>%
  ungroup()

# Put HIGH at top and LOW at bottom
reference_bins$Effect_Bin <- factor(
  reference_bins$Effect_Bin,
  levels = c("High", "Medium", "Low")
)

# =========================
# 3. Join bins back onto all downsampled results
# =========================
binned_df <- all_results_df %>%
  left_join(
    reference_bins %>% select(Dataset, Metabolite, Effect_Bin, Ref_Abs_Effect),
    by = c("Dataset", "Metabolite")
  )

binned_df$Dataset <- factor(
  binned_df$Dataset,
  levels = c("C18_pos", "C18_neg", "HILIC_pos", "HILIC_neg")
)

binned_df$PoolSize <- factor(
  binned_df$PoolSize,
  levels = c("5", "50", "100")
)

# enforce full range and order
n_levels <- sort(unique(all_results_df$N_remaining), decreasing = TRUE)

binned_df$N_remaining <- factor(
  binned_df$N_remaining,
  levels = n_levels
)

binned_df$Effect_Bin <- factor(
  binned_df$Effect_Bin,
  levels = c("High", "Medium", "Low")
)

# =========================
# 4. Detection retention by bin (PER ITERATION)
# =========================
bin_iter_df <- binned_df %>%
  group_by(Dataset, PoolSize, N_remaining, Iteration, Effect_Bin) %>%
  summarise(
    frac_sig = mean(Diet_sig, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# 5. Summarise for plotting (mean + IQR)
# =========================
bin_detection_summary <- bin_iter_df %>%
  group_by(Dataset, PoolSize, N_remaining, Effect_Bin) %>%
  summarise(
    frac_sig_mean = mean(frac_sig, na.rm = TRUE),
    frac_sig_q25 = quantile(frac_sig, 0.25, na.rm = TRUE),
    frac_sig_q75 = quantile(frac_sig, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

bin_detection_summary$Dataset <- factor(
  bin_detection_summary$Dataset,
  levels = c("C18_pos", "C18_neg", "HILIC_pos", "HILIC_neg")
)

bin_detection_summary$PoolSize <- factor(
  bin_detection_summary$PoolSize,
  levels = c("5", "50", "100")
)

bin_detection_summary$N_remaining <- factor(
  bin_detection_summary$N_remaining,
  levels = n_levels
)

bin_detection_summary$Effect_Bin <- factor(
  bin_detection_summary$Effect_Bin,
  levels = c("High", "Medium", "Low")
)

# =========================
# 6. Plot (with ribbon)
# =========================
fig_bin_detection <- ggplot(
  bin_detection_summary,
  aes(x = N_remaining, y = frac_sig_mean, group = PoolSize)
) +
  geom_ribbon(
    aes(
      ymin = frac_sig_q25,
      ymax = frac_sig_q75,
      fill = PoolSize
    ),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(aes(color = PoolSize), linewidth = 1) +
  geom_point(aes(color = PoolSize), size = 2.2) +
  facet_grid(Effect_Bin ~ Dataset) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  labs(
    x = "Replicates remaining per diet",
    y = "Fraction still significant",
    title = "Retention of metabolite signal by effect-size bin",
    subtitle = "Line = mean across iterations; ribbon = interquartile range (25–75%)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    
    # 👇 NEW
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(0.8, "lines"),
    strip.background = element_rect(fill = "white", color = "black")
  )

print(fig_bin_detection)
# =========================
# 7. Save
# =========================
write.csv(
  bin_iter_df,
  "DownSamplingXPool_Figures/EffectSizeBin_DetectionRetention_byIteration.csv",
  row.names = FALSE
)

write.csv(
  bin_detection_summary,
  "DownSamplingXPool_Figures/EffectSizeBin_DetectionRetention_summary.csv",
  row.names = FALSE
)

ggsave(
  "DownSamplingXPool_Figures/Fig_EffectSizeBin_DetectionRetention.tif",
  fig_bin_detection,
  width = 12,
  height = 9,
  dpi = 300,
  compression = "lzw"
)

