# ============================================================
# Beta-Binomial GLMM: Determinants of Metabolite Signal Retention
# Under Replicate Downsampling
#
# Reads the per-metabolite, per-iteration linear model results
# produced by DownSamplingByPoolSizeLM_capped_resampling.R and
# models the probability that a metabolite's diet effect remains
# significant (FDR < 0.05) as a function of its reference effect
# size, the variability of its effect size estimate across
# iterations, pool size, and replicate number, using a
# beta-binomial GLMM with Dataset and Metabolite random intercepts.
#
# The full-replication reference cell (N_remaining = 8) is excluded
# from the regression, since it has only a single iteration and no
# variability to estimate.
#
# The pool size x replicate number interaction is tested explicitly
# via a likelihood ratio test against an additive model. Model
# coefficients are also converted to predicted retention
# probabilities and odds ratios for a "typical" metabolite, for use
# in the manuscript text.
#
# Output: model objects (model_df_betabinomial_capped_resampling.rds,
# beta_binomial_signal_retention_model_capped_resampling.rds) and
# summary tables (Retention_Predicted_Probabilities.csv,
# Retention_OddsRatios.csv), all written to this folder.
# ============================================================
# set working directory
# This script expects to be run with the working directory set to
# Experiment 2/Signal_Retention/ within the repository. It reads
# ../LM_Results/LM_population_removal_all_iterations.csv and writes
# its model output and prediction tables here.

rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

library(dplyr)
library(glmmTMB)
library(DHARMa)

# =========================
# Load data
# =========================
results_df <- read.csv("../LM_Results/LM_population_removal_all_iterations.csv", stringsAsFactors = FALSE)

results_df <- results_df %>%
  mutate(
    PoolSize    = factor(as.character(PoolSize), levels = c("5", "50", "100")),
    N_remaining = as.numeric(as.character(N_remaining))
  )

# =========================
# Reference effect size (pool size = 100, full replication)
# =========================
full_n <- max(results_df$N_remaining, na.rm = TRUE)

reference_effects <- results_df %>%
  filter(PoolSize == "100", N_remaining == full_n) %>%
  group_by(Dataset, Metabolite) %>%
  summarise(ref_abs_effect = mean(abs(Diet_coef), na.rm = TRUE), .groups = "drop")

# =========================
# Build success/trial counts
# =========================
retention_df <- results_df %>%
  mutate(sig = Diet_FDR < 0.05) %>%
  group_by(Dataset, PoolSize, N_remaining, Metabolite) %>%
  summarise(
    n_iter = n(),
    n_sig  = sum(sig, na.rm = TRUE),
    .groups = "drop"
  )

# Sanity check: n_iter should be a uniform 64 across all
# N_remaining levels except N_remaining = 8 (which is always 1).
print(
  retention_df %>%
    group_by(N_remaining) %>%
    summarise(n_iter_unique_values = n_distinct(n_iter), n_iter_example = first(n_iter))
)

# =========================
# EXCLUDE the full-replication / reference cell
# =========================
retention_df <- retention_df %>%
  filter(N_remaining != full_n)

cat("Retained N_remaining levels:", paste(sort(unique(retention_df$N_remaining)), collapse = ", "), "\n")
cat("Rows after exclusion:", nrow(retention_df), "\n")

# =========================
# Variability in effect size estimates
# =========================
variance_df <- results_df %>%
  filter(N_remaining != full_n) %>%
  group_by(Dataset, PoolSize, N_remaining, Metabolite) %>%
  summarise(coef_sd = sd(Diet_coef, na.rm = TRUE), .groups = "drop")

# =========================
# Combine into modeling dataset
# =========================
model_df <- retention_df %>%
  left_join(reference_effects, by = c("Dataset", "Metabolite")) %>%
  left_join(variance_df, by = c("Dataset", "PoolSize", "N_remaining", "Metabolite")) %>%
  mutate(
    coef_sd    = ifelse(is.na(coef_sd), 0, coef_sd),
    Dataset    = factor(Dataset),
    PoolSize   = factor(PoolSize, levels = c("100", "50", "5")),  # 100 = reference
    Metabolite = factor(Metabolite),
    n_fail     = n_iter - n_sig
  )

summary(model_df)

# ============================================================
# Beta-binomial GLMM, with PoolSize x N_remaining interaction
# ============================================================

bb_model <- glmmTMB(
  cbind(n_sig, n_fail) ~ ref_abs_effect + coef_sd + PoolSize * N_remaining +
    (1 | Dataset) + (1 | Metabolite),
  dispformula = ~ Dataset,
  data   = model_df,
  family = betabinomial(link = "logit")
)

summary(bb_model)

# ============================================================
# Explicit test of the interaction term
# ============================================================

bb_model_additive <- glmmTMB(
  cbind(n_sig, n_fail) ~ ref_abs_effect + coef_sd + PoolSize + N_remaining +
    (1 | Dataset) + (1 | Metabolite),
  dispformula = ~ Dataset,
  data   = model_df,
  family = betabinomial(link = "logit")
)

anova(bb_model_additive, bb_model)  # LRT: interaction vs. additive model

# ============================================================
# Model diagnostics
# ============================================================

sim_res <- simulateResiduals(bb_model, n = 500)
plot(sim_res)
testDispersion(sim_res)
testZeroInflation(sim_res)
plotResiduals(sim_res, form = model_df$N_remaining)
plotResiduals(sim_res, form = model_df$ref_abs_effect)

# ============================================================

# ============================================================
# Predicted retention probabilities and odds ratios
#
# The raw bb_model coefficients are on the logit scale and aren't
# directly interpretable for the manuscript text. This section
# converts them into predicted retention probabilities for a
# "typical" metabolite (median reference effect size, median
# estimate variability) across PoolSize x N_remaining, plus odds
# ratios for each fixed effect.
#
# Predictions use re.form = NA (population-level, marginalizing
# over the Dataset and Metabolite random intercepts) so they
# reflect an "average" dataset/metabolite rather than any specific
# one.
# ============================================================

# ---- typical metabolite: median reference effect size / variability ----
med_ref_abs_effect <- median(model_df$ref_abs_effect, na.rm = TRUE)
med_coef_sd         <- median(model_df$coef_sd, na.rm = TRUE)

cat("Median reference effect size:", round(med_ref_abs_effect, 4), "\n")
cat("Median effect size variability (coef_sd):", round(med_coef_sd, 4), "\n")

# ---- prediction grid ----
pred_grid <- expand.grid(
  PoolSize        = factor(c("5", "50", "100"), levels = c("100", "50", "5")),
  N_remaining     = seq(min(model_df$N_remaining), max(model_df$N_remaining)),
  ref_abs_effect  = med_ref_abs_effect,
  coef_sd         = med_coef_sd
)

pred_grid$Predicted_Retention_Prob <- predict(
  bb_model,
  newdata = pred_grid,
  type    = "response",
  re.form = NA
)

pred_grid <- pred_grid %>%
  select(PoolSize, N_remaining, Predicted_Retention_Prob) %>%
  arrange(PoolSize, N_remaining)

print(pred_grid)

# ---- odds ratios from fixed-effect coefficients ----
fe <- fixef(bb_model)$cond

# Pool size vs. Pool size = 100, at N_remaining = 5 (mid-range)
n_ref <- 5
OR_pool5_vs_100 <- exp(
  fe["PoolSize5"] + fe["PoolSize5:N_remaining"] * n_ref
)
OR_pool50_vs_100 <- exp(
  fe["PoolSize50"] + fe["PoolSize50:N_remaining"] * n_ref
)

# Odds ratio per additional replicate, by pool size (main effect + interaction)
# NB: fe[...] returns a *named* scalar, so combining with c("label" = fe[...])
# concatenates names (e.g. "Pool 5.N_remaining") instead of using "label" --
# unname() strips that so the later ["Pool 5"] lookups below actually match.
OR_per_replicate <- c(
  "Pool 5"   = exp(unname(fe["N_remaining"]) + unname(fe["PoolSize5:N_remaining"])),
  "Pool 50"  = exp(unname(fe["N_remaining"]) + unname(fe["PoolSize50:N_remaining"])),
  "Pool 100" = exp(unname(fe["N_remaining"]))
)

# Odds ratio per 0.1-unit increase in reference effect size
OR_ref_abs_effect <- exp(fe["ref_abs_effect"] * 0.1)

# Odds ratio per 0.01-unit increase in effect size variability (coef_sd)
OR_coef_sd <- exp(fe["coef_sd"] * 0.01)

odds_ratio_summary <- data.frame(
  Term = c(
    "Pool size 5 vs 100 (at N_remaining = 5)",
    "Pool size 50 vs 100 (at N_remaining = 5)",
    "Per additional replicate, Pool size 5",
    "Per additional replicate, Pool size 50",
    "Per additional replicate, Pool size 100",
    "Per 0.1-unit increase in reference effect size",
    "Per 0.01-unit increase in effect size variability (coef_sd)"
  ),
  Odds_Ratio = c(
    OR_pool5_vs_100,
    OR_pool50_vs_100,
    OR_per_replicate["Pool 5"],
    OR_per_replicate["Pool 50"],
    OR_per_replicate["Pool 100"],
    OR_ref_abs_effect,
    OR_coef_sd
  )
)

print(odds_ratio_summary)

# ============================================================
# Save output
# ============================================================
saveRDS(model_df, "model_df_betabinomial_capped_resampling.rds")
saveRDS(bb_model, "beta_binomial_signal_retention_model_capped_resampling.rds")

write.csv(
  pred_grid,
  "Retention_Predicted_Probabilities.csv",
  row.names = FALSE
)
write.csv(
  odds_ratio_summary,
  "Retention_OddsRatios.csv",
  row.names = FALSE
)
