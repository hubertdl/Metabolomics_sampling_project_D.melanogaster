# ============================================================
# Signal detection analysis — Experiment 2 (HSD vs STD diet)
#
# Experimental design:
#   Diet       : HSD (high sugar) vs STD (standard)  [fixed]
#   Pool size  : 5, 50, 100 individuals               [key factor]
#   Population : 1–8 per diet (independent replicates)[random]
#   pop_id     : Diet_Population (e.g. HSD_1, STD_3)  [unique ID]
#   Panel      : C18_neg  ( 62 metabolites)
#                C18_pos  ( 69 metabolites)
#                Hilic_neg ( 52 metabolites)
#                Hilic_pos ( 18 metabolites)
#
# All flies same inbred background, same age — diet is the
# only biological variable of interest.
# 48 samples total, no missing combinations.
#
# Core question:
#   Does pool size affect the ability to detect real diet-
#   related metabolomic changes?
#
# Analysis structure:
#   Section A — Per-metabolite diet signal detection
#               For each metabolite × pool size, fit:
#                 value ~ Diet + (1|pop_id)
#               pop_id = Diet_Population, random intercept
#               absorbs between-replicate variation within diet.
#               With 8 replicates per diet the random effect is
#               properly estimable (unlike Exp 1 with 3 reps).
#
#               FDR correction applied within each pool size
#               separately. n=100 treated as ground truth.
#               NOTE: input data are log-transformed and
#               metabolite-centered (mean of each metabolite
#               subtracted across all samples). This preserves
#               absolute differences between samples for each
#               metabolite, making it appropriate for detecting
#               diet effects on individual metabolite levels.
#               Classify diet effects at n=5 and n=50 as:
#                 True positive  : FDR sig at n=100, sig + same dir
#                 False negative : FDR sig at n=100, not sig
#                 False positive : not sig at n=100, sig at smaller n
#                 True negative  : not sig at either
#
#               Summarise as sensitivity and FDR per panel × pool.
#
#   Section B — Effect size attenuation
#               For metabolites with FDR-significant diet effect
#               at n=100, compare effect sizes (Diet coefficient)
#               at n=5 and n=50 against n=100.
#               Faceted by panel, coloured by pool size comparison.
#
#   Section C — Supplementary: per-metabolite pool-size
#               sensitivity (|t| for log(n) term from:
#               value ~ log(PoolSize) + Diet + (1|pop_id))
#
#   Section D — Save figures and supplementary tables
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(patchwork)
library(writexl)

# -- 0. Configuration -----------------------------------------
DATA_DIR <- ""
OUT_DIR  <- ""

PANELS <- c("C18_neg", "C18_pos", "Hilic_neg", "Hilic_pos")

PANEL_LABELS_SHORT <- c(
  C18_neg   = "C18 Neg",
  C18_pos   = "C18 Pos",
  Hilic_neg = "HILIC Neg",
  Hilic_pos = "HILIC Pos"
)

N_LEVELS  <- c("5", "50", "100")
META_COLS <- c("Sample", "Diet", "Population", "PoolSize")

DIET_COLS  <- c("HSD" = "#D55E00", "STD" = "#0072B2")
POOL_COLS  <- c("5" = "#E69F00", "50" = "#56B4E9", "100" = "#009E73")
PANEL_COLS <- c(
  "C18_neg"   = "#E69F00",
  "C18_pos"   = "#56B4E9",
  "Hilic_neg" = "#009E73",
  "Hilic_pos" = "#CC79A7"
)


# -- 1. Read & prepare all panels -----------------------------
read_panel <- function(panel_name, data_dir) {
  path <- file.path(data_dir,
                    paste0(panel_name, "_normalized_metabolite.csv"))
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      Diet       = factor(Diet,     levels = c("STD", "HSD")),
      Population = factor(Population),
      PoolSize   = factor(PoolSize, levels = N_LEVELS),
      panel      = panel_name,
      # Unique population ID — HSD_1 and STD_1 are unrelated
      pop_id     = factor(paste(Diet, Population, sep = "_"))
    ) %>%
    select(-any_of(c("RunOrder", "Tube", "Cell", "ID"))) %>%
    select(all_of(META_COLS), pop_id, panel, everything())
}

panels_list <- map(PANELS, read_panel, data_dir = DATA_DIR) %>%
  set_names(PANELS)

cat("Panels loaded:\n")
walk2(panels_list, PANELS, ~ {
  met_n <- ncol(.x) - length(META_COLS) - 2
  cat(sprintf("  %-12s : %2d samples x %3d metabolites\n",
              .y, nrow(.x), met_n))
})
cat("\n")


# ============================================================
# SECTION A — Per-metabolite diet signal detection
# ============================================================
# Model: value ~ Diet + (1|pop_id)
#
# pop_id (e.g. HSD_1, STD_3) is the unit of replication.
# Random intercept absorbs between-population variance so
# the Diet fixed effect reflects consistent diet differences
# across replicate populations.
#
# With 16 observations per pool size (8 pops × 2 diets),
# the random effect is properly estimable.
#
# FDR correction applied within each pool size so that the
# ground truth at n=100 is conservative and reliable.

fit_diet_models <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "pop_id", "panel"))
  
  cat(sprintf("\n--- %s: fitting diet models (%d metabolites) ---\n",
              PANEL_LABELS_SHORT[panel_name], length(met_cols)))
  
  map_dfr(met_cols, function(met) {
    map_dfr(N_LEVELS, function(n_val) {
      
      d <- dat %>%
        filter(PoolSize == n_val) %>%
        select(all_of(META_COLS), pop_id, value = all_of(met)) %>%
        filter(!is.na(value))
      
      if (length(unique(d$Diet)) < 2 || nrow(d) < 6) return(NULL)
      
      tryCatch({
        # Each pop_id has exactly one observation per pool size,
        # so a random intercept per pop_id is always singular.
        # The correct model treats each population as an
        # independent replicate and tests the Diet fixed effect
        # across those 16 observations (8 per diet).
        # Population is NOT included as a covariate because it
        # is perfectly nested within Diet — including it would
        # absorb all between-population variance leaving no
        # residual df for the Diet test.
        fit <- lm(value ~ Diet, data = d)
        
        coefs    <- summary(fit)$coefficients
        diet_row <- grep("^DietHSD", rownames(coefs), value = TRUE)[1]
        if (is.na(diet_row)) return(NULL)
        
        tibble(
          metabolite    = met,
          panel         = panel_name,
          pool_size     = n_val,
          diet_estimate = coefs[diet_row, "Estimate"],
          diet_se       = coefs[diet_row, "Std. Error"],
          diet_t        = coefs[diet_row, "t value"],
          diet_p        = coefs[diet_row, "Pr(>|t|)"]
        )
      },
      error = function(e) NULL)
    })
  })
}

cat(strrep("=", 60), "\n")
cat("SECTION A — Per-metabolite diet signal detection\n")
cat(strrep("=", 60), "\n")

model_results <- map_dfr(PANELS, function(pn) {
  fit_diet_models(panels_list[[pn]], pn)
}) %>%
  mutate(
    panel     = factor(panel,     levels = PANELS),
    pool_size = factor(pool_size, levels = N_LEVELS)
  )

# FDR correction within each panel × pool size
model_results <- model_results %>%
  group_by(panel, pool_size) %>%
  mutate(diet_p_adj = p.adjust(diet_p, method = "fdr")) %>%
  ungroup()

cat("\nSignificant diet effects per pool size (FDR < 0.05):\n")
model_results %>%
  group_by(pool_size, panel) %>%
  summarise(
    n_tested = n(),
    n_sig    = sum(diet_p_adj < 0.05, na.rm = TRUE),
    pct_sig  = round(100 * n_sig / n_tested, 1),
    .groups  = "drop"
  ) %>%
  print(n = Inf)


# -- Classify concordance relative to n=100 -------------------
classify_concordance <- function(model_results) {
  
  ref <- model_results %>%
    filter(pool_size == "100") %>%
    select(metabolite, panel,
           ref_estimate = diet_estimate,
           ref_p_adj    = diet_p_adj)
  
  model_results %>%
    filter(pool_size != "100") %>%
    select(metabolite, panel, pool_size,
           comp_estimate = diet_estimate,
           comp_p_adj    = diet_p_adj) %>%
    left_join(ref, by = c("metabolite", "panel")) %>%
    filter(!is.na(ref_p_adj), !is.na(comp_p_adj)) %>%
    mutate(
      ref_sig  = ref_p_adj  < 0.05,
      comp_sig = comp_p_adj < 0.05,
      same_dir = sign(ref_estimate) == sign(comp_estimate),
      category = case_when(
        ref_sig  & comp_sig & same_dir  ~ "True positive",
        ref_sig  & comp_sig & !same_dir ~ "False positive",
        ref_sig  & !comp_sig            ~ "False negative",
        !ref_sig & comp_sig             ~ "False positive",
        !ref_sig & !comp_sig            ~ "True negative"
      ),
      category = factor(category,
                        levels = c("True positive", "False negative",
                                   "False positive", "True negative"))
    )
}

concordance_all <- classify_concordance(model_results)

cat("\nConcordance summary (diet effect, FDR-corrected):\n")
concordance_all %>%
  group_by(panel, pool_size, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(panel, pool_size) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print(n = Inf)


# -- Compute sensitivity and FDR ------------------------------
compute_sens_fdr <- function(concordance_all) {
  concordance_all %>%
    group_by(panel, pool_size, category) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = category, values_from = n,
                values_fill = 0) %>%
    mutate(
      across(any_of(c("True positive", "False negative",
                      "False positive", "True negative")),
             as.numeric),
      TP = `True positive`,
      FN = `False negative`,
      FP = `False positive`,
      sensitivity = ifelse((TP + FN) > 0,
                           TP / (TP + FN), NA_real_),
      fdr         = ifelse((FP + TP) > 0,
                           FP / (FP + TP), NA_real_),
      panel_label = factor(PANEL_LABELS_SHORT[as.character(panel)],
                           levels = PANEL_LABELS_SHORT),
      pool_label  = factor(paste0("n = ", pool_size),
                           levels = c("n = 5", "n = 50"))
    )
}

# Named colour vector keyed by short label for plotting
PANEL_COLS_LABEL <- setNames(PANEL_COLS, PANEL_LABELS_SHORT)

sens_fdr <- compute_sens_fdr(concordance_all)

cat("\nSensitivity and FDR:\n")
print(sens_fdr %>%
        select(panel, pool_size, TP, FN, FP,
               sensitivity, fdr),
      n = Inf)


# -- Plot sensitivity and FDR ---------------------------------
plot_sens_fdr <- function(sens_fdr) {
  
  dodge <- position_dodge(width = 0.25)
  
  p_sens <- ggplot(
    sens_fdr,
    aes(x      = pool_label,
        y      = sensitivity,
        colour = panel_label,
        group  = panel_label)) +
    geom_line(linewidth = 0.9, alpha = 0.7, position = dodge) +
    geom_point(size = 4, position = dodge) +
    scale_colour_manual(values = PANEL_COLS_LABEL,
                        name   = "Panel") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(title    = "Sensitivity",
         subtitle = "Proportion of true diet effects (at n=100) recovered\n(TP / (TP + FN)); significance: FDR < 0.05",
         x        = "Pool size vs n = 100",
         y        = "Sensitivity") +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(size = 9, colour = "grey40"),
          legend.position = "right")
  
  p_fdr <- ggplot(
    sens_fdr,
    aes(x      = pool_label,
        y      = fdr,
        colour = panel_label,
        group  = panel_label)) +
    geom_line(linewidth = 0.9, alpha = 0.7, position = dodge) +
    geom_point(size = 4, position = dodge) +
    scale_colour_manual(values = PANEL_COLS_LABEL,
                        name   = "Panel") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(title    = "False discovery rate",
         subtitle = "Proportion of detected effects that are spurious\n(FP / (FP + TP)); significance: FDR < 0.05",
         x        = "Pool size vs n = 100",
         y        = "False discovery rate") +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(size = 9, colour = "grey40"),
          legend.position = "right")
  
  list(sens = p_sens, fdr = p_fdr)
}

sens_fdr_plots <- plot_sens_fdr(sens_fdr)
p_sens <- sens_fdr_plots$sens
p_fdr  <- sens_fdr_plots$fdr


# ============================================================
# SECTION B — Effect size attenuation
# ============================================================
# For ALL metabolites, compare diet effect sizes at n=5 and
# n=50 against n=100. Metabolites significant at n=100
# (FDR < 0.05) are highlighted; all others shown in background.
# This gives an honest picture of effect size reliability
# across the full metabolome, not just pre-selected hits.

build_attenuation_data <- function(model_results) {
  
  ref <- model_results %>%
    filter(pool_size == "100") %>%
    select(metabolite, panel,
           ref_estimate = diet_estimate,
           ref_p_adj    = diet_p_adj) %>%
    mutate(sig_at_100 = ref_p_adj < 0.05)
  
  model_results %>%
    filter(pool_size != "100") %>%
    select(metabolite, panel, pool_size,
           comp_estimate = diet_estimate) %>%
    inner_join(ref, by = c("metabolite", "panel")) %>%
    mutate(
      panel_label = factor(PANEL_LABELS_SHORT[as.character(panel)],
                           levels = PANEL_LABELS_SHORT),
      pool_label  = factor(
        paste0("n = ", pool_size, " vs n = 100"),
        levels = c("n = 5 vs n = 100", "n = 50 vs n = 100"))
    )
}

atten_data <- build_attenuation_data(model_results)

plot_attenuation <- function(atten_data) {
  
  # Correlation computed across all metabolites
  cors <- atten_data %>%
    group_by(panel_label, pool_label) %>%
    summarise(
      r     = round(cor(ref_estimate, comp_estimate,
                        use = "complete.obs"), 2),
      r_sig = round(cor(ref_estimate[sig_at_100],
                        comp_estimate[sig_at_100],
                        use = "complete.obs"), 2),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("r = ", r, " (all)\nr = ", r_sig, " (sig)"))
  
  lim <- max(abs(c(atten_data$ref_estimate,
                   atten_data$comp_estimate)),
             na.rm = TRUE) * 1.05
  
  # Split into background (all) and foreground (significant)
  d_all <- atten_data
  d_sig <- atten_data %>% filter(sig_at_100)
  
  ggplot(d_all,
         aes(x = ref_estimate, y = comp_estimate)) +
    # Background: all metabolites, grey
    geom_point(alpha = 0.25, size = 1.5, colour = "grey60") +
    # Foreground: significant at n=100, coloured by panel
    geom_point(data = d_sig,
               aes(colour = panel_label),
               alpha = 0.75, size = 2.0) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
    geom_text(
      data = cors,
      aes(label = label, x = -Inf, y = Inf),
      hjust = -0.05, vjust = 1.3,
      inherit.aes = FALSE,
      size = 2.8, colour = "grey30"
    ) +
    scale_colour_manual(values = PANEL_COLS_LABEL,
                        name   = "Panel\n(sig. at n=100)") +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    facet_grid(panel_label ~ pool_label) +
    labs(
      x = "Diet effect size at n = 100",
      y = "Diet effect size at smaller pool size"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.text      = element_text(face = "bold"),
      legend.position = "none"
    )
}

p_attenuation <- plot_attenuation(atten_data)


# ============================================================
# SECTION C — Supplementary: per-metabolite sensitivity
# ============================================================
# Model: value ~ log(PoolSize) + Diet + (1|pop_id)
# log(PoolSize) coefficient = how strongly each metabolite
# tracks pool size after accounting for diet and replicate
# variation.

fit_sensitivity_models <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "pop_id", "panel"))
  
  map_dfr(met_cols, function(met) {
    
    d <- dat %>%
      select(all_of(META_COLS), pop_id, value = all_of(met)) %>%
      filter(!is.na(value)) %>%
      mutate(n_log = log(as.numeric(as.character(PoolSize))))
    
    if (nrow(d) < 8) return(NULL)
    
    tryCatch({
      # All 48 observations used here. Each pop_id appears 3
      # times (once per pool size), so (1|pop_id) is estimable.
      # Random intercept absorbs stable between-population
      # differences; n_log captures systematic change with pool
      # size; Diet accounts for diet effect so n_log reflects
      # pool-size sensitivity independently of diet.
      fit <- suppressWarnings(
        lmer(value ~ n_log + Diet + (1 | pop_id),
             data = d, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
      )
      # Fall back to lm if random effect is still singular
      if (isSingular(fit)) {
        fit <- lm(value ~ n_log + Diet, data = d)
      }
      coefs <- summary(fit)$coefficients
      if (!"n_log" %in% rownames(coefs)) return(NULL)
      tibble(
        metabolite = met,
        panel      = panel_name,
        estimate   = coefs["n_log", "Estimate"],
        t_value    = coefs["n_log", "t value"],
        p_value    = coefs["n_log", "Pr(>|t|)"],
        abs_t      = abs(coefs["n_log", "t value"])
      )
    },
    error = function(e) NULL)
  }) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "fdr"),
      sig   = p_adj < 0.05
    )
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SECTION C — Per-metabolite pool-size sensitivity\n")
cat(strrep("=", 60), "\n")

sensitivity_results <- map_dfr(PANELS, function(pn) {
  cat(sprintf("  %s...\n", PANEL_LABELS_SHORT[pn]))
  fit_sensitivity_models(panels_list[[pn]], pn)
}) %>%
  mutate(
    panel     = factor(panel, levels = PANELS),
    direction = ifelse(estimate > 0,
                       "Increases with n", "Decreases with n")
  )

sensitivity_results %>%
  group_by(panel) %>%
  summarise(
    n_tested     = n(),
    n_sig        = sum(sig, na.rm = TRUE),
    pct_sig      = round(100 * n_sig / n_tested, 1),
    median_abs_t = round(median(abs_t, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print()

plot_sensitivity <- function(sensitivity_results) {
  
  d <- sensitivity_results %>%
    mutate(panel_label = factor(PANEL_LABELS_SHORT[panel],
                                levels = PANEL_LABELS_SHORT))
  
  ggplot(d, aes(x = panel_label, y = abs_t)) +
    geom_violin(aes(fill = panel), alpha = 0.2, colour = NA) +
    geom_jitter(aes(colour = direction),
                width = 0.1, alpha = 0.45, size = 1.3) +
    geom_boxplot(width = 0.1, outlier.shape = NA,
                 fill = "white", alpha = 0.8, linewidth = 0.5) +
    geom_hline(yintercept = 2, linetype = "dashed",
               colour = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = PANEL_COLS, guide = "none") +
    scale_colour_manual(
      values = c("Increases with n" = "#0072B2",
                 "Decreases with n" = "#D55E00"),
      name = "Direction"
    ) +
    labs(
      title    = "Per-metabolite pool-size sensitivity",
      subtitle = "Each point = one metabolite; dashed line = |t| = 2",
      x        = "Panel",
      y        = "|t-value| for pool size effect"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 9, colour = "grey40"),
      legend.position = "right"
    )
}

p_sensitivity <- plot_sensitivity(sensitivity_results)


# ============================================================
# SECTION D — Assemble and save figures
# ============================================================
#
# Main figure — 3-panel layout:
#   Top row    : sensitivity (A) | FDR (B)
#   Bottom row : effect size attenuation (C), full width
#
# Top row: did you detect the right diet effects?
# Bottom row: were the effect sizes reliable?

fig_main <- (p_sens | p_fdr) / p_attenuation +
  plot_layout(heights = c(1, 1.8)) +
  plot_annotation(
    title      = "Diet signal detection vs pool size",
    tag_levels = "A"
  ) &
  theme(
    plot.tag          = element_text(face = "bold", size = 13),
    plot.tag.position = c(0, 1)
  )

ggsave(file.path(OUT_DIR, "Figure_1.png"),
       fig_main, width = 14, height = 14, dpi = 300)

ggsave(file.path(OUT_DIR, "Supplementary_Figure_1.png"),
       p_sensitivity, width = 9, height = 6, dpi = 300)

print(fig_main)
print(p_sensitivity)
cat("\nFigures saved to", OUT_DIR, "\n")


# ============================================================
# SECTION E — Supplementary tables
# ============================================================

supp_concordance <- concordance_all %>%
  arrange(panel, pool_size, metabolite) %>%
  rename(
    Panel         = panel,
    Pool_size     = pool_size,
    Metabolite    = metabolite,
    Category      = category
  ) %>%
  select(Panel, Pool_size, Metabolite,
         ref_estimate, ref_p_adj,
         comp_estimate, comp_p_adj,
         Category)

supp_sens_fdr <- sens_fdr %>%
  arrange(panel, pool_size) %>%
  select(panel, pool_size, TP, FN, FP,
         sensitivity, fdr) %>%
  rename(
    Panel       = panel,
    Pool_size   = pool_size,
    Sensitivity = sensitivity,
    FDR         = fdr
  )

supp_model <- model_results %>%
  arrange(panel, pool_size, diet_p_adj) %>%
  rename(
    Panel      = panel,
    Pool_size  = pool_size,
    Metabolite = metabolite,
    Estimate   = diet_estimate,
    SE         = diet_se,
    t_value    = diet_t,
    P_value    = diet_p,
    P_adj_FDR  = diet_p_adj
  ) %>%
  select(Panel, Pool_size, Metabolite,
         Estimate, SE, t_value, P_value, P_adj_FDR)

supp_sensitivity <- sensitivity_results %>%
  arrange(panel, p_adj) %>%
  rename(
    Panel       = panel,
    Metabolite  = metabolite,
    Estimate    = estimate,
    t_value     = t_value,
    P_value     = p_value,
    P_adj_FDR   = p_adj,
    Abs_t       = abs_t,
    Significant = sig
  )

write_xlsx(
  list(
    "Concordance"      = supp_concordance,
    "Sensitivity_FDR"  = supp_sens_fdr,
    "Model_results"    = supp_model,
    "Met_sensitivity"  = supp_sensitivity
  ),
  path = file.path(OUT_DIR,
                   "Supplementary_Tables_Signal_Detection.xlsx")
)

cat("Supplementary tables saved to", OUT_DIR, "\n")