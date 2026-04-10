# ============================================================
# Metabolomic profile similarity vs pool size
#
# Experimental design:
#   Strain     : ORWT (inbred) vs CRB (outbred)  [fixed]
#   Age        : 14 vs 40 days                    [fixed]
#   Pool size  : 5, 50, 100 individuals           [key factor]
#   Rep        : A, B, C (replicate populations)  [blocking]
#   pop_id     : strain_rep (e.g. ORWT_A, CRB_B)   [unique pop ID for blocking]
#   Panel      : C18_pos  ( 19 metabolites)
#                Hilic_neg (149 metabolites)
#                Hilic_pos ( 53 metabolites)
#
# NOTE: C18_neg excluded — only 2 metabolites survived QC.
#
# Analysis structure:
#   Section A — PCA biplots (visualisation)
#   Section B — PERMANOVA: does pool size explain metabolomic
#               variance? (vegan::adonis2, permutations blocked
#               by replicate)
#   Section C — Betadisper: multivariate dispersion per pool
#               size as a measure of diminishing returns
#               (vegan::betadisper)
#   Section D — Save figures and supplementary tables
# ============================================================

library(tidyverse)
library(patchwork)
library(vegan)
library(writexl)

# -- 0. Configuration -----------------------------------------
DATA_DIR <- ""

OUT_DIR  <- ""
  
PANELS <- c("C18_pos", "Hilic_neg", "Hilic_pos")

PANEL_LABELS_SHORT <- c(
  C18_pos   = "C18 Pos",
  Hilic_neg = "HILIC Neg",
  Hilic_pos = "HILIC Pos"
)

N_LEVELS <- c("5", "50", "100")
META_COLS <- c("Sample", "Strain", "Population", "PoolSize", "Age")

# Colour-blind safe (Okabe-Ito)
STRAIN_COLS <- c("ORWT" = "#D55E00", "CRB"  = "#0072B2")
POOL_COLS   <- c("5"    = "#E69F00", "50"   = "#56B4E9", "100" = "#009E73")
PANEL_COLS  <- c(
  "C18_pos"   = "#56B4E9",
  "Hilic_neg" = "#009E73",
  "Hilic_pos" = "#CC79A7"
)


# -- 1. Read & prepare all panels -----------------------------
read_panel <- function(panel_name, data_dir) {
  path <- file.path(data_dir, paste0(panel_name, "_normalized_sample.csv"))
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      Strain = factor(Strain),
      Population = factor(Population),
      Age    = factor(Age),
      PoolSize = factor(PoolSize, levels = N_LEVELS),
      panel  = panel_name,
      # Unique population ID: rep labels (A/B/C) are independent
      # within each strain and must not be treated as shared units
      # across strains. pop_id gives each replicate population a
      # globally unique identifier for permutation blocking.
      pop_id = factor(paste(Strain, Population, sep = "_"))
    ) %>%
    select(all_of(META_COLS), pop_id, panel, everything())
}

panels_list <- map(PANELS, read_panel, data_dir = DATA_DIR) %>%
  set_names(PANELS)

cat("Panels loaded:\n")
walk2(panels_list, PANELS, ~ {
  met_n <- ncol(.x) - length(META_COLS) - 1
  cat(sprintf("  %-12s : %2d samples x %3d metabolites\n",
              .y, nrow(.x), met_n))
})
cat("\n")


# ============================================================
# SECTION A -- PCA biplots (visualisation only)
# ============================================================
# PCA run separately per panel and per strain.
# Strains kept separate because their metabolomes are not
# directly comparable — a joint PCA would confound pool-size
# effects with between-strain genetic differences.
# Data were log-transformed and within-sample scaled during
# preprocessing (correcting for loading variation). Across-sample
# centering is still applied here so that PCA captures between-sample
# variation rather than mean abundance differences across metabolites.
# Scaling is left FALSE as metabolites are already on a comparable
# scale after preprocessing.

run_pca <- function(dat) {
  met_cols <- setdiff(names(dat), c(META_COLS, "pop_id", "panel"))
  X <- dat[, met_cols] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  X <- X[, colSums(is.na(X)) == 0]
  if (ncol(X) < 2) return(NULL)
  pca     <- prcomp(X, center = TRUE, scale. = FALSE)
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  scores  <- as_tibble(pca$x) %>%
    bind_cols(dat[, c(META_COLS, "panel")])
  list(scores = scores, var_exp = var_exp)
}

pca_results <- map(PANELS, function(pn) {
  list(
    ORWT = run_pca(filter(panels_list[[pn]], Strain == "ORWT")),
    CRB  = run_pca(filter(panels_list[[pn]], Strain == "CRB"))
  )
}) %>% set_names(PANELS)

plot_pca <- function(scores, var_exp, title_text,
                     show_legend = FALSE) {
  p <- ggplot(scores,
              aes(x = PC1, y = PC2, colour = PoolSize, shape = Age)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = POOL_COLS, name = "Pool size (n)") +
    scale_shape_manual(
      values = c("14" = 16, "40" = 17),
      name   = "Age (days)"
    ) +
    labs(
      title = title_text,
      x     = sprintf("PC1 (%.1f%%)", 100 * var_exp[1]),
      y     = sprintf("PC2 (%.1f%%)", 100 * var_exp[2])
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 9))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

pca_plot_list <- map(PANELS, function(pn) {
  res <- pca_results[[pn]]
  p_orwt <- if (!is.null(res$ORWT))
    plot_pca(res$ORWT$scores, res$ORWT$var_exp,
             paste0(PANEL_LABELS_SHORT[pn], " - ORWT (inbred)"),
             show_legend = TRUE)
  else ggplot() + theme_void()
  p_crb  <- if (!is.null(res$CRB))
    plot_pca(res$CRB$scores, res$CRB$var_exp,
             paste0(PANEL_LABELS_SHORT[pn], " - CRB (outbred)"),
             show_legend = TRUE)
  else ggplot() + theme_void()
  list(orwt = p_orwt, crb = p_crb)
}) %>% set_names(PANELS)

plot_grid_leg <- unlist(
  map(PANELS, ~ list(pca_plot_list[[.x]]$orwt,
                     pca_plot_list[[.x]]$crb)),
  recursive = FALSE
)

fig1 <- wrap_plots(plot_grid_leg, ncol = 2, byrow = TRUE) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title      = "PCA of metabolomic profiles by pool size and age",
    tag_levels = "A"
  ) &
  theme(
    legend.position   = "right",
    plot.tag          = element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1)
  )


# ============================================================
# SECTION B -- PERMANOVA
# ============================================================
# For each panel, test whether pool size explains a significant
# proportion of total metabolomic variance, after accounting
# for age and strain.
#
# Model: dist ~ n + strain + age  (marginal / Type III tests)
# Permutations (999) blocked by replicate population to respect
# the non-independence of samples from the same replicate.
#
# Key output: R² for pool size (n) — the proportion of total
# metabolomic variance attributable to pool size.

run_permanova <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "panel"))
  
  X <- dat[, met_cols] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  # Remove zero-variance metabolites
  X <- X[, apply(X, 2, var, na.rm = TRUE) > 0]
  
  dist_mat <- dist(X, method = "euclidean")
  
  set.seed(42)
  perm <- how(nperm = 999, blocks = dat$pop_id)
  
  result <- adonis2(
    dist_mat ~ PoolSize + Strain + Age,
    data         = dat,
    permutations = perm,
    by           = "margin"   # marginal (Type III) tests
  )
  
  cat("\n--- Panel:", PANEL_LABELS_SHORT[panel_name], "---\n")
  print(result)
  
  as_tibble(result, rownames = "term") %>%
    filter(term %in% c("PoolSize", "Strain", "Age")) %>%
    mutate(panel = panel_name) %>%
    rename(df = Df, SS = SumOfSqs, R2 = R2,
           F_stat = F, p_value = `Pr(>F)`)
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SECTION B — PERMANOVA\n")
cat("dist ~ n + strain + age  (marginal tests, blocked by pop_id)\n")
cat(strrep("=", 60), "\n")

permanova_results <- map_dfr(PANELS, function(pn) {
  run_permanova(panels_list[[pn]], pn)
}) %>%
  mutate(
    panel = factor(panel, levels = PANELS),
    term  = factor(term,  levels = c("PoolSize", "Strain", "Age"))
  )

print(permanova_results, n = Inf)


# -- B1. Plot PERMANOVA R² ------------------------------------
plot_permanova <- function(permanova_results) {
  
  d <- permanova_results %>%
    mutate(
      term_label = recode(as.character(term),
                          "PoolSize" = "Pool size",
                          "Strain"   = "Strain",
                          "Age"      = "Age"),
      term_label = factor(term_label,
                          levels = c("Pool size", "Strain", "Age"))
    )
  
  ggplot(d, aes(x = panel, y = R2, fill = panel)) +
    geom_col(width = 0.6, colour = "white") +
    geom_text(
      aes(label = ifelse(p_value < 0.001, "***",
                         ifelse(p_value < 0.01,  "**",
                                ifelse(p_value < 0.05,  "*", "ns")))),
      vjust = -0.4, size = 4
    ) +
    scale_fill_manual(values = PANEL_COLS,
                      labels = PANEL_LABELS_SHORT,
                      name   = "Panel") +
    scale_x_discrete(labels = PANEL_LABELS_SHORT) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.15))
    ) +
    facet_wrap(~ term_label, nrow = 1) +
    labs(
      title = "PERMANOVA: proportion of metabolomic variance explained",
      x     = NULL,
      y     = "R\u00b2 (variance explained)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
}

p_permanova <- plot_permanova(permanova_results)


# ============================================================
# SECTION C -- Pairwise distances between pool sizes
# ============================================================
# For each panel, compute the Euclidean distance in full
# metabolomic space between all pairs of pool sizes (5-50,
# 5-100, 50-100), within each strain x age x replicate
# combination. This directly answers: how different are
# profiles between pool sizes?
#
# Using full metabolomic space (consistent with PERMANOVA and
# betadisper) rather than PCA space.

compute_pairwise_distances <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "panel"))
  
  dat %>%
    mutate(n_int = as.integer(as.character(PoolSize))) %>%
    group_by(Strain, Age, Population) %>%
    group_modify(~ {
      df <- .x %>% arrange(n_int)
      ns <- df$n_int
      if (length(ns) < 2) return(tibble())
      
      X <- df[, met_cols] %>%
        mutate(across(everything(), as.numeric)) %>%
        as.matrix()
      
      pair_combns <- combn(seq_len(nrow(df)), 2, simplify = FALSE)
      
      map_dfr(pair_combns, function(idx) {
        n1  <- ns[idx[1]]
        n2  <- ns[idx[2]]
        v1  <- as.numeric(X[idx[1], ])
        v2  <- as.numeric(X[idx[2], ])
        tibble(
          n1   = n1,
          n2   = n2,
          pair = paste0(n1, "-", n2),
          dist = sqrt(sum((v1 - v2)^2, na.rm = TRUE))
        )
      })
    }) %>%
    ungroup() %>%
    mutate(
      panel = panel_name,
      pair  = factor(pair, levels = c("5-100", "5-50", "50-100"))
    )
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SECTION C - PAIRWISE DISTANCES\n")
cat("Euclidean distances in full metabolomic space\n")
cat(strrep("=", 60), "\n")

pairwise_dist_results <- map_dfr(PANELS, function(pn) {
  compute_pairwise_distances(panels_list[[pn]], pn)
}) %>%
  mutate(panel = factor(panel, levels = PANELS))

print(
  pairwise_dist_results %>%
    group_by(panel, pair) %>%
    summarise(mean_dist = round(mean(dist), 2),
              se_dist   = round(sd(dist) / sqrt(n()), 2),
              .groups   = "drop"),
  n = Inf
)


# -- C1. Plot pairwise distances ------------------------------
# Faceted by panel (rows) x age (cols), coloured by strain.
# Each panel row has its own free y-axis scale so differences
# within each panel are visible regardless of absolute magnitude.
plot_pairwise_distances <- function(pairwise_dist_results) {
  
  summ <- pairwise_dist_results %>%
    group_by(panel, Strain, Age, pair) %>%
    summarise(
      mean_d = mean(dist, na.rm = TRUE),
      se_d   = sd(dist, na.rm = TRUE) / sqrt(sum(!is.na(dist))),
      .groups = "drop"
    ) %>%
    mutate(
      age_label   = paste0(Age, " days"),
      age_label   = factor(age_label, levels = c("14 days", "40 days")),
      Strain      = factor(Strain, levels = c("ORWT", "CRB")),
      panel_label = factor(PANEL_LABELS_SHORT[panel],
                           levels = PANEL_LABELS_SHORT)
    )
  
  raw <- pairwise_dist_results %>%
    mutate(
      age_label   = paste0(Age, " days"),
      age_label   = factor(age_label, levels = c("14 days", "40 days")),
      Strain      = factor(Strain, levels = c("ORWT", "CRB")),
      panel_label = factor(PANEL_LABELS_SHORT[panel],
                           levels = PANEL_LABELS_SHORT)
    )
  
  ggplot(summ,
         aes(x = pair, y = mean_d,
             colour = Strain)) +
    geom_errorbar(
      aes(ymin = mean_d - se_d, ymax = mean_d + se_d),
      width    = 0.15, linewidth = 0.7,
      position = position_dodge(width = 0.4)
    ) +
    geom_point(
      data = raw,
      aes(x = pair, y = dist, colour = Strain),
      position = position_jitterdodge(jitter.width = 0.08,
                                      dodge.width  = 0.4),
      alpha = 0.35, size = 1.5
    ) +
    geom_point(size = 3.5,
               position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = STRAIN_COLS,
                        labels = c(ORWT = "ORWT (inbred)",
                                   CRB  = "CRB (outbred)"),
                        name   = "Strain") +
    facet_grid(
      panel_label ~ age_label,
      scales = "free_y"
    ) +
    labs(
      title    = "Pairwise metabolomic distances between pool sizes",
      subtitle = "Euclidean distance in full metabolomic space (mean \u00b1 SE); points = individual replicates",
      x        = "Pool-size pair",
      y        = "Euclidean distance"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 9, colour = "grey40"),
      strip.text      = element_text(face = "bold"),
      legend.position = "right"
    )
}

p_pairwise <- plot_pairwise_distances(pairwise_dist_results)


# ============================================================
# SECTION D -- Betadisper (diminishing returns)
# ============================================================
# betadisper measures the average distance of each sample from
# its group centroid in multivariate space (groups = pool size).
# A large dispersion at n=5 that shrinks at n=50 and n=100
# indicates that small pools produce more variable, less
# reproducible profiles.
#
# IMPORTANT: betadisper is run separately within each
# strain x age combination. This ensures the centroid for each
# pool-size group is computed among samples that only differ in
# pool size, not in strain or age. Without this, distances to
# centroid are dominated by strain/age differences within each
# pool-size group rather than by pool-size effects on
# reproducibility.
#
# Permutation test (permutest) checks whether dispersion
# differs significantly across pool sizes within each subgroup.

run_betadisper <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "panel"))
  
  # Run betadisper within each strain x age combination
  subgroups <- dat %>%
    group_by(Strain, Age) %>%
    group_split()
  
  map_dfr(subgroups, function(sub) {
    
    strain_label <- as.character(sub$Strain[1])
    age_label    <- as.character(sub$Age[1])
    
    X <- sub[, met_cols] %>%
      mutate(across(everything(), as.numeric)) %>%
      as.matrix()
    X <- X[, apply(X, 2, var, na.rm = TRUE) > 0]
    
    # Need at least 2 pool sizes to run betadisper
    if (length(unique(sub$PoolSize)) < 2) return(NULL)
    
    dist_mat <- dist(X, method = "euclidean")
    
    bd <- betadisper(dist_mat, group = sub$PoolSize)
    
    set.seed(42)
    perm_test <- permutest(bd, permutations = 999)
    
    cat(sprintf("\n--- Panel: %s | Strain: %s | Age: %s days ---\n",
                PANEL_LABELS_SHORT[panel_name], strain_label, age_label))
    cat("Betadisper permutation test:\n")
    print(perm_test)
    
    tibble(
      panel         = panel_name,
      Strain        = strain_label,
      Age           = age_label,
      Sample        = sub$Sample,
      Population    = sub$Population,
      PoolSize      = sub$PoolSize,
      dist_centroid = bd$distances
    )
  })
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SECTION D — BETADISPER\n")
cat("Dispersion per pool size, within strain x age subgroups\n")
cat(strrep("=", 60), "\n")

betadisper_results <- map_dfr(PANELS, function(pn) {
  run_betadisper(panels_list[[pn]], pn)
}) %>%
  mutate(
    panel = factor(panel, levels = PANELS),
    PoolSize = factor(PoolSize, levels = N_LEVELS)
  )



# -- C1. Plot betadisper distances ----------------------------
# Faceted by panel (rows) x age (cols), coloured by strain.
# Free y-axis scales so each panel's diminishing returns
# pattern is visible on its own scale.
# Lines retained here as pool size IS an ordered sequence.
plot_betadisper <- function(betadisper_results) {
  
  d <- betadisper_results %>%
    mutate(
      age_label   = paste0(Age, " days"),
      age_label   = factor(age_label, levels = c("14 days", "40 days")),
      Strain      = factor(Strain, levels = c("ORWT", "CRB")),
      panel_label = factor(PANEL_LABELS_SHORT[panel],
                           levels = PANEL_LABELS_SHORT)
    )
  
  summ <- d %>%
    group_by(panel_label, Strain, age_label, PoolSize) %>%
    summarise(
      mean_d = mean(dist_centroid, na.rm = TRUE),
      se_d   = sd(dist_centroid, na.rm = TRUE) /
        sqrt(sum(!is.na(dist_centroid))),
      .groups = "drop"
    )
  
  ggplot(summ,
         aes(x = PoolSize, y = mean_d,
             colour = Strain, group = Strain)) +
    geom_line(linewidth = 1.0,
              position = position_dodge(width = 0.2)) +
    geom_errorbar(
      aes(ymin = mean_d - se_d, ymax = mean_d + se_d),
      width    = 0.15, linewidth = 0.7,
      position = position_dodge(width = 0.2)
    ) +
    geom_point(
      data = d,
      aes(x = PoolSize, y = dist_centroid, colour = Strain),
      position = position_jitterdodge(jitter.width = 0.08,
                                      dodge.width  = 0.2),
      alpha = 0.4, size = 1.5
    ) +
    geom_point(size = 3.5,
               position = position_dodge(width = 0.2)) +
    scale_colour_manual(values = STRAIN_COLS,
                        labels = c(ORWT = "ORWT (inbred)",
                                   CRB  = "CRB (outbred)"),
                        name   = "Strain") +
    scale_x_discrete(labels = function(x) paste0("n = ", x)) +
    facet_grid(
      panel_label ~ age_label,
      scales = "free_y"
    ) +
    labs(
      title    = "Multivariate dispersion by pool size",
      subtitle = "Mean distance to group centroid (\u00b1 SE); points = individual samples",
      x        = "Pool size (n individuals)",
      y        = "Distance to centroid"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      strip.text    = element_text(face = "bold"),
      legend.position = "right"
    )
}

p_betadisper <- plot_betadisper(betadisper_results)


# ============================================================
# SECTION E -- Save figures and supplementary tables
# ============================================================

# Figure 1: PCA biplots
ggsave(file.path(OUT_DIR, "Fig1_PCA_biplots.tif"),
       fig1, width = 12, height = 10, dpi = 400)

# Figure 2: Pairwise distances between pool sizes
ggsave(file.path(OUT_DIR, "Fig2_Pairwise_distances.tif"),
       p_pairwise, width = 10, height = 11, dpi = 400)

# Figure 3: PERMANOVA (A) + Betadisper (B)
fig3 <- (p_permanova / p_betadisper) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title      = "Pool-size effects on metabolomic profiles",
    tag_levels = "A"
  ) &
  theme(
    plot.tag          = element_text(face = "bold", size = 13),
    plot.tag.position = c(0, 1)
  )

ggsave(file.path(OUT_DIR, "Fig3_PERMANOVA_Betadisper.tif"),
       fig3, width = 12, height = 14, dpi = 400)

print(fig1); print(p_pairwise); print(fig3)
cat("\nAll figures saved to", OUT_DIR, "\n")


# -- Supplementary tables -------------------------------------
supp_pairwise <- pairwise_dist_results %>%
  arrange(panel, Strain, Age, Population, pair) %>%
  select(-n1, -n2) %>%
  rename(
    Panel          = panel,
    Strain         = Strain,
    Age_days       = Age,
    Replicate      = Population,
    Pool_size_pair = pair,
    Euclidean_dist = dist
  )

supp_permanova <- permanova_results %>%
  arrange(panel, term) %>%
  rename(
    Panel   = panel,
    Term    = term,
    DF      = df,
    SS      = SS,
    R2      = R2,
    F_stat  = F_stat,
    P_value = p_value
  )

supp_betadisper <- betadisper_results %>%
  arrange(panel, Strain, Age, PoolSize, Population) %>%
  rename(
    Panel         = panel,
    Strain        = Strain,
    Age_days      = Age,
    Replicate     = Population,
    Pool_size     = PoolSize,
    Dist_centroid = dist_centroid
  )

write_xlsx(
  list(
    "Pairwise_distances" = supp_pairwise,
    "PERMANOVA"          = supp_permanova,
    "Betadisper"         = supp_betadisper
  ),
  path = file.path(OUT_DIR, "Supplementary_Tables_Similarity.xlsx")
)

cat("Supplementary tables saved to", OUT_DIR, "\n")