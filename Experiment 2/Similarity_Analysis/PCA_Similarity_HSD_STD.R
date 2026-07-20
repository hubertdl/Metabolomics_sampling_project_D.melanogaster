# ============================================================
# Metabolomic profile similarity vs pool size — Experiment 2
#
# Experimental design:
#   Diet       : HSD (high sugar) vs STD (standard)  [fixed]
#   Pool size  : 5, 50, 100 individuals               [key factor]
#   Population : 1–8 per diet (independent replicates)[blocking]
#   pop_id     : Diet_Population (e.g. HSD_1, STD_3)  [unique pop ID]
#   Panel      : C18_neg  ( 62 metabolites)
#                C18_pos  ( 69 metabolites)
#                Hilic_neg ( 52 metabolites)
#                Hilic_pos ( 18 metabolites)
#
# NOTE: All flies are the same age and from the same inbred
# ancestral strain. Populations are nested within diet —
# HSD_1 and STD_1 are independent and share no relationship.
#
# Key question:
#   Does pool size distort metabolomic profiles, and does this
#   differ between dietary conditions (HSD vs STD)?
#
# Analysis structure:
#   Section A — PCA biplots (joint across diets — same genetic
#               background so direct comparison is valid)
#   Section B — Pairwise Euclidean distances between pool sizes
#               within each diet x population combination
#   Section C — PERMANOVA: does pool size explain metabolomic
#               variance? (vegan::adonis2, permutations blocked
#               by unique population ID)
#   Section D — Betadisper: multivariate dispersion per pool
#               size, run within each diet separately
#   Section E — Save figures and supplementary tables
# ============================================================

library(tidyverse)
library(patchwork)
library(vegan)
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

# Colour-blind safe (Okabe-Ito)
DIET_COLS <- c("HSD" = "#D55E00", "STD" = "#0072B2")
POOL_COLS <- c("5"   = "#E69F00", "50"  = "#56B4E9", "100" = "#009E73")
PANEL_COLS <- c(
  "C18_neg"  = "#E69F00",
  "C18_pos"  = "#56B4E9",
  "Hilic_neg" = "#009E73",
  "Hilic_pos" = "#CC79A7"
)


# -- 1. Read & prepare all panels -----------------------------
read_panel <- function(panel_name, data_dir) {
  path <- file.path(data_dir, paste0(panel_name, "_normalized_sample.csv"))
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      Diet       = factor(Diet, levels = c("STD", "HSD")),
      Population = factor(Population),
      PoolSize   = factor(PoolSize, levels = N_LEVELS),
      panel      = panel_name,
      # Unique population ID: populations are nested within diet,
      # so Population 1 in HSD is independent of Population 1 in
      # STD. pop_id ensures correct permutation blocking.
      pop_id     = factor(paste(Diet, Population, sep = "_"))
    ) %>%
    select(all_of(META_COLS), pop_id, panel, everything(),
           -any_of(c("RunOrder", "Tube", "Cell", "ID")))
}

panels_list <- map(PANELS, read_panel, data_dir = DATA_DIR) %>%
  set_names(PANELS)

cat("Panels loaded:\n")
walk2(panels_list, PANELS, ~ {
  met_n <- ncol(.x) - length(META_COLS) - 2  # -2 for pop_id and panel
  cat(sprintf("  %-12s : %2d samples x %3d metabolites\n",
              .y, nrow(.x), met_n))
})
cat("\n")


# ============================================================
# SECTION A -- PCA biplots (joint across diets)
# ============================================================
# PCA run jointly across both diet groups because all flies
# share the same inbred genetic background. A joint PCA allows
# direct visualisation of both diet separation and pool size
# effects in the same space.
# Data were log-transformed and within-sample scaled during
# preprocessing. Across-sample centering applied here so PCA
# captures between-sample variation rather than mean abundance
# differences across metabolites.

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
    bind_cols(dat[, c(META_COLS, "pop_id", "panel")])
  list(scores = scores, var_exp = var_exp)
}

pca_results <- map(PANELS, function(pn) {
  run_pca(panels_list[[pn]])
}) %>% set_names(PANELS)

plot_pca <- function(scores, var_exp, title_text,
                     show_legend = FALSE) {
  p <- ggplot(scores,
              aes(x = PC1, y = PC2,
                  colour = PoolSize, shape = Diet)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = POOL_COLS,
                        name   = "Pool size (n)") +
    scale_shape_manual(
      values = c("STD" = 16, "HSD" = 17),
      name   = "Diet"
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
  if (is.null(res)) return(ggplot() + theme_void())
  plot_pca(res$scores, res$var_exp,
           PANEL_LABELS_SHORT[pn],
           show_legend = TRUE)
}) %>% set_names(PANELS)

fig1 <- wrap_plots(pca_plot_list, ncol = 2, byrow = TRUE) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title      = "PCA of metabolomic profiles by pool size and diet",
    tag_levels = "A"
  ) &
  theme(
    legend.position   = "right",
    plot.tag          = element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1)
  )


# ============================================================
# SECTION B -- Pairwise distances between pool sizes
# ============================================================
# Euclidean distances in full metabolomic space computed between
# all pool-size pairs (5-100, 5-50, 50-100) within each
# diet x population combination. One distance per pair per
# replicate population. Treated as descriptive given the formal
# PERMANOVA test in Section C.

compute_pairwise_distances <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "pop_id", "panel"))
  
  dat %>%
    mutate(n_int = as.integer(as.character(PoolSize))) %>%
    group_by(Diet, Population) %>%
    group_modify(~ {
      df <- .x %>% arrange(n_int)
      ns <- df$n_int
      if (length(ns) < 2) return(tibble())
      
      X <- df[, met_cols] %>%
        mutate(across(everything(), as.numeric)) %>%
        as.matrix()
      
      pair_combns <- combn(seq_len(nrow(df)), 2, simplify = FALSE)
      
      map_dfr(pair_combns, function(idx) {
        n1 <- ns[idx[1]]; n2 <- ns[idx[2]]
        v1 <- as.numeric(X[idx[1], ])
        v2 <- as.numeric(X[idx[2], ])
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
cat("SECTION B — PAIRWISE DISTANCES\n")
cat("Euclidean distances in full metabolomic space\n")
cat(strrep("=", 60), "\n")

pairwise_dist_results <- map_dfr(PANELS, function(pn) {
  compute_pairwise_distances(panels_list[[pn]], pn)
}) %>%
  mutate(panel = factor(panel, levels = PANELS))

print(
  pairwise_dist_results %>%
    group_by(panel, Diet, pair) %>%
    summarise(mean_dist = round(mean(dist), 2),
              se_dist   = round(sd(dist) / sqrt(n()), 2),
              .groups   = "drop"),
  n = Inf
)

plot_pairwise_distances <- function(pairwise_dist_results) {
  
  summ <- pairwise_dist_results %>%
    group_by(panel, Diet, pair) %>%
    summarise(
      mean_d = mean(dist, na.rm = TRUE),
      se_d   = sd(dist, na.rm = TRUE) / sqrt(sum(!is.na(dist))),
      .groups = "drop"
    ) %>%
    mutate(
      panel_label = factor(PANEL_LABELS_SHORT[panel],
                           levels = PANEL_LABELS_SHORT),
      Diet = factor(Diet, levels = c("STD", "HSD"))
    )
  
  raw <- pairwise_dist_results %>%
    mutate(
      panel_label = factor(PANEL_LABELS_SHORT[panel],
                           levels = PANEL_LABELS_SHORT),
      Diet = factor(Diet, levels = c("STD", "HSD"))
    )
  
  ggplot(summ,
         aes(x = pair, y = mean_d, colour = Diet)) +
    geom_errorbar(
      aes(ymin = mean_d - se_d, ymax = mean_d + se_d),
      width    = 0.15, linewidth = 0.7,
      position = position_dodge(width = 0.4)
    ) +
    geom_point(
      data = raw,
      aes(x = pair, y = dist, colour = Diet),
      position = position_jitterdodge(jitter.width = 0.08,
                                      dodge.width  = 0.4),
      alpha = 0.35, size = 1.5
    ) +
    geom_point(size = 3.5,
               position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = DIET_COLS,
                        name   = "Diet") +
    facet_wrap(~ panel_label, ncol = 2, scales = "free_y") +
    labs(
      title    = "Pairwise metabolomic distances between pool sizes",
      subtitle = "Euclidean distance in full metabolomic space (mean \u00b1 SE); points = individual populations",
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
# SECTION C -- PERMANOVA
# ============================================================
# For each panel, test whether pool size and diet explain
# significant proportions of total metabolomic variance.
# Model: dist ~ PoolSize + Diet  (marginal Type III tests)
# Permutations (999) blocked by pop_id (Diet_Population) to
# respect independence structure: populations are nested within
# diet, so HSD_1 and STD_1 are unrelated blocking units.

run_permanova <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "pop_id", "panel"))
  
  X <- dat[, met_cols] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  X <- X[, apply(X, 2, var, na.rm = TRUE) > 0]
  
  dist_mat <- dist(X, method = "euclidean")
  
  set.seed(42)
  perm <- how(nperm = 999, blocks = dat$pop_id)
  
  result <- adonis2(
    dist_mat ~ PoolSize + Diet,
    data         = dat,
    permutations = perm,
    by           = "margin"
  )
  
  cat("\n--- Panel:", PANEL_LABELS_SHORT[panel_name], "---\n")
  print(result)
  
  as_tibble(result, rownames = "term") %>%
    filter(term %in% c("PoolSize", "Diet")) %>%
    mutate(panel = panel_name) %>%
    rename(df = Df, SS = SumOfSqs, R2 = R2,
           F_stat = F, p_value = `Pr(>F)`)
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SECTION C — PERMANOVA\n")
cat("dist ~ PoolSize + Diet  (marginal, blocked by pop_id)\n")
cat(strrep("=", 60), "\n")

permanova_results <- map_dfr(PANELS, function(pn) {
  run_permanova(panels_list[[pn]], pn)
}) %>%
  mutate(
    panel = factor(panel, levels = PANELS),
    term  = factor(term, levels = c("PoolSize", "Diet"))
  )

print(permanova_results, n = Inf)

plot_permanova <- function(permanova_results) {
  
  d <- permanova_results %>%
    mutate(
      term_label = recode(as.character(term),
                          "PoolSize" = "Pool size",
                          "Diet"     = "Diet"),
      term_label = factor(term_label,
                          levels = c("Pool size", "Diet")),
      panel_label = PANEL_LABELS_SHORT[as.character(panel)]
    )
  
  ggplot(d, aes(x = panel_label, y = R2, fill = panel)) +
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
# SECTION D -- Betadisper
# ============================================================
# Multivariate dispersion per pool size, run separately within
# each diet group. This ensures distances to centroid reflect
# variability among replicate populations at each pool size,
# not differences between diets.

run_betadisper <- function(dat, panel_name) {
  
  met_cols <- setdiff(names(dat), c(META_COLS, "pop_id", "panel"))
  
  subgroups <- dat %>%
    group_by(Diet) %>%
    group_split()
  
  map_dfr(subgroups, function(sub) {
    
    diet_label <- as.character(sub$Diet[1])
    
    X <- sub[, met_cols] %>%
      mutate(across(everything(), as.numeric)) %>%
      as.matrix()
    X <- X[, apply(X, 2, var, na.rm = TRUE) > 0]
    
    if (length(unique(sub$PoolSize)) < 2) return(NULL)
    
    dist_mat <- dist(X, method = "euclidean")
    bd       <- betadisper(dist_mat, group = sub$PoolSize)
    
    set.seed(42)
    perm_test <- permutest(bd, permutations = 999)
    
    cat(sprintf("\n--- Panel: %s | Diet: %s ---\n",
                PANEL_LABELS_SHORT[panel_name], diet_label))
    cat("Betadisper permutation test:\n")
    print(perm_test)
    
    tibble(
      panel         = panel_name,
      Diet          = diet_label,
      Sample        = sub$Sample,
      Population    = sub$Population,
      PoolSize      = sub$PoolSize,
      dist_centroid = bd$distances
    )
  })
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SECTION D — BETADISPER\n")
cat("Dispersion per pool size, within each diet\n")
cat(strrep("=", 60), "\n")

betadisper_results <- map_dfr(PANELS, function(pn) {
  run_betadisper(panels_list[[pn]], pn)
}) %>%
  mutate(
    panel    = factor(panel, levels = PANELS),
    PoolSize = factor(PoolSize, levels = N_LEVELS),
    Diet     = factor(Diet, levels = c("STD", "HSD"))
  )

plot_betadisper <- function(betadisper_results) {
  
  d <- betadisper_results %>%
    mutate(
      panel_label = factor(PANEL_LABELS_SHORT[panel],
                           levels = PANEL_LABELS_SHORT)
    )
  
  summ <- d %>%
    group_by(panel_label, Diet, PoolSize) %>%
    summarise(
      mean_d = mean(dist_centroid, na.rm = TRUE),
      se_d   = sd(dist_centroid, na.rm = TRUE) /
        sqrt(sum(!is.na(dist_centroid))),
      .groups = "drop"
    )
  
  ggplot(summ,
         aes(x = PoolSize, y = mean_d,
             colour = Diet, group = Diet)) +
    geom_line(linewidth = 1.0,
              position = position_dodge(width = 0.2)) +
    geom_errorbar(
      aes(ymin = mean_d - se_d, ymax = mean_d + se_d),
      width    = 0.15, linewidth = 0.7,
      position = position_dodge(width = 0.2)
    ) +
    geom_point(
      data = d,
      aes(x = PoolSize, y = dist_centroid, colour = Diet),
      position = position_jitterdodge(jitter.width = 0.08,
                                      dodge.width  = 0.2),
      alpha = 0.4, size = 1.5
    ) +
    geom_point(size = 3.5,
               position = position_dodge(width = 0.2)) +
    scale_colour_manual(values = DIET_COLS, name = "Diet") +
    scale_x_discrete(labels = function(x) paste0("n = ", x)) +
    facet_wrap(~ panel_label, ncol = 2, scales = "free_y") +
    labs(
      title    = "Multivariate dispersion by pool size",
      subtitle = "Mean distance to group centroid (\u00b1 SE); points = individual populations",
      x        = "Pool size (n individuals)",
      y        = "Distance to centroid"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 9, colour = "grey40"),
      strip.text      = element_text(face = "bold"),
      legend.position = "right"
    )
}

p_betadisper <- plot_betadisper(betadisper_results)


# ============================================================
# SECTION E -- Save figures and supplementary tables
# ============================================================

# Figure 1: PCA biplots (2x2 grid, one per panel)
ggsave(file.path(OUT_DIR, "Fig1_PCA_biplots.png"),
       fig1, width = 12, height = 10, dpi = 300)

# Figure 2: Pairwise distances
ggsave(file.path(OUT_DIR, "Fig2_Pairwise_distances.png"),
       p_pairwise, width = 10, height = 10, dpi = 300)

# Figure 3: PERMANOVA (A) + Betadisper (B)
fig3 <- (p_permanova / p_betadisper) +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    title      = "Pool-size effects on metabolomic profiles",
    tag_levels = "A"
  ) &
  theme(
    plot.tag          = element_text(face = "bold", size = 13),
    plot.tag.position = c(0, 1)
  )

ggsave(file.path(OUT_DIR, "Fig3_PERMANOVA_Betadisper.png"),
       fig3, width = 12, height = 14, dpi = 300)

print(fig1); print(p_pairwise); print(fig3)
cat("\nAll figures saved to", OUT_DIR, "\n")


# -- Supplementary tables -------------------------------------
supp_pairwise <- pairwise_dist_results %>%
  arrange(panel, Diet, Population, pair) %>%
  select(-n1, -n2) %>%
  rename(
    Panel          = panel,
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
  arrange(panel, Diet, PoolSize, Population) %>%
  rename(
    Panel         = panel,
    Pool_size     = PoolSize,
    Dist_centroid = dist_centroid
  )

write_xlsx(
  list(
    "Pairwise_distances" = supp_pairwise,
    "PERMANOVA"          = supp_permanova,
    "Betadisper"         = supp_betadisper
  ),
  path = file.path(OUT_DIR, "Supplementary_Tables_Similarity_Exp2.xlsx")
)

cat("Supplementary tables saved to", OUT_DIR, "\n")