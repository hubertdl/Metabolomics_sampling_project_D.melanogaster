# ============================================================
# Functional Module Analysis and Signal Preservation
# Across Downsampling Conditions
#
# This script assigns metabolites to functional modules and 
# evaluates the biological interpretation of diet-associated 
# metabolomic signals, as well as their preservation under 
# replicate downsampling and varying pool sizes.
#
# Workflow:
# 1. Import combined linear model results from downsampling 
#    analyses across all LC-MS panels.
# 2. Assign metabolites to functional modules based on curated 
#    biochemical classifications and pathway context.
# 3. Define a reference dataset (PoolSize = 100, full replication) 
#    to identify significant diet-associated metabolites.
# 4. Export a supplementary table linking metabolites to 
#    functional modules and identifiers.
# 5. Identify top diet-responsive metabolites within each module 
#    based on effect size and statistical significance.
# 6. Generate a module-faceted dot plot showing direction, 
#    magnitude, and significance of diet effects.
# 7. Quantify module-level signal preservation by calculating 
#    the proportion of reference metabolites recovered across 
#    replicate downsampling and pool size conditions.
# 8. Generate heatmaps showing the percentage of metabolites 
#    retained within each functional module across sampling 
#    conditions.
# 9. Combine plots into a multi-panel figure for visualization 
#    of biological signal and its robustness to sampling design.
#
# This analysis links statistical signal detection to biological 
# interpretation, demonstrating how sampling design influences 
# the preservation of functional metabolic pathways.
# ============================================================

# set working directory
setwd("")

# =========================
# Load packages
# =========================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# =========================
# Clear environment
# =========================
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

# =========================
# Load data
# =========================
results_df <- readRDS("./LM_downsample_results_META.rds") # requires large output file from downsampling analysis

# =========================
# Constants
# =========================
module_order <- c(
  "Fatty acid oxidation overload",
  "Lipogenesis and lipid storage",
  "Membrane and sphingolipid remodeling",
  "TCA cycle and central carbon metabolism",
  "Sugar metabolism and PPP",
  "Fatty acid overflow and omega-oxidation",
  "Redox imbalance and oxidative stress",
  "One-carbon and methylation metabolism",
  "Nucleotide turnover and energy balance",
  "Amino acid and nitrogen metabolism"
)

pool_order <- c("100", "50", "5")

# =========================
# Helper: assign functional modules
# =========================
assign_module <- function(df) {
  df %>%
    mutate(
      Module = case_when(
        # Acylcarnitines / FAO
        str_detect(Metabolite_meta, regex("carnitine", ignore_case = TRUE)) |
          Metabolism == "Carnitine" ~ "Fatty acid oxidation overload",
        
        # Lipid storage / synthesis
        Metabolite_meta %in% c(
          "Oleoyl glycerol",
          "Glycerol 1-myristate",
          "Hexadecanoic acid (Palmitic acid)",
          "Hexadecenoic acid",
          "Stearic acid (FA 18:0)",
          "Oleic acid (FA 18:1)",
          "Linoleic acid (FA 18:2)",
          "alpha Linolenic acid (FA 18:3)",
          "Tetradecanoic acid (Myristic acid)",
          "Decanoic acid (Capric acid)",
          "Dodecanoic acid (Lauric acid)",
          "Glycerol 3-phosphate",
          "Acetyl coA"
        ) ~ "Lipogenesis and lipid storage",
        
        # Choline / membrane / sphingolipid
        Metabolite_meta %in% c(
          "Sphinganine",
          "phosphocholine",
          "Choline",
          "Glycerophosphocholine",
          "Cytidine diphosphate choline (CDPcholine)",
          "CDP-Ethanolamine",
          "Glycerylphosphorylethanolamine"
        ) ~ "Membrane and sphingolipid remodeling",
        
        # TCA / central carbon
        Metabolite_meta %in% c(
          "Citric acid (Citrate)",
          "cis-Aconitic acid",
          "Isocitric acid (Isocitrate)",
          "2-Oxoglutaric acid (Ketoglutaric acid)",
          "Succinic acid",
          "Fumaric acid",
          "L-Malic acid (L-Malate)",
          "Oxaloacetate",
          "oxalosuccinic acid",
          "Phosphoenolpyruvic acid",
          "pyruvate"
        ) ~ "TCA cycle and central carbon metabolism",
        
        # PPP / sugar handling
        Metabolite_meta %in% c(
          "glucose",
          "glucose 6 phosphate",
          "fructose 6 phosphate",
          "fructose 1, 6 phosphate",
          "6-phosphoglucono-D-lactone",
          "6-Phosphogluconic acid",
          "D-Ribose 5-phosphate",
          "Erythrose-4-Phosphate",
          "Sedoheptulose 7-phosphate",
          "seduheptulose",
          "5-Dehydro-D-gluconate",
          "D-Gluconic acid"
        ) ~ "Sugar metabolism and PPP",
        
        # Overflow oxidation / dicarboxylic acids
        Metabolite_meta %in% c(
          "Suberic acid",
          "Azelaic Acid",
          "10-Hydroxydecanoic acid",
          "12-Hydroxydodecanoic acid ",
          "Glutaric acid",
          "3-Methyl-3-Hydroxyglutaric Acid",
          "N-Methyl-L-glutarate"
        ) ~ "Fatty acid overflow and omega-oxidation",
        
        # Redox stress
        Metabolite_meta %in% c(
          "NADP", "NADPH", "NADH", "NAD",
          "Glutathione", "Oxidized glutathione", "Ophthalmic acid",
          "Pyroglutamic acid (5-oxoproline)", "FAD", "FADH",
          "Niacinamide (Nicotinamide)", "Nicotinamide mononucleotide",
          "2-Hydroxyisobutyrate/2-Hydroxybutyrate"
        ) ~ "Redox imbalance and oxidative stress",
        
        # Methyl / one-carbon
        Metabolite_meta %in% c(
          "Methionine",
          "5'-Methylthioadenosine",
          "5'-Deoxyadenosine",
          "N6,N6,N6-Trimethyl-L-lysine",
          "S-Adenosylhomocysteine"
        ) ~ "One-carbon and methylation metabolism",
        
        # Nucleotide turnover
        str_detect(Metabolism, regex("nucleotide|purine|pyrimidine", ignore_case = TRUE)) ~
          "Nucleotide turnover and energy balance",
        
        # Amino acids / nitrogen
        str_detect(Metabolism, regex("amino acid|urea cycle|imidazole", ignore_case = TRUE)) ~
          "Amino acid and nitrogen metabolism",
        
        TRUE ~ "Other"
      )
    )
}

# =========================
# Shared preprocessing
# =========================
full_n <- max(results_df$N_remaining, na.rm = TRUE)

results_df <- results_df %>%
  mutate(
    PoolSize = as.character(PoolSize),
    N_remaining = as.numeric(as.character(N_remaining))
  ) %>%
  assign_module()

# =========================
# Reference dataset
# =========================
reference_df <- results_df %>%
  filter(
    PoolSize == "100",
    N_remaining == full_n
  )

#write.csv(reference_df, "Reference_LM_results.csv", row.names = FALSE)

# =========================
# Export module table
# =========================

module_table <- results_df %>%
  select(Metabolite_meta, Metabolism, HMDB.ID) %>%  # include HMDB here if already present
  distinct() %>%
  assign_module() %>%
  filter(Module != "Other")

write.csv(
  module_table,
  "Supplementary_Module_Metabolites.csv",
  row.names = FALSE
)

# =========================
# Figure 1: module-faceted dot plot
# =========================
dot_df <- reference_df %>%
  mutate(
    abs_coef = abs(Diet_coef),
    logFDR = -log10(Diet_FDR),
    Direction = ifelse(Diet_coef > 0, "Higher in HSD", "Lower in HSD")
  ) %>%
  filter(
    Diet_FDR < 0.05,
    Module != "Other"
  ) %>%
  mutate(score = abs_coef * logFDR) %>%
  group_by(Module) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  mutate(
    Module = factor(Module, levels = module_order)
  ) %>%
  group_by(Module) %>%
  mutate(Metabolite_plot = reorder(Metabolite_meta, Diet_coef)) %>%
  ungroup()

fig_module_dot <- ggplot(dot_df, aes(x = Diet_coef, y = Metabolite_plot)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(aes(size = logFDR, fill = Direction), shape = 21, color = "black", stroke = 0.3) +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Higher in HSD" = "#C23B22", "Lower in HSD" = "#2C7FB8")) +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(3, 9)) +
  labs(
    x = "Diet effect size (HSD vs STD coefficient)",
    y = NULL,
    fill = NULL,
    title = "Top diet-responsive metabolites grouped by functional module"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 6, color = NA)
    )
  )

fig_module_dot

ggsave(
  "Reference_module_dotplot.png",
  fig_module_dot,
  width = 15,
  height = 10,
  dpi = 400
)

# =========================
# Figure 2: module preservation heatmap
# =========================
reference_hits <- results_df %>%
  filter(
    PoolSize == "100",
    N_remaining == full_n,
    Diet_FDR < 0.05,
    Module != "Other"
  ) %>%
  distinct(Metabolite_meta, Module)

module_counts <- reference_hits %>%
  count(Module, name = "n")

module_order_labeled <- module_counts %>%
  mutate(Module_label = paste0(Module, " (n = ", n, ")")) %>%
  arrange(match(Module, module_order)) %>%
  pull(Module_label)

recovery_df <- results_df %>%
  filter(Module != "Other") %>%
  inner_join(reference_hits, by = c("Metabolite_meta", "Module")) %>%
  mutate(Recovered_iter = Diet_FDR < 0.05) %>%
  group_by(PoolSize, N_remaining, Module, Metabolite_meta) %>%
  summarise(
    Recovery_freq = mean(Recovered_iter, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(PoolSize, N_remaining, Module) %>%
  summarise(
    Percent_recovered = 100 * mean(Recovery_freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(module_counts, by = "Module") %>%
  mutate(
    Module_label = paste0(Module, " (n = ", n, ")"),
    Module_label = factor(Module_label, levels = rev(module_order_labeled)),
    PoolSize = factor(PoolSize, levels = pool_order),
    N_remaining = factor(N_remaining, levels = sort(unique(N_remaining), decreasing = TRUE))
  )

fig_module_heatmap <- ggplot(recovery_df, aes(x = N_remaining, y = Module_label, fill = Percent_recovered)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(
    ~ PoolSize,
    nrow = 1,
    labeller = labeller(PoolSize = function(x) paste0("Pool size ", x))
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08306B",
    limits = c(0, 100),
    name = "% recovered"
  ) +
  labs(
    x = "Replicates remaining",
    y = NULL,
    title = "Preservation of metabolic modules across replicate downsampling and pool size",
    subtitle = "Mean percentage of iterations in which reference metabolites were recovered"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
  )


fig_module_heatmap

ggsave(
  "Module_preservation_HeatMap.png",
  fig_module_heatmap,
  width = 15,
  height = 10,
  dpi = 400
)


# multipanel figure

# remove titles
fig_module_dot_clean <- fig_module_dot +
  labs(title = NULL)

fig_module_heatmap_clean <- fig_module_heatmap +
  labs(title = NULL)

# combine with panel labels
fig_combined <- fig_module_dot_clean + fig_module_heatmap_clean +
  plot_layout(ncol = 2, widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")

fig_combined

ggsave(
  "Combined_Figure.tif",
  fig_combined,
  width = 25,
  height = 10,
  dpi = 400
)
