


library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(lme4)
library(lmerTest)
library(broom.mixed)

setwd("C:/Users/huber/Box/Phillips Lab/Projects/Metabolomics_Sampling_Project/PoolSize&GeneticVariation/LMM")


# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)


# ---------------------------------------------------------------
# 1. Load and stack the four panel files
#    Layout: one row per sample, with Strain/Population/PoolSize/Age
#    as columns, followed by one column per metabolite.
# ---------------------------------------------------------------
panel_files <- list.files("Mean_centered_by_Metabolite", pattern = "\\.csv$",
                          full.names = TRUE, ignore.case = TRUE)

panel_label <- function(path) {
  fname <- basename(path)
  case_when(
    str_detect(fname, regex("C18.*neg",   ignore_case = TRUE)) ~ "C18_neg",
    str_detect(fname, regex("C18.*pos",   ignore_case = TRUE)) ~ "C18_pos",
    str_detect(fname, regex("Hilic.*neg", ignore_case = TRUE)) ~ "HILIC_neg",
    str_detect(fname, regex("Hilic.*pos", ignore_case = TRUE)) ~ "HILIC_pos",
    TRUE ~ fname
  )
}

id_cols <- c("Sample", "Strain", "Population", "PoolSize", "Age")

read_panel <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  
  df %>%
    pivot_longer(-all_of(id_cols), names_to = "metabolite", values_to = "value") %>%
    mutate(panel = panel_label(path))
}

dat1 <- map_dfr(panel_files, read_panel) %>%
  mutate(
    strain_pop = paste(Strain, Population, sep = "_"),
    Age        = factor(Age, levels = c(14, 40)),   # 14 = reference
    PoolSize   = as.numeric(PoolSize)
  )

# Sanity checks
message("Rows per panel:")
print(dat1 %>% count(panel))
message("Distinct strain_pop groups (should be 6, e.g. CRB_A...ORWT_C): ",
        n_distinct(dat1$strain_pop))
message("Distinct metabolites per panel:")
print(dat1 %>% distinct(panel, metabolite) %>% count(panel))

# ---------------------------------------------------------------
# 4. Fit one LMM per metabolite, within each panel
# ---------------------------------------------------------------
fit_one <- function(d) {
  m <- tryCatch(
    lmerTest::lmer(value ~ Age + log(PoolSize) + Strain + (1 | strain_pop), data = d),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(m)) return(tibble())
  broom.mixed::tidy(m, effects = "fixed") %>%
    filter(term == "Age40")
}

age_results <- dat1 %>%
  group_by(panel, metabolite) %>%
  group_modify(~ fit_one(.x)) %>%
  ungroup() %>%
  rename(age_estimate = estimate, age_se = std.error,
         age_statistic = statistic, age_pvalue = p.value) %>%
  select(-term, -effect)

n_total <- dat1 %>% distinct(panel, metabolite) %>% nrow()
n_fit   <- nrow(age_results)
message(n_total - n_fit, " of ", n_total, " metabolite models failed to fit or converge")

# ---------------------------------------------------------------
# 5. BH-FDR correction within panel, summary, and targeted checks
# ---------------------------------------------------------------
age_results <- age_results %>%
  group_by(panel) %>%
  mutate(age_fdr = p.adjust(age_pvalue, method = "BH")) %>%
  ungroup()

age_summary <- age_results %>%
  filter(age_fdr < 0.05) %>%
  group_by(panel) %>%
  summarise(n_sig = n(), n_increase = sum(age_estimate > 0),
            n_decrease = sum(age_estimate < 0), .groups = "drop")
print(age_summary)

carnitines <- age_results %>%
  filter(str_detect(metabolite, regex("carnitine", ignore_case = TRUE)), age_fdr < 0.05)
print(carnitines)

glycolysis <- age_results %>%
  filter(str_detect(metabolite, regex("phosphate|glucose|fructose|G6P|F6P|G1P|F1P", ignore_case = TRUE)),
         age_fdr < 0.05)
print(glycolysis)

write_csv(age_results, "age_results_all_metabolites.csv")
write_csv(age_summary, "age_summary_by_panel.csv")
