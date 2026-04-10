# ============================================================
# Downsampling Analysis: Effects of Pool Size and Replication
# on Metabolite Signal Detection
#
# This script evaluates how pool size and biological replication 
# jointly influence the detection of diet-associated metabolite 
# signals using a replicate downsampling framework.
#
# Workflow:
# 1. Import normalized metabolomics datasets (mean-centered by metabolite) 
#    across four LC-MS panels (C18 positive, C18 negative, HILIC positive, 
#    HILIC negative).
# 2. Define replicate populations within each diet (STD and HSD).
# 3. Systematically remove replicate populations in a balanced manner 
#    across diets (0–5 replicates removed per diet), generating all 
#    combinations of replicate removal.
# 4. For each downsampling iteration and pool size (5, 50, 100), fit 
#    per-metabolite linear models (metabolite ~ Diet).
# 5. Extract diet effect sizes, test statistics, and p-values, and apply 
#    false discovery rate (FDR) correction within each dataset × pool size 
#    × iteration.
# 6. Quantify signal detection, including the number of significant 
#    metabolites, effect size distributions, and summary statistics.
# 7. Aggregate results across iterations to assess:
#       - signal retention under reduced replication
#       - variability in effect size estimates
#       - stability of metabolite detection
# 8. Perform all computations in parallel to efficiently evaluate all 
#    replicate removal combinations.
# 9. Export full results, iteration summaries, metabolite-level stability, 
#    and aggregated summaries for downstream analysis and visualization.
#
# This analysis quantifies how reductions in pool size and biological 
# replication impact statistical power and the reliability of metabolomic 
# inference.
# ============================================================


# set workng directoy
setwd("")

# =========================
# Load packages
# =========================
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(future.apply)
library(parallel)

# =========================
# Parallel setup
# =========================
n_workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = n_workers)

# =========================
# Clear environment
# =========================
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

# =========================
# Read in normalized data
# =========================
files <- list.files(
  "../Normalization/Mean_centered_by_Metabolite",
  pattern = "\\.csv$",
  full.names = TRUE
)

data_list <- lapply(files, read.csv, header = TRUE)

names(data_list) <- gsub(
  "_normalized_metabolite",
  "",
  paste0("df_", tools::file_path_sans_ext(basename(files)))
)

names(data_list) <- gsub("C18_pos", "C18pos", names(data_list))
names(data_list) <- gsub("C18_neg", "C18neg", names(data_list))
names(data_list) <- gsub("Hilic_pos", "Hpos", names(data_list))
names(data_list) <- gsub("Hilic_neg", "Hneg", names(data_list))

# =========================
# Convert columns to factors
# Population = replicate number within diet
# =========================
data_list <- lapply(data_list, function(df) {
  df$PoolSize   <- factor(df$PoolSize, levels = c("5", "50", "100"))
  df$Population <- factor(df$Population)
  df$Diet       <- factor(df$Diet, levels = c("STD", "HSD"))
  df
})

list2env(data_list, envir = .GlobalEnv)

# =========================
# Store datasets in named list
# =========================
dataset_list <- list(
  C18_pos   = df_C18pos,
  C18_neg   = df_C18neg,
  HILIC_pos = df_Hpos,
  HILIC_neg = df_Hneg
)

# =========================
# Define replicate IDs within each diet
# Remove the same NUMBER from each diet, but not necessarily matched IDs
# =========================
std_reps <- sort(unique(as.character(df_C18pos$Population[df_C18pos$Diet == "STD"])))
hsd_reps <- sort(unique(as.character(df_C18pos$Population[df_C18pos$Diet == "HSD"])))

if (length(std_reps) != length(hsd_reps)) {
  stop("STD and HSD do not have the same number of replicate IDs.")
}

total_reps_per_diet <- length(std_reps)

# choose how many replicates to remove from EACH diet
remove_values <- 0:5

# =========================
# Function to run LMs within each pool size
# =========================
run_lm_dataset_by_pool <- function(df, dataset_name, removed_std, removed_hsd, iteration, n_remove, total_reps_per_diet) {
  
  met_cols <- 9:ncol(df)
  pool_levels <- levels(df$PoolSize)
  
  out_list <- list()
  counter <- 1
  
  for (pool in pool_levels) {
    
    df_pool <- droplevels(df %>% filter(PoolSize == pool))
    
    if (nrow(df_pool) == 0) next
    
    pool_results <- vector("list", length(met_cols))
    
    for (i in seq_along(met_cols)) {
      met <- met_cols[i]
      met_name <- colnames(df)[met]
      
      formula <- as.formula(
        paste0("`", met_name, "` ~ Diet")
      )
      
      out <- tryCatch({
        model <- lm(formula, data = df_pool)
        anova_res <- anova(model)
        coef_table <- summary(model)$coefficients
        
        diet_coef <- if ("DietHSD" %in% rownames(coef_table)) coef_table["DietHSD", "Estimate"] else NA_real_
        diet_se   <- if ("DietHSD" %in% rownames(coef_table)) coef_table["DietHSD", "Std. Error"] else NA_real_
        diet_t    <- if ("DietHSD" %in% rownames(coef_table)) coef_table["DietHSD", "t value"] else NA_real_
        diet_coef_p <- if ("DietHSD" %in% rownames(coef_table)) coef_table["DietHSD", "Pr(>|t|)"] else NA_real_
        
        mean_std <- mean(df_pool[df_pool$Diet == "STD", met_name], na.rm = TRUE)
        mean_hsd <- mean(df_pool[df_pool$Diet == "HSD", met_name], na.rm = TRUE)
        
        data.frame(
          Metabolite = met_name,
          Mean_STD = mean_std,
          Mean_HSD = mean_hsd,
          Diet_coef = diet_coef,
          Diet_abs_coef = abs(diet_coef),
          Diet_coef_SE = diet_se,
          Diet_t = diet_t,
          Diet_coef_p = diet_coef_p,
          Diet_F = anova_res["Diet", "F value"],
          Diet_p = anova_res["Diet", "Pr(>F)"],
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          Metabolite = met_name,
          Mean_STD = NA_real_,
          Mean_HSD = NA_real_,
          Diet_coef = NA_real_,
          Diet_abs_coef = NA_real_,
          Diet_coef_SE = NA_real_,
          Diet_t = NA_real_,
          Diet_coef_p = NA_real_,
          Diet_F = NA_real_,
          Diet_p = NA_real_,
          stringsAsFactors = FALSE
        )
      })
      
      pool_results[[i]] <- out
    }
    
    pool_results <- bind_rows(pool_results)
    pool_results$Diet_FDR <- p.adjust(pool_results$Diet_p, method = "fdr")
    pool_results$Diet_sig <- pool_results$Diet_FDR < 0.05
    
    pool_results$Dataset <- dataset_name
    pool_results$PoolSize <- pool
    pool_results$Iteration <- iteration
    pool_results$N_removed <- n_remove
    pool_results$N_remaining <- total_reps_per_diet - n_remove
    pool_results$Removed_STD <- paste(removed_std, collapse = ",")
    pool_results$Removed_HSD <- paste(removed_hsd, collapse = ",")
    
    out_list[[counter]] <- pool_results
    counter <- counter + 1
  }
  
  bind_rows(out_list)
}

# =========================
# Build all jobs first
# Each job = one combination of STD removals and HSD removals
# =========================
job_list <- list()
job_counter <- 1

for (n_remove in remove_values) {
  
  cat("Preparing removal level:", n_remove, "\n")
  
  if (n_remove == 0) {
    std_removed_sets <- list(character(0))
    hsd_removed_sets <- list(character(0))
  } else {
    std_removed_sets <- combn(std_reps, n_remove, simplify = FALSE)
    hsd_removed_sets <- combn(hsd_reps, n_remove, simplify = FALSE)
  }
  
  iter <- 1
  
  for (std_set in std_removed_sets) {
    for (hsd_set in hsd_removed_sets) {
      
      job_list[[job_counter]] <- list(
        n_remove = n_remove,
        iteration = iter,
        removed_std = std_set,
        removed_hsd = hsd_set
      )
      
      job_counter <- job_counter + 1
      iter <- iter + 1
    }
  }
}

cat("Total jobs:", length(job_list), "\n")

# optional: jobs per removal level
print(table(sapply(job_list, function(x) x$n_remove)))

# =========================
# Run jobs in parallel
# =========================
results_parallel <- future_lapply(
  job_list,
  function(job) {
    
    n_remove <- job$n_remove
    iter <- job$iteration
    std_set <- job$removed_std
    hsd_set <- job$removed_hsd
    
    reduced_datasets <- lapply(dataset_list, function(df) {
      if (n_remove == 0) {
        return(df)
      } else {
        droplevels(
          df %>%
            filter(
              !(Diet == "STD" & Population %in% std_set),
              !(Diet == "HSD" & Population %in% hsd_set)
            )
        )
      }
    })
    
    iter_results <- lapply(names(reduced_datasets), function(ds_name) {
      run_lm_dataset_by_pool(
        df = reduced_datasets[[ds_name]],
        dataset_name = ds_name,
        removed_std = std_set,
        removed_hsd = hsd_set,
        iteration = iter,
        n_remove = n_remove,
        total_reps_per_diet = total_reps_per_diet
      )
    })
    
    iter_results <- bind_rows(iter_results)
    
    iter_summary <- iter_results %>%
      group_by(Dataset, PoolSize, Iteration, N_removed, N_remaining, Removed_STD, Removed_HSD) %>%
      summarise(
        Diet_sig_n = sum(Diet_sig, na.rm = TRUE),
        mean_Diet_F = mean(Diet_F, na.rm = TRUE),
        median_Diet_F = median(Diet_F, na.rm = TRUE),
        mean_Diet_coef = mean(Diet_coef, na.rm = TRUE),
        median_Diet_coef = median(Diet_coef, na.rm = TRUE),
        mean_abs_Diet_coef = mean(Diet_abs_coef, na.rm = TRUE),
        median_abs_Diet_coef = median(Diet_abs_coef, na.rm = TRUE),
        .groups = "drop"
      )
    
    list(
      all_results = iter_results,
      summary_results = iter_summary
    )
  },
  future.seed = TRUE
)

# =========================
# Combine outputs
# =========================
all_results_df <- bind_rows(lapply(results_parallel, `[[`, "all_results"))
summary_results_df <- bind_rows(lapply(results_parallel, `[[`, "summary_results"))

# =========================
# Stability summary by metabolite
# =========================
metabolite_stability <- all_results_df %>%
  group_by(Dataset, PoolSize, N_removed, N_remaining, Metabolite) %>%
  summarise(
    Diet_sig_freq = mean(Diet_sig, na.rm = TRUE),
    mean_Diet_FDR = mean(Diet_FDR, na.rm = TRUE),
    mean_Diet_F = mean(Diet_F, na.rm = TRUE),
    mean_Diet_coef = mean(Diet_coef, na.rm = TRUE),
    mean_abs_Diet_coef = mean(Diet_abs_coef, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# Mean summary by removal level
# =========================
summary_by_removal <- summary_results_df %>%
  group_by(Dataset, PoolSize, N_removed, N_remaining) %>%
  summarise(
    mean_Diet_sig = mean(Diet_sig_n, na.rm = TRUE),
    sd_Diet_sig = sd(Diet_sig_n, na.rm = TRUE),
    mean_Diet_F = mean(mean_Diet_F, na.rm = TRUE),
    sd_Diet_F = sd(mean_Diet_F, na.rm = TRUE),
    mean_Diet_coef = mean(mean_Diet_coef, na.rm = TRUE),
    sd_Diet_coef = sd(mean_Diet_coef, na.rm = TRUE),
    mean_abs_Diet_coef = mean(mean_abs_Diet_coef, na.rm = TRUE),
    sd_abs_Diet_coef = sd(mean_abs_Diet_coef, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# Save results
# =========================
if (!dir.exists("DownSamplingXPool_Results")) {
  dir.create("DownSamplingXPool_Results")
}

write.csv(
  all_results_df,
  file = "DownSamplingXPool_Results/LM_population_removal_all_iterations.csv",
  row.names = FALSE
)

write.csv(
  summary_results_df,
  file = "DownSamplingXPool_Results/LM_population_removal_iteration_summary.csv",
  row.names = FALSE
)

write.csv(
  metabolite_stability,
  file = "DownSamplingXPool_Results/LM__population_removal_metabolite_stability.csv",
  row.names = FALSE
)

write.csv(
  summary_by_removal,
  file = "DownSamplingXPool_Results/LM_population_removal_overall_summary.csv",
  row.names = FALSE
)

# =========================
# Shut down parallel workers
# =========================
future::plan(future::sequential)


