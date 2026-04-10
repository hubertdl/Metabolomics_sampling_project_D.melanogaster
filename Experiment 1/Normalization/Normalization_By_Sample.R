# ============================================================
# PCA Analysis of Metabolomic Data by Strain, Age, and Pool Size
# (Mean-Centered by sample)

# This script performs preprocessing and principal component 
# analysis (PCA) on metabolomics data across four LC-MS panels 
# (C18 positive, C18 negative, HILIC positive, HILIC negative),
# with analyses conducted separately for each genetic background.
#
# Workflow:
# 1. Import raw metabolite abundance data and sample metadata.
# 2. Perform quality control filtering by removing metabolites 
#    with high technical variability (coefficient of variation > 0.30 
#    in pooled QC samples).
# 3. Log10-transform metabolite abundances (with pseudocount) to 
#    stabilize variance.
# 4. Normalize data by mean-centering within each sample to control 
#    for differences in overall signal intensity while preserving 
#    within-sample metabolite composition.
# 5. Split data by strain (ORWT and CRB) to avoid confounding 
#    pool-size effects with genetic background differences.
# 6. Perform PCA separately within each strain for each LC-MS panel.
# 7. Generate PCA plots with samples colored by pool size and 
#    shaped by age.
# 8. Combine PCA plots into a multi-panel figure across panels 
#    and strains.
# 9. Export normalized datasets and high-resolution PCA figures.
#
# This analysis evaluates how pool size influences metabolomic 
# profile structure within genetic backgrounds while accounting 
# for age-related variation.
# ============================================================

# set working directory
setwd("")

# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

#load packages For PCA
library(ggplot2)
library(ggvenn)
library(dbplyr)

#read in data
df <- read.csv("Data_raw.csv", header = TRUE)

meta <- read.csv("Meta.csv", header = TRUE)

meta$PoolSize <- factor(meta$PoolSize)
meta$Age <- factor(meta$Age)



# check QC for issues

qc_cols <- c("QC_01","QC_02","QC_03", "QC_04", "QC_05", "QC_06")

# calculate CV
df$QC_CV <- apply(df[, qc_cols], 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

# inspect distribution

summary(df$QC_CV)
hist(df$QC_CV, breaks = 50)

# filter for CV > .30
df_filtered <- df[df$QC_CV < 0.30, ] # 6 metabolites removed due to high CV in QC samples

# visualization of this
ggplot(df, aes(x = QC_CV)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0.3, color = "red") +
  theme_classic()

#remove QC columns from df_filtered

df_filtered <- df_filtered[, 1:(ncol(df_filtered) - 7)]

#create separate DF by run type and strain

dfs <- split(df_filtered, list(df_filtered$Column, df_filtered$Ion.Mode), drop = TRUE)

list2env(
  setNames(dfs,
           paste0("df_", gsub("\\.", "", names(dfs)))),
  envir = .GlobalEnv
)

#### HILIC Positive

# change format for normalization

library(dplyr)

met_names <- df_HILICPositive$Metabolite
abund_mat <- df_HILICPositive[, 7:ncol(df_HILICPositive)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Sample <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Sample, everything())

df_Hpos <- meta %>%
  left_join(abund_t, by = "Sample")


## transformation ####
df_Hpos_log <- df_Hpos

# Log-transform metabolite abundance data (log base 10)
df_Hpos_log[, 6:ncol(df_Hpos_log)] <- log10(df_Hpos_log[, 6:ncol(df_Hpos_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_Hpos_log)


## mean centering by SAMPLE ####

#set up to mean center logged data
df_Hpos_norm <-df_Hpos_log

## mean centering by sample

df_Hpos_norm[, 6:ncol(df_Hpos_norm)] <- t(scale(
  t(df_Hpos_norm[, 6:ncol(df_Hpos_norm)]),
  center = TRUE,
  scale = FALSE
))


# Check the transformed dataframe
head(df_Hpos_norm)

## PCA ####

Run_PCA_Hpos <- function() {
  
  # split normalized data by strain
  df_Hpos_ORWT <- df_Hpos_norm[df_Hpos_norm$Strain == "ORWT", ]
  df_Hpos_CRB  <- df_Hpos_norm[df_Hpos_norm$Strain == "CRB", ]
  
  # keep only complete cases for metabolite matrix
  keep_ORWT <- complete.cases(df_Hpos_ORWT[, 6:ncol(df_Hpos_ORWT)])
  keep_CRB  <- complete.cases(df_Hpos_CRB[, 6:ncol(df_Hpos_CRB)])
  
  df_Hpos_ORWT <- df_Hpos_ORWT[keep_ORWT, ]
  df_Hpos_CRB  <- df_Hpos_CRB[keep_CRB, ]
  
  # PCA input matrices
  pca_input_ORWT <- as.matrix(df_Hpos_ORWT[, 6:ncol(df_Hpos_ORWT)])
  pca_input_CRB  <- as.matrix(df_Hpos_CRB[, 6:ncol(df_Hpos_CRB)])
  
  # run PCA separately for each strain
  pca_Hpos_ORWT <- prcomp(pca_input_ORWT, scale = FALSE)
  pca_Hpos_CRB  <- prcomp(pca_input_CRB, scale = FALSE)
  
  # percent variance explained
  pca.var.ORWT <- pca_Hpos_ORWT$sdev^2
  pca.var.per.ORWT <- round(pca.var.ORWT / sum(pca.var.ORWT) * 100, 1)
  
  pca.var.CRB <- pca_Hpos_CRB$sdev^2
  pca.var.per.CRB <- round(pca.var.CRB / sum(pca.var.CRB) * 100, 1)
  
  # score dataframes
  pca.data.ORWT <- data.frame(
    Sample   = rownames(pca_Hpos_ORWT$x),
    X        = pca_Hpos_ORWT$x[, 1],
    Y        = pca_Hpos_ORWT$x[, 2],
    Age      = df_Hpos_ORWT$Age,
    PoolSize = df_Hpos_ORWT$PoolSize,
    Strain   = df_Hpos_ORWT$Strain
  )
  
  pca.data.CRB <- data.frame(
    Sample   = rownames(pca_Hpos_CRB$x),
    X        = pca_Hpos_CRB$x[, 1],
    Y        = pca_Hpos_CRB$x[, 2],
    Age      = df_Hpos_CRB$Age,
    PoolSize = df_Hpos_CRB$PoolSize,
    Strain   = df_Hpos_CRB$Strain
  )
  
  # plot ORWT
  p_Hpos_ORWT <- ggplot(
    pca.data.ORWT,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    xlab(paste("PC1 - ", pca.var.per.ORWT[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.ORWT[2], "%", sep = "")) +
    ggtitle("HILIC Positive - ORWT")
  
  # plot CRB
  p_Hpos_CRB <- ggplot(
    pca.data.CRB,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    xlab(paste("PC1 - ", pca.var.per.CRB[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.CRB[2], "%", sep = "")) +
    ggtitle("HILIC Positive - CRB")
  
  list(ORWT = p_Hpos_ORWT, CRB = p_Hpos_CRB)
}

Hpos_plots <- Run_PCA_Hpos()

Hpos_plots$ORWT
Hpos_plots$CRB

#### HILIC Negative

# change format for normalization

library(dplyr)

met_names <- df_HILICNegative$Metabolite
abund_mat <- df_HILICNegative[, 7:ncol(df_HILICNegative)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Sample <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Sample, everything())

df_Hneg <- meta %>%
  left_join(abund_t, by = "Sample")


## transformation ####
df_Hneg_log <- df_Hneg

# Log-transform metabolite abundance data (log base 10)
df_Hneg_log[, 6:ncol(df_Hneg_log)] <- log10(df_Hneg_log[, 6:ncol(df_Hneg_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_Hneg_log)


## mean centering by SAMPLE ####

#set up to mean center logged data
df_Hneg_norm <-df_Hneg_log

## mean centering by sample

df_Hneg_norm[, 6:ncol(df_Hneg_norm)] <- t(scale(
  t(df_Hneg_norm[, 6:ncol(df_Hneg_norm)]),
  center = TRUE,
  scale = FALSE
))


# Check the transformed dataframe
head(df_Hneg_norm)


## PCA ####

Run_PCA_Hneg <- function() {
  
  # split normalized data by strain
  df_Hneg_ORWT <- df_Hneg_norm[df_Hneg_norm$Strain == "ORWT", ]
  df_Hneg_CRB  <- df_Hneg_norm[df_Hneg_norm$Strain == "CRB", ]
  
  # keep only complete cases for metabolite matrix
  keep_ORWT <- complete.cases(df_Hneg_ORWT[, 6:ncol(df_Hneg_ORWT)])
  keep_CRB  <- complete.cases(df_Hneg_CRB[, 6:ncol(df_Hneg_CRB)])
  
  df_Hneg_ORWT <- df_Hneg_ORWT[keep_ORWT, ]
  df_Hneg_CRB  <- df_Hneg_CRB[keep_CRB, ]
  
  # PCA input matrices
  pca_input_ORWT <- as.matrix(df_Hneg_ORWT[, 6:ncol(df_Hneg_ORWT)])
  pca_input_CRB  <- as.matrix(df_Hneg_CRB[, 6:ncol(df_Hneg_CRB)])
  
  # run PCA separately for each strain
  pca_Hneg_ORWT <- prcomp(pca_input_ORWT, scale = FALSE)
  pca_Hneg_CRB  <- prcomp(pca_input_CRB, scale = FALSE)
  
  # percent variance explained
  pca.var.ORWT <- pca_Hneg_ORWT$sdev^2
  pca.var.per.ORWT <- round(pca.var.ORWT / sum(pca.var.ORWT) * 100, 1)
  
  pca.var.CRB <- pca_Hneg_CRB$sdev^2
  pca.var.per.CRB <- round(pca.var.CRB / sum(pca.var.CRB) * 100, 1)
  
  # score dataframes
  pca.data.ORWT <- data.frame(
    Sample   = rownames(pca_Hneg_ORWT$x),
    X        = pca_Hneg_ORWT$x[, 1],
    Y        = pca_Hneg_ORWT$x[, 2],
    Age      = df_Hneg_ORWT$Age,
    PoolSize = df_Hneg_ORWT$PoolSize,
    Strain   = df_Hneg_ORWT$Strain
  )
  
  pca.data.CRB <- data.frame(
    Sample   = rownames(pca_Hneg_CRB$x),
    X        = pca_Hneg_CRB$x[, 1],
    Y        = pca_Hneg_CRB$x[, 2],
    Age      = df_Hneg_CRB$Age,
    PoolSize = df_Hneg_CRB$PoolSize,
    Strain   = df_Hneg_CRB$Strain
  )
  
  # plot ORWT
  p_Hneg_ORWT <- ggplot(
    pca.data.ORWT,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    xlab(paste("PC1 - ", pca.var.per.ORWT[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.ORWT[2], "%", sep = "")) +
    ggtitle("HILIC Positive - ORWT")
  
  # plot CRB
  p_Hneg_CRB <- ggplot(
    pca.data.CRB,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    xlab(paste("PC1 - ", pca.var.per.CRB[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.CRB[2], "%", sep = "")) +
    ggtitle("HILIC Positive - CRB")
  
  list(ORWT = p_Hneg_ORWT, CRB = p_Hneg_CRB)
}

Hneg_plots <- Run_PCA_Hneg()

Hneg_plots <- Run_PCA_Hneg()

Hneg_plots$ORWT
Hneg_plots$CRB


#### C18 Negative

# change format for normalization

library(dplyr)

met_names <- df_C18Negative$Metabolite
abund_mat <- df_C18Negative[, 7:ncol(df_C18Negative)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Sample <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Sample, everything())

df_C18neg <- meta %>%
  left_join(abund_t, by = "Sample")


## transformation ####
df_C18neg_log <- df_C18neg

# Log-transform metabolite abundance data (log base 10)
df_C18neg_log[, 6:ncol(df_C18neg_log)] <- log10(df_C18neg_log[, 6:ncol(df_C18neg_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_C18neg_log)


## mean centering by SAMPLE ####

#set up to mean center logged data
df_C18neg_norm <-df_C18neg_log

## mean centering by sample

df_C18neg_norm[, 6:ncol(df_C18neg_norm)] <- t(scale(
  t(df_C18neg_norm[, 6:ncol(df_C18neg_norm)]),
  center = TRUE,
  scale = FALSE
))


# Check the transformed dataframe
head(df_C18neg_norm)


## PCA ####

Run_PCA_C18neg <- function() {
  
  # split normalized data by strain
  df_C18neg_ORWT <- df_C18neg_norm[df_C18neg_norm$Strain == "ORWT", ]
  df_C18neg_CRB  <- df_C18neg_norm[df_C18neg_norm$Strain == "CRB", ]
  
  # keep only complete cases for metabolite matrix
  keep_ORWT <- complete.cases(df_C18neg_ORWT[, 6:ncol(df_C18neg_ORWT)])
  keep_CRB  <- complete.cases(df_C18neg_CRB[, 6:ncol(df_C18neg_CRB)])
  
  df_C18neg_ORWT <- df_C18neg_ORWT[keep_ORWT, ]
  df_C18neg_CRB  <- df_C18neg_CRB[keep_CRB, ]
  
  # PCA input matrices
  pca_input_ORWT <- as.matrix(df_C18neg_ORWT[, 6:ncol(df_C18neg_ORWT)])
  pca_input_CRB  <- as.matrix(df_C18neg_CRB[, 6:ncol(df_C18neg_CRB)])
  
  # run PCA separately for each strain
  pca_C18neg_ORWT <- prcomp(pca_input_ORWT, scale = FALSE)
  pca_C18neg_CRB  <- prcomp(pca_input_CRB, scale = FALSE)
  
  # percent variance explained
  pca.var.ORWT <- pca_C18neg_ORWT$sdev^2
  pca.var.per.ORWT <- round(pca.var.ORWT / sum(pca.var.ORWT) * 100, 1)
  
  pca.var.CRB <- pca_C18neg_CRB$sdev^2
  pca.var.per.CRB <- round(pca.var.CRB / sum(pca.var.CRB) * 100, 1)
  
  # score dataframes
  pca.data.ORWT <- data.frame(
    Sample   = rownames(pca_C18neg_ORWT$x),
    X        = pca_C18neg_ORWT$x[, 1],
    Y        = pca_C18neg_ORWT$x[, 2],
    Age      = df_C18neg_ORWT$Age,
    PoolSize = df_C18neg_ORWT$PoolSize,
    Strain   = df_C18neg_ORWT$Strain
  )
  
  pca.data.CRB <- data.frame(
    Sample   = rownames(pca_C18neg_CRB$x),
    X        = pca_C18neg_CRB$x[, 1],
    Y        = pca_C18neg_CRB$x[, 2],
    Age      = df_C18neg_CRB$Age,
    PoolSize = df_C18neg_CRB$PoolSize,
    Strain   = df_C18neg_CRB$Strain
  )
  
  # plot ORWT
  p_C18neg_ORWT <- ggplot(
    pca.data.ORWT,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    xlab(paste("PC1 - ", pca.var.per.ORWT[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.ORWT[2], "%", sep = "")) +
    ggtitle("HILIC Positive - ORWT")
  
  # plot CRB
  p_C18neg_CRB <- ggplot(
    pca.data.CRB,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    xlab(paste("PC1 - ", pca.var.per.CRB[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.CRB[2], "%", sep = "")) +
    ggtitle("HILIC Positive - CRB")
  
  list(ORWT = p_C18neg_ORWT, CRB = p_C18neg_CRB)
}

C18neg_plots <- Run_PCA_C18neg()

C18neg_plots$ORWT
C18neg_plots$CRB


#### C18 Positive

# change format for normalization

library(dplyr)

met_names <- df_C18Positive$Metabolite
abund_mat <- df_C18Positive[, 7:ncol(df_C18Positive)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Sample <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Sample, everything())

df_C18pos <- meta %>%
  left_join(abund_t, by = "Sample")


## transformation ####
df_C18pos_log <- df_C18pos

# Log-transform metabolite abundance data (log base 10)
df_C18pos_log[, 6:ncol(df_C18pos_log)] <- log10(df_C18pos_log[, 6:ncol(df_C18pos_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_C18pos_log)


## mean centering by SAMPLE ####

#set up to mean center logged data
df_C18pos_norm <-df_C18pos_log

## mean centering by sample

df_C18pos_norm[, 6:ncol(df_C18pos_norm)] <- t(scale(
  t(df_C18pos_norm[, 6:ncol(df_C18pos_norm)]),
  center = TRUE,
  scale = FALSE
))


# Check the transformed dataframe
head(df_C18pos_norm)


## PCA ####

Run_PCA_C18pos <- function() {
  
  # split normalized data by strain
  df_C18pos_ORWT <- df_C18pos_norm[df_C18pos_norm$Strain == "ORWT", ]
  df_C18pos_CRB  <- df_C18pos_norm[df_C18pos_norm$Strain == "CRB", ]
  
  # keep only complete cases for metabolite matrix
  keep_ORWT <- complete.cases(df_C18pos_ORWT[, 6:ncol(df_C18pos_ORWT)])
  keep_CRB  <- complete.cases(df_C18pos_CRB[, 6:ncol(df_C18pos_CRB)])
  
  df_C18pos_ORWT <- df_C18pos_ORWT[keep_ORWT, ]
  df_C18pos_CRB  <- df_C18pos_CRB[keep_CRB, ]
  
  # PCA input matrices
  pca_input_ORWT <- as.matrix(df_C18pos_ORWT[, 6:ncol(df_C18pos_ORWT)])
  pca_input_CRB  <- as.matrix(df_C18pos_CRB[, 6:ncol(df_C18pos_CRB)])
  
  # run PCA separately for each strain
  pca_C18pos_ORWT <- prcomp(pca_input_ORWT, scale = FALSE)
  pca_C18pos_CRB  <- prcomp(pca_input_CRB, scale = FALSE)
  
  # percent variance explained
  pca.var.ORWT <- pca_C18pos_ORWT$sdev^2
  pca.var.per.ORWT <- round(pca.var.ORWT / sum(pca.var.ORWT) * 100, 1)
  
  pca.var.CRB <- pca_C18pos_CRB$sdev^2
  pca.var.per.CRB <- round(pca.var.CRB / sum(pca.var.CRB) * 100, 1)
  
  # score dataframes
  pca.data.ORWT <- data.frame(
    Sample   = rownames(pca_C18pos_ORWT$x),
    X        = pca_C18pos_ORWT$x[, 1],
    Y        = pca_C18pos_ORWT$x[, 2],
    Age      = df_C18pos_ORWT$Age,
    PoolSize = df_C18pos_ORWT$PoolSize,
    Strain   = df_C18pos_ORWT$Strain
  )
  
  pca.data.CRB <- data.frame(
    Sample   = rownames(pca_C18pos_CRB$x),
    X        = pca_C18pos_CRB$x[, 1],
    Y        = pca_C18pos_CRB$x[, 2],
    Age      = df_C18pos_CRB$Age,
    PoolSize = df_C18pos_CRB$PoolSize,
    Strain   = df_C18pos_CRB$Strain
  )
  
  # plot ORWT
  p_C18pos_ORWT <- ggplot(
    pca.data.ORWT,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    xlab(paste("PC1 - ", pca.var.per.ORWT[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.ORWT[2], "%", sep = "")) +
    ggtitle("HILIC Positive - ORWT")
  
  # plot CRB
  p_C18pos_CRB <- ggplot(
    pca.data.CRB,
    aes(x = X, y = Y, color = PoolSize, shape = Age)
  ) +
    guides(
      shape = guide_legend(title = "Age"),
      color = guide_legend(title = "Pool Size")
    ) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      legend.key = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    geom_point(size = 5) +
    xlab(paste("PC1 - ", pca.var.per.CRB[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per.CRB[2], "%", sep = "")) +
    ggtitle("HILIC Positive - CRB")
  
  list(ORWT = p_C18pos_ORWT, CRB = p_C18pos_CRB)
}

C18pos_plots <- Run_PCA_C18pos()

C18pos_plots$ORWT
C18pos_plots$CRB



#create plot object
Hneg_ORWT_plot <- Hneg_plots$ORWT + theme(plot.title = element_blank())
Hneg_CRB_plot  <- Hneg_plots$CRB  + theme(plot.title = element_blank())
Hpos_ORWT_plot <- Hpos_plots$ORWT + theme(plot.title = element_blank())
Hpos_CRB_plot  <- Hpos_plots$CRB  + theme(plot.title = element_blank())
C18neg_ORWT_plot <- C18neg_plots$ORWT + theme(plot.title = element_blank())
C18neg_CRB_plot  <- C18neg_plots$CRB  + theme(plot.title = element_blank())
C18pos_ORWT_plot <- C18pos_plots$ORWT + theme(plot.title = element_blank())
C18pos_CRB_plot  <- C18pos_plots$CRB  + theme(plot.title = element_blank())


#### Multi-panel plot ####

library(ggpubr)

multi_plot <- ggarrange(Hneg_ORWT_plot, Hneg_CRB_plot,
                        Hpos_ORWT_plot, Hpos_CRB_plot, 
                        C18neg_ORWT_plot, C18neg_CRB_plot,
                        C18pos_ORWT_plot, C18pos_CRB_plot,
          labels = c("A", "B", "C", "D","E","F", "G"),
          ncol = 2, nrow = 4,
          common.legend = TRUE,  # Ensures a shared legend
          legend = "right")     # Adjust legend position (options: "top", "bottom", "left", "right")

multi_plot


final_plot <- annotate_figure(
  multi_plot,
  bottom = text_grob(
    "A: ORWT (HILIC-)   B: CRB (HILIC-)   C: ORWT (HILIC+)   D: CRB (HILIC+)  E: ORWT (C18-)   F: CRB (C18-)  G: ORWT (C18+)  H: CRB (C18+)",
    size = 10, face = "bold"
  )
)

final_plot <- final_plot +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

final_plot



#### write normalized data to file ####
# create directory if it doesn't exist
if (!dir.exists("Mean_centered_by_Sample")) {
  dir.create("Mean_centered_by_Sample")
}


# save PCA

ggsave(
  filename = "Mean_centered_by_Sample/PCA_multi_panel_mean_centered_by_sample.tif",
  plot = final_plot,
  device = "tiff",
  width = 10,
  height = 12,
  dpi = 600
)

write.csv(df_C18pos_norm, file = "Mean_centered_by_Sample/C18_pos_normalized_sample.csv", row.names = FALSE)
write.csv(df_C18neg_norm, file = "Mean_centered_by_Sample/C18_neg_normalized_sample.csv", row.names = FALSE)
write.csv(df_Hpos_norm, file = "Mean_centered_by_Sample/Hilic_pos_normalized_sample.csv", row.names = FALSE)
write.csv(df_Hneg_norm, file = "Mean_centered_by_Sample/Hilic_neg_normalized_sample.csv", row.names = FALSE)


