# ============================================================
# PCA Analysis of Metabolomic Data (Mean-Centered by Metabolite)
#
# This script performs preprocessing and principal component 
# analysis (PCA) on metabolomics data across four LC-MS panels 
# (C18 positive, C18 negative, HILIC positive, HILIC negative),
# using mean-centering by metabolite.
#
# Workflow:
# 1. Import raw metabolite abundance data and sample metadata.
# 2. Perform quality control filtering by removing metabolites 
#    with high technical variability (coefficient of variation > 0.30 
#    in pooled QC samples).
# 3. Log10-transform metabolite abundances (with pseudocount) to 
#    stabilize variance.
# 4. Normalize data by mean-centering each metabolite across samples 
#    to remove metabolite-specific intensity differences while 
#    preserving variation among samples.
# 5. Perform PCA on normalized metabolite data for each LC-MS panel.
# 6. Generate PCA plots with samples colored by pool size and 
#    shaped by diet.
# 7. Combine individual PCA plots into a multi-panel figure.
# 8. Export normalized datasets and high-resolution PCA figures.
#
# This analysis emphasizes between-sample variation in metabolomic 
# profiles and is used to assess how pool size influences overall 
# metabolomic structure.
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

# add X to tube name to match df
meta$Tube <- paste0("X", meta$Tube)
meta$PoolSize <- factor(meta$PoolSize)
meta$Population <- factor(meta$Population)



# check QC for issues

qc_cols <- c("QC.01","QC.02","QC.03")

# calculate CV
df$QC_CV <- apply(df[, qc_cols], 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

# inspect distribution

summary(df$QC_CV)
hist(df$QC_CV, breaks = 50)

# filter for CV > .30
df_filtered <- df[df$QC_CV < 0.30, ] # 9 metabolites removed due to high CV in QC samples

# visualization of this
ggplot(df, aes(x = QC_CV)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0.3, color = "red") +
  theme_classic()

#remove QC columns from df_filtered

df_filtered <- df_filtered[, 1:(ncol(df_filtered) - 4)]

#create separate DF by run type

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
abund_t$Tube <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Tube, everything())

df_Hpos <- meta %>%
  left_join(abund_t, by = "Tube")


## transformation ####
df_Hpos_log <- df_Hpos

# Log-transform metabolite abundance data (log base 10)
df_Hpos_log[, 9:ncol(df_Hpos_log)] <- log10(df_Hpos_log[, 9:ncol(df_Hpos_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_Hpos_log)


## mean centering by METABOLITE ####

#set up to mean center logged data
df_Hpos_norm <-df_Hpos_log

## mean centering by metabolite

df_Hpos_norm[, 9:ncol(df_Hpos_norm)] <- scale(
  df_Hpos_norm[, 9:ncol(df_Hpos_norm)],
  center = TRUE,
  scale = FALSE
)

# Check the transformed dataframe
head(df_Hpos_norm)


## PCA ####

pca_df_Hpos <-as.matrix(na.omit(df_Hpos_norm[,9:ncol(df_Hpos_norm)]))
pca_df_Hpos <- prcomp(pca_df_Hpos, scale=FALSE) 

pca.var_df_Hpos <- pca_df_Hpos$sdev^2
pca.var.per_df_Hpos <- round(pca.var_df_Hpos/sum(pca.var_df_Hpos)*100, 1)

# diet = shape, pool size = color

Run_PCA_Hpos <- function() {
  pca.data_df_Hpos <- data.frame(Sample = rownames(pca_df_Hpos$x),
                                 X = pca_df_Hpos$x[,1],
                                 Y = pca_df_Hpos$x[,2])
  
  p_df_Hpos <- ggplot(pca.data_df_Hpos, aes(x = X, y = Y, label = df_Hpos_norm$Diet, 
                                            color = df_Hpos$PoolSize, shape = df_Hpos$Diet)) +
    guides(shape = guide_legend(title = "Diet"), 
           color = guide_legend(title = "Pool Size")) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          legend.key = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(color = "black")) +
    geom_point(size=5) + 
    xlab(paste("PC1 - ", pca.var.per_df_Hpos[1], "%", sep = "")) + 
    ylab(paste("PC2 - ", pca.var.per_df_Hpos[2], "%", sep = "")) + 
    theme(legend.key.size = unit(1, 'cm'))
  
  p_df_Hpos
}

Run_PCA_Hpos()


#### HILIC Negative

# change format for normalization

library(dplyr)

met_names <- df_HILICNegative$Metabolite
abund_mat <- df_HILICNegative[, 7:ncol(df_HILICNegative)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Tube <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Tube, everything())

df_Hneg <- meta %>%
  left_join(abund_t, by = "Tube")


## transformation ####
df_Hneg_log <- df_Hneg

# Log-transform metabolite abundance data (log base 10)
df_Hneg_log[, 9:ncol(df_Hneg_log)] <- log10(df_Hneg_log[, 9:ncol(df_Hneg_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_Hneg_log)


## mean centering ####

#set up to mean center logged data
df_Hneg_norm <-df_Hneg_log

## mean centering by metabolite

df_Hneg_norm[, 9:ncol(df_Hneg_norm)] <- scale(
  df_Hneg_norm[, 9:ncol(df_Hneg_norm)],
  center = TRUE,
  scale = FALSE
)


# Check the transformed dataframe
head(df_Hneg_norm)


## PCA ####

pca_df_Hneg <-as.matrix(na.omit(df_Hneg_norm[,9:ncol(df_Hneg_norm)]))
pca_df_Hneg <- prcomp(pca_df_Hneg, scale=FALSE) 

pca.var_df_Hneg <- pca_df_Hneg$sdev^2
pca.var.per_df_Hneg <- round(pca.var_df_Hneg/sum(pca.var_df_Hneg)*100, 1)

# diet = shape, pool size = color

Run_PCA_Hneg <- function() {
  pca.data_df_Hneg <- data.frame(Sample = rownames(pca_df_Hneg$x),
                                 X = pca_df_Hneg$x[,1],
                                 Y = pca_df_Hneg$x[,2])
  
  p_df_Hneg <- ggplot(pca.data_df_Hneg, aes(x = X, y = Y, label = df_Hneg_norm$Diet, 
                                            color = df_Hneg$PoolSize, shape = df_Hneg$Diet)) +
    guides(shape = guide_legend(title = "Diet"), 
           color = guide_legend(title = "Pool Size")) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          legend.key = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(color = "black")) +
    geom_point(size=5) + 
    xlab(paste("PC1 - ", pca.var.per_df_Hneg[1], "%", sep = "")) + 
    ylab(paste("PC2 - ", pca.var.per_df_Hneg[2], "%", sep = "")) + 
    theme(legend.key.size = unit(1, 'cm'))
  
  p_df_Hneg
}

Run_PCA_Hneg()



#### C18 Negative

# change format for normalization

library(dplyr)

met_names <- df_C18Negative$Metabolite
abund_mat <- df_C18Negative[, 7:ncol(df_C18Negative)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Tube <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Tube, everything())

df_C18neg <- meta %>%
  left_join(abund_t, by = "Tube")


## transformation ####
df_C18neg_log <- df_C18neg

# Log-transform metabolite abundance data (log base 10)
df_C18neg_log[, 9:ncol(df_C18neg_log)] <- log10(df_C18neg_log[, 9:ncol(df_C18neg_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_C18neg_log)


## mean centering ####

#set up to mean center logged data
df_C18neg_norm <-df_C18neg_log

## mean centering by metabolite

df_C18neg_norm[, 9:ncol(df_C18neg_norm)] <- scale(
  df_C18neg_norm[, 9:ncol(df_C18neg_norm)],
  center = TRUE,
  scale = FALSE
)


# Check the transformed dataframe
head(df_C18neg_norm)


## PCA ####

pca_df_C18neg <-as.matrix(na.omit(df_C18neg_norm[,9:ncol(df_C18neg_norm)]))
pca_df_C18neg <- prcomp(pca_df_C18neg, scale=FALSE) 

pca.var_df_C18neg <- pca_df_C18neg$sdev^2
pca.var.per_df_C18neg <- round(pca.var_df_C18neg/sum(pca.var_df_C18neg)*100, 1)

# diet = shape, pool size = color

Run_PCA_C18neg <- function() {
  pca.data_df_C18neg <- data.frame(Sample = rownames(pca_df_C18neg$x),
                                   X = pca_df_C18neg$x[,1],
                                   Y = pca_df_C18neg$x[,2])
  
  p_df_C18neg <- ggplot(pca.data_df_C18neg, aes(x = X, y = Y, label = df_C18neg_norm$Diet, 
                                                color = df_C18neg$PoolSize, shape = df_C18neg$Diet)) +
    guides(shape = guide_legend(title = "Diet"), 
           color = guide_legend(title = "Pool Size")) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          legend.key = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(color = "black")) +
    geom_point(size=5) + 
    xlab(paste("PC1 - ", pca.var.per_df_C18neg[1], "%", sep = "")) + 
    ylab(paste("PC2 - ", pca.var.per_df_C18neg[2], "%", sep = "")) + 
    theme(legend.key.size = unit(1, 'cm'))
  
  p_df_C18neg
}

Run_PCA_C18neg()


#### C18 Positive

# change format for normalization

library(dplyr)

met_names <- df_C18Positive$Metabolite
abund_mat <- df_C18Positive[, 7:ncol(df_C18Positive)]

abund_t <- as.data.frame(t(abund_mat), stringsAsFactors = FALSE)
colnames(abund_t) <- met_names
abund_t$Tube <- rownames(abund_t)
rownames(abund_t) <- NULL

abund_t <- abund_t %>%
  select(Tube, everything())

df_C18pos <- meta %>%
  left_join(abund_t, by = "Tube")


## transformation ####
df_C18pos_log <- df_C18pos

# Log-transform metabolite abundance data (log base 10)
df_C18pos_log[, 9:ncol(df_C18pos_log)] <- log10(df_C18pos_log[, 9:ncol(df_C18pos_log)] + 1)  # Adding 1 to avoid log(0) issues

# Check the transformed dataframe
head(df_C18pos_log)


## mean centering ####

#set up to mean center logged data
df_C18pos_norm <-df_C18pos_log

## mean centering by metabolite

df_C18pos_norm[, 9:ncol(df_C18pos_norm)] <- scale(
  df_C18pos_norm[, 9:ncol(df_C18pos_norm)],
  center = TRUE,
  scale = FALSE
)

# Check the transformed dataframe
head(df_C18pos_norm)


## PCA ####

pca_df_C18pos <-as.matrix(na.omit(df_C18pos_norm[,9:ncol(df_C18pos_norm)]))
pca_df_C18pos <- prcomp(pca_df_C18pos, scale=FALSE) 

pca.var_df_C18pos <- pca_df_C18pos$sdev^2
pca.var.per_df_C18pos <- round(pca.var_df_C18pos/sum(pca.var_df_C18pos)*100, 1)

# diet = shape, pool size = color

Run_PCA_C18pos <- function() {
  pca.data_df_C18pos <- data.frame(Sample = rownames(pca_df_C18pos$x),
                                   X = pca_df_C18pos$x[,1],
                                   Y = pca_df_C18pos$x[,2])
  
  p_df_C18pos <- ggplot(pca.data_df_C18pos, aes(x = X, y = Y, label = df_C18pos_norm$Diet, 
                                                color = df_C18pos$PoolSize, shape = df_C18pos$Diet)) +
    guides(shape = guide_legend(title = "Diet"), 
           color = guide_legend(title = "Pool Size")) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          legend.key = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(color = "black")) +
    geom_point(size=5) + 
    xlab(paste("PC1 - ", pca.var.per_df_C18pos[1], "%", sep = "")) + 
    ylab(paste("PC2 - ", pca.var.per_df_C18pos[2], "%", sep = "")) + 
    theme(legend.key.size = unit(1, 'cm'))
  
  p_df_C18pos
}

Run_PCA_C18pos()


#create plot object
C18pos_plot <- Run_PCA_C18pos()
C18neg_plot <- Run_PCA_C18neg()
Hpos_plot <- Run_PCA_Hpos()
Hneg_plot <- Run_PCA_Hneg()


#### Multi-panel plot ####

library(ggpubr)

multi_plot <- ggarrange(C18pos_plot, C18neg_plot, Hpos_plot, Hneg_plot,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE,  # Ensures a shared legend
          legend = "right")     # Adjust legend position (options: "top", "bottom", "left", "right")

multi_plot


final_plot <- annotate_figure(
  multi_plot,
  bottom = text_grob(
    "A: C18 Positive   B: C18 Negative   C: HILIC Positive   D: HILIC Negative",
    size = 10
  )
)

final_plot <- final_plot +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

final_plot




#### write normalized data to file ####
# create directory if it doesn't exist
if (!dir.exists("Mean_centered_by_Metabolite")) {
  dir.create("Mean_centered_by_Metabolite")
}


# save PCA

ggsave(
  filename = "Mean_centered_by_Metabolite/PCA_multi_panel_mean_centered_by_metabolite.tif",
  plot = final_plot,
  device = "tiff",
  width = 10,
  height = 8,
  dpi = 600
)

write.csv(df_C18pos_norm, file = "Mean_centered_by_Metabolite/C18_pos_normalized_metabolite.csv", row.names = FALSE)
write.csv(df_C18neg_norm, file = "Mean_centered_by_Metabolite/C18_neg_normalized_metabolite.csv", row.names = FALSE)
write.csv(df_Hpos_norm, file = "Mean_centered_by_Metabolite/Hilic_pos_normalized_metabolite.csv", row.names = FALSE)
write.csv(df_Hneg_norm, file = "Mean_centered_by_Metabolite/Hilic_neg_normalized_metabolite.csv", row.names = FALSE)


