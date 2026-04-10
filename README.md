# Metabolomics Sampling Project in *Drosophila melanogaster*

This repository contains all code and supplementary tables used to evaluate how pool size and biological replication influence metabolomic signal detection and statistical inference.

## Overview

Metabolomics studies often require pooling multiple individuals and balancing biological replication due to cost and sample limitations. This project investigates how these design choices affect the detection, stability, and interpretation of metabolite signals.

Using controlled experiments in *Drosophila melanogaster*, we systematically vary pool size (5, 50, 100 individuals) and biological replication, and apply statistical modeling and downsampling approaches to quantify their effects on metabolomic inference.

## Contents

### Scripts
All R scripts used for:
- Data preprocessing and quality control (CV filtering, normalization)
- Principal component analysis (PCA)
- Linear modeling of metabolite abundance
- Replicate downsampling and signal detection analysis
- Figure generation for main and supplementary results
- Functional module assignment and biological interpretation

### Supplementary Tables
- Metabolite annotations and functional module assignments
- Summary statistics from linear models and downsampling analyses
- Additional outputs used in figures and interpretation

## Data Processing and Analysis

The analysis workflow includes:
1. Filtering metabolites based on technical variability (coefficient of variation > 0.30 removed)
2. Log10 transformation of metabolite abundances
3. Normalization by mean-centering (by metabolite or sample, depending on analysis)
4. Principal component analysis across LC–MS panels
5. Linear modeling to identify diet-associated metabolites
6. Downsampling of replicate populations to evaluate signal retention
7. Quantification of true/false positives, effect size stability, and detection frequency
8. Functional grouping of metabolites into curated biochemical modules

## Reproducibility

All scripts are organized to reproduce the analyses and figures presented in the manuscript. Intermediate files are generated within scripts where necessary.

## Requirements

Analyses were performed in R (version 4.x). Required packages include:

- dplyr  
- tidyr  
- ggplot2  
- lme4  
- lmerTest  
- emmeans  
- ggpubr / patchwork  
- future.apply  

## Authors

David L. Hubert & Mark A. Phillips 
Oregon State University  

## Notes

This repository is intended to accompany the manuscript and provide full transparency and reproducibility of all analyses. Scripts are annotated and structured to reflect the analytical workflow described in the Methods section.
