#!/usr/bin/env Rscript
################################################################################
# run_all.R
# Master Script to Reproduce All Analyses and Figures
#
# Arfania, H., Kayler, Z., Strawn, D., Brooks, E., & Laan, M. (2026)
# Environmental Science & Technology
#
# Usage: source("R/run_all.R")
################################################################################

rm(list = ls())
graphics.off()

cat("================================================================\n")
cat("  Arfania et al. (2026) ES&T - Full Analysis Pipeline\n")
cat("================================================================\n\n")

start_time <- Sys.time()

source("R/00_setup.R")
source("R/01_data_preparation.R")
source("R/02_isotope_analysis.R")
source("R/03_pca_analysis.R")
source("R/04_bayesian_model.R")
source("R/05_conceptual_figure.R")

elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 1)

cat("\n================================================================\n")
cat("  Analysis Complete\n")
cat(sprintf("  Total time: %s seconds\n", elapsed))
cat("================================================================\n")

sessionInfo()
