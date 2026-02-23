#!/usr/bin/env Rscript
################################################################################
# run_all.R
# Master Script to Reproduce All Analyses and Figures
#
# Arfania, H., Kayler, Z., Strawn, D., Brooks, E., & Laan, M. (2026)
# "Functional Decoupling of Active and Legacy Phosphorus Pools
#  Revealed by Phosphate Oxygen Isotopes (δ¹⁸O_P)"
# Environmental Science & Technology
#
# Usage:
#   source("R/run_all.R")
#   # or from terminal: Rscript R/run_all.R
################################################################################

cat("\n")
cat("================================================================\n")
cat("  Arfania et al. (2026) ES&T - Full Analysis Pipeline\n")
cat("  Legacy Phosphorus Isotope Analysis\n")
cat("================================================================\n\n")

start_time <- Sys.time()

# Run scripts in order
source("R/00_setup.R")
source("R/01_data_preparation.R")
source("R/02_isotope_analysis.R")
source("R/03_pca_analysis.R")
source("R/04_bayesian_model.R")
source("R/05_conceptual_figure.R")

# Summary
elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 1)

cat("================================================================\n")
cat("  Analysis Complete\n")
cat("================================================================\n\n")

cat("Figures generated:\n")
cat("  Figure 1: Isotope deviations across fractions and landscapes\n")
cat("  Figure 2: PCA biplots (A: by fraction; B: by position/tillage)\n")
cat("  Figure 3: Bayesian transformation model (4 panels)\n")
cat("  Figure 4: Two-domain conceptual model\n\n")

cat("Output directory: output/\n")
cat(sprintf("  Figures: %d files\n", length(list.files("output/figures"))))
cat(sprintf("  Tables:  %d files\n", length(list.files("output/tables"))))
cat(sprintf("  PCA:     %d files\n", length(list.files("output/PCA_results"))))
cat(sprintf("\nTotal time: %s seconds\n", elapsed))
cat("================================================================\n")
