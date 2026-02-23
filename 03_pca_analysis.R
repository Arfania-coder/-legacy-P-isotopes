#!/usr/bin/env Rscript
################################################################################
# 00_setup.R
# Package Installation and Configuration
#
# Part of: Arfania et al. (2026) ES&T
# "Functional Decoupling of Active and Legacy Phosphorus Pools
#  Revealed by Phosphate Oxygen Isotopes"
################################################################################

cat("\n=== Setting up environment ===\n")

# Required packages
packages <- c(
  # Core data manipulation and visualization
  "tidyverse", "ggplot2", "scales",
  # Multivariate analysis
  "FactoMineR", "factoextra",
  # Layout and annotation
  "gridExtra", "grid", "ggrepel", "patchwork",
  # Additional visualization
  "RColorBrewer", "viridis", "ggridges",
  # Statistical comparisons
  "ggpubr", "emmeans",
  # Correlation and heatmaps
  "ggcorrplot", "reshape2",
  # Network visualization
  "igraph", "ggraph", "tidygraph",
  # Dendrogram
  "ggdendro"
)

# Install missing packages
installed <- installed.packages()[, "Package"]
to_install <- packages[!packages %in% installed]

if (length(to_install) > 0) {
  cat(sprintf("  Installing %d packages: %s\n", 
              length(to_install), paste(to_install, collapse = ", ")))
  install.packages(to_install, repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Load all packages
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

cat(sprintf("  Loaded %d packages\n", length(packages)))

# Global options
options(scipen = 999, digits = 3)

# Output directory setup (relative to project root)
OUTPUT_DIR <- "output"
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "PCA_results"), showWarnings = FALSE, recursive = TRUE)

# Color schemes used throughout the paper
colors_fraction <- c(
  "DI"     = "#2166AC",
  "NaHCO3" = "#762A83",
  "NaOH"   = "#5AAE61",
  "HNO3"   = "#D6604D"
)

colors_position <- c(
  "Topslope"  = "#1B9E77",
  "Toeslope"  = "#D95F02"
)

colors_landscape <- c(
  "Top-CT" = "#8C510A", "Top-NT" = "#01665E",
  "Toe-CT" = "#DFC27D", "Toe-NT" = "#80CDC1"
)

colors_date <- c(
  "December" = "#377EB8",
  "January"  = "#4DAF4A",
  "February" = "#E41A1C"
)

# Publication theme for ggplot2
theme_est <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      axis.title    = element_text(size = 13, face = "bold", color = "black"),
      axis.text     = element_text(size = 11, color = "black"),
      legend.title  = element_text(size = 12, face = "bold", color = "black"),
      legend.text   = element_text(size = 10, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey90")
    )
}

# Helper function to save figures in multiple formats
save_figure <- function(plot, filename, width = 10, height = 7, dpi_png = 300, dpi_tiff = 600) {
  base <- file.path(OUTPUT_DIR, "figures", filename)
  
  ggsave(paste0(base, ".pdf"),  plot, width = width, height = height)
  ggsave(paste0(base, ".png"),  plot, width = width, height = height, dpi = dpi_png)
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = dpi_tiff,
         compression = "lzw")
  
  cat(sprintf("  Saved: %s (.pdf, .png, .tiff)\n", filename))
}

cat("  Setup complete\n\n")
