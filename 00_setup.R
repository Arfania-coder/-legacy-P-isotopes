#!/usr/bin/env Rscript
################################################################################
# 01_data_preparation.R
# Load and preprocess phosphate oxygen isotope data
#
# Part of: Arfania et al. (2026) ES&T
################################################################################

cat("=== Loading and preparing data ===\n")

# Load data
df <- read.csv("data/isotope_data.csv", stringsAsFactors = FALSE)

# Set factor levels
df$Date <- factor(df$Date, levels = c("Dec", "Jan", "Feb"),
                  labels = c("December", "January", "February"))
df$Position <- factor(df$Position, levels = c("Top", "Toe"),
                      labels = c("Topslope", "Toeslope"))
df$Treatment <- factor(df$Treatment, levels = c("CT", "NT"),
                       labels = c("Conv. Tillage", "No-Till"))
df$Fraction <- factor(df$Fraction, levels = c("DI", "NaHCO3", "NaOH", "HNO3"))

# Derived variables
df$Landscape <- interaction(df$Position, df$Treatment, sep = " - ")
df$Sample_ID <- paste(df$Position, df$Treatment, df$Fraction, df$Date, sep = "_")

cat(sprintf("  Loaded %d observations\n", nrow(df)))
cat(sprintf("  Fractions: %s\n", paste(levels(df$Fraction), collapse = ", ")))
cat(sprintf("  Treatments: %s\n", paste(levels(df$Treatment), collapse = ", ")))
cat(sprintf("  Positions: %s\n", paste(levels(df$Position), collapse = ", ")))

# Summary statistics by fraction (Table S1)
fraction_summary <- df %>%
  group_by(Fraction) %>%
  summarise(
    n             = n(),
    mean_d18Op    = round(mean(d18Op_soil), 2),
    sd_d18Op      = round(sd(d18Op_soil), 2),
    mean_dev      = round(mean(d18Op_deviation), 2),
    sd_dev        = round(sd(d18Op_deviation), 2),
    pct_in_equil  = round(100 * mean(abs(d18Op_deviation) <= 2), 1),
    .groups = "drop"
  )

cat("\n  Summary by fraction:\n")
print(fraction_summary)

write.csv(fraction_summary, 
          file.path(OUTPUT_DIR, "tables/Table_S1_Fraction_Summary.csv"),
          row.names = FALSE)

# Summary by fraction, position, and treatment
detailed_summary <- df %>%
  group_by(Fraction, Position, Treatment) %>%
  summarise(
    n        = n(),
    mean_dev = round(mean(d18Op_deviation), 2),
    se_dev   = round(sd(d18Op_deviation) / sqrt(n()), 2),
    .groups  = "drop"
  )

write.csv(detailed_summary, 
          file.path(OUTPUT_DIR, "tables/Table_S2_Detailed_Summary.csv"),
          row.names = FALSE)

cat("  Data preparation complete\n\n")
