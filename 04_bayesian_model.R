#!/usr/bin/env Rscript
################################################################################
# 02_isotope_analysis.R
# Figure 1: Phosphate Oxygen Isotope Deviations from Equilibrium
#
# Reproduces Figure 1 from Arfania et al. (2026) ES&T
# Shows Δδ¹⁸O_P across landscape positions and tillage systems,
# colored by Hedley P fraction with monthly sampling shapes.
################################################################################

cat("=== Figure 1: Isotope Deviation Analysis ===\n")

# ---------------------------------------------------------------------------
# Summary statistics for mean points and error bars
# ---------------------------------------------------------------------------

# Create a simplified landscape-fraction grouping for plotting
df$LandPos <- interaction(
  gsub("slope", "", as.character(df$Position)),
  gsub("Conv. Tillage", "CT", gsub("No-Till", "NT", as.character(df$Treatment))),
  sep = "-"
)
df$LandPos <- factor(df$LandPos, levels = c("Top-CT", "Top-NT", "Toe-CT", "Toe-NT"))

summary_stats <- df %>%
  group_by(LandPos, Fraction) %>%
  summarise(
    mean_dev = mean(d18Op_deviation),
    se_dev   = sd(d18Op_deviation) / sqrt(n()),
    .groups  = "drop"
  )

# ---------------------------------------------------------------------------
# Figure 1: All fractions scatter plot
# ---------------------------------------------------------------------------

fig1 <- ggplot(df, aes(x = LandPos, y = d18Op_deviation)) +
  
  # Equilibrium zone (±2‰)
  annotate("rect", xmin = 0.4, xmax = 4.6, ymin = -2, ymax = 2,
           alpha = 0.1, fill = "lightblue") +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "grey50", alpha = 0.5) +
  
  # Individual observations
  geom_point(aes(color = Fraction, shape = Date),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
             size = 3, alpha = 0.7) +
  
  # Mean ± SE error bars
  geom_errorbar(data = summary_stats,
                mapping = aes(x = LandPos,
                              ymin = mean_dev - se_dev,
                              ymax = mean_dev + se_dev,
                              color = Fraction),
                width = 0.12,
                position = position_dodge(width = 0.65),
                inherit.aes = FALSE) +
  
  # Mean diamonds
  geom_point(data = summary_stats,
             mapping = aes(x = LandPos, y = mean_dev, color = Fraction),
             shape = 18, size = 4,
             position = position_dodge(width = 0.65),
             inherit.aes = FALSE) +
  
  # Scales
  scale_color_manual(
    values = colors_fraction,
    name   = "Hedley P Fraction",
    labels = c(
      expression(H[2]*O ~ "-" ~ P[i]),
      expression(NaHCO[3] ~ "-" ~ P[i]),
      expression(NaOH ~ "-" ~ P[i]),
      expression(HNO[3] ~ "-" ~ P[i])
    )
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Sampling Month") +
  scale_y_continuous(breaks = seq(-4, 8, 2), limits = c(-4, 8)) +
  scale_x_discrete(labels = c("Top-CT" = "Top\nCT",
                              "Top-NT" = "Top\nNT",
                              "Toe-CT" = "Toe\nCT",
                              "Toe-NT" = "Toe\nNT")) +
  
  # Labels
  labs(
    x = "Landscape Position – Tillage System",
    y = expression(Delta*delta^{18}*O[P] ~ "(‰ from equilibrium)")
  ) +
  
  # Equilibrium zone annotation
  annotate("text", x = 2.5, y = 1.5,
           label = "Equilibrium Zone (±2‰)",
           size = 3.5, fontface = "italic", color = "steelblue") +
  
  # Theme
  theme_est() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(fig1)
save_figure(fig1, "Figure1_Isotope_Deviations", width = 10, height = 7)

cat("  Figure 1 complete\n\n")
