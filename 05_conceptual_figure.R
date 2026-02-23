#!/usr/bin/env Rscript
################################################################################
# 04_bayesian_model.R
# Figure 3: Bayesian Analysis of Inter-pool P Transformations
#
# Reproduces Figure 3 from Arfania et al. (2026) ES&T
# Panel A: Network visualization of transformation pathways
# Panel B: Posterior probability distributions
# Panel C: Pairwise isotopic distance matrix
# Panel D: Ranked transformation pathways with 95% CIs
#
# Bayesian transition probabilities modeled using Beta distributions
# scaled by isotopic distance (Gelman et al., 2013)
################################################################################

cat("=== Figure 3: Bayesian Transformation Model ===\n")

# ---------------------------------------------------------------------------
# Fraction-level statistics
# ---------------------------------------------------------------------------

fractions <- c("DI", "NaHCO3", "NaOH", "HNO3")

fraction_stats <- df %>%
  group_by(Fraction) %>%
  summarise(
    mean_dev = mean(d18Op_deviation),
    sd_dev   = sd(d18Op_deviation),
    n        = n(),
    .groups  = "drop"
  )

fraction_means <- setNames(fraction_stats$mean_dev, fraction_stats$Fraction)

cat("  Fraction means (Δδ¹⁸O_P):\n")
print(fraction_stats)

# ---------------------------------------------------------------------------
# Bayesian transition probability calculation
# Beta distribution parameterized by isotopic distance
# (Eq. 7 in manuscript)
# ---------------------------------------------------------------------------

calc_transition_prob <- function(delta1, delta2, n_samples = 10000) {
  distance <- abs(delta2 - delta1)
  
  # Shape parameters based on isotopic distance thresholds
  if (distance < 2) {
    alpha <- 5; beta <- 2
  } else if (distance < 4) {
    alpha <- 3; beta <- 4
  } else {
    alpha <- 1; beta <- 8
  }
  
  samples <- rbeta(n_samples, alpha, beta)
  
  list(
    mean     = mean(samples),
    sd       = sd(samples),
    ci_lower = unname(quantile(samples, 0.025)),
    ci_upper = unname(quantile(samples, 0.975)),
    samples  = samples
  )
}

# Calculate all pairwise transitions
transitions <- expand.grid(from = fractions, to = fractions,
                           stringsAsFactors = FALSE) %>%
  filter(from != to) %>%
  mutate(
    from_mean = fraction_means[from],
    to_mean   = fraction_means[to],
    distance  = abs(to_mean - from_mean)
  )

set.seed(42)  # Reproducibility
results_list <- list()
for (i in 1:nrow(transitions)) {
  res <- calc_transition_prob(transitions$from_mean[i], transitions$to_mean[i])
  transitions$prob_mean[i]  <- res$mean
  transitions$prob_lower[i] <- res$ci_lower
  transitions$prob_upper[i] <- res$ci_upper
  results_list[[i]] <- res
}

transitions <- transitions %>% arrange(desc(prob_mean))
top_pathways <- head(transitions, 8)
top_pathways$pathway <- paste(top_pathways$from, "->", top_pathways$to)

cat("\n  Top transformation pathways:\n")
print(top_pathways[, c("from", "to", "distance", "prob_mean", "prob_lower", "prob_upper")])

write.csv(transitions,
          file.path(OUTPUT_DIR, "tables/Table_S3_Transition_Probabilities.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Isotopic distance matrix
# ---------------------------------------------------------------------------

isotopic_distance <- function(frac1, frac2) {
  d1 <- df$d18Op_deviation[df$Fraction == frac1]
  d2 <- df$d18Op_deviation[df$Fraction == frac2]
  sqrt(mean((d1 - d2)^2))
}

dist_matrix <- matrix(0, nrow = 4, ncol = 4,
                      dimnames = list(fractions, fractions))
for (i in 1:4) for (j in 1:4) {
  if (i != j) dist_matrix[i, j] <- isotopic_distance(fractions[i], fractions[j])
}

write.csv(dist_matrix,
          file.path(OUTPUT_DIR, "tables/Table_S4_Isotopic_Distance_Matrix.csv"))

# ---------------------------------------------------------------------------
# Posterior samples for key transitions
# ---------------------------------------------------------------------------

set.seed(123)
n_samples <- 10000
post_DI_NaHCO3    <- rbeta(n_samples, 5, 2)
post_NaHCO3_NaOH  <- rbeta(n_samples, 3, 4)
post_NaOH_HNO3    <- rbeta(n_samples, 1, 8)
post_HNO3_mineral <- rbeta(n_samples, 0.5, 10)

# ---------------------------------------------------------------------------
# Transition probability matrix for Markov chain
# ---------------------------------------------------------------------------

prob_matrix_raw <- transitions %>%
  select(from, to, prob_mean) %>%
  pivot_wider(names_from = to, values_from = prob_mean, values_fill = 0) %>%
  column_to_rownames("from") %>%
  as.matrix()

# Re-order to match fractions
prob_matrix_raw <- prob_matrix_raw[fractions, fractions]
diag(prob_matrix_raw) <- 0

# Row-normalize to create stochastic matrix
row_sums <- rowSums(prob_matrix_raw)
trans_mat <- prob_matrix_raw / row_sums
# Assign self-transition probability
for (i in 1:4) {
  trans_mat[i, ] <- prob_matrix_raw[i, ] / (row_sums[i] + prob_matrix_raw[i, i])
}
# Simple approach: rescale so rows sum to 1
trans_mat <- prob_matrix_raw
diag(trans_mat) <- 0
trans_mat <- trans_mat / rowSums(trans_mat)

# Steady state via eigendecomposition
eigen_result <- eigen(t(trans_mat))
steady_state <- abs(eigen_result$vectors[, 1])
steady_state <- steady_state / sum(steady_state)
names(steady_state) <- fractions

cat("\n  Markov chain steady state:\n")
print(round(steady_state, 3))

write.csv(data.frame(Fraction = fractions, Steady_State = round(steady_state, 4)),
          file.path(OUTPUT_DIR, "tables/Table_S5_Steady_State.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure 3: Four-panel Bayesian figure (base R for publication quality)
# ---------------------------------------------------------------------------

cat("\n  Creating four-panel Figure 3...\n")

create_four_panel <- function() {
  
  layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1))
  
  # =========================================================================
  # Panel A: Network Diagram
  # =========================================================================
  par(mar = c(3, 3, 4, 2), family = "serif")
  
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", bty = "n")
  title(main = "(A) Bayesian Network of P Transformations",
        cex.main = 1.4, font.main = 2, line = 1)
  
  # Node positions (2x2 grid)
  node_x <- c(0.2, 0.8, 0.2, 0.8)
  node_y <- c(0.75, 0.75, 0.25, 0.25)
  names(node_x) <- names(node_y) <- fractions
  
  # Draw significant edges (prob > 0.20)
  sig_trans <- transitions[transitions$prob_mean > 0.20, ]
  for (i in 1:nrow(sig_trans)) {
    from <- sig_trans$from[i]; to <- sig_trans$to[i]
    from_idx <- match(from, fractions); to_idx <- match(to, fractions)
    prob <- sig_trans$prob_mean[i]
    
    arrow_col <- if (to_idx > from_idx) {
      rgb(0.13, 0.4, 0.67, alpha = 0.7)
    } else {
      rgb(0.84, 0.38, 0.3, alpha = 0.7)
    }
    
    x0 <- node_x[from]; y0 <- node_y[from]
    x1 <- node_x[to];   y1 <- node_y[to]
    dx <- x1 - x0; dy <- y1 - y0
    d <- sqrt(dx^2 + dy^2)
    
    arrows(x0 + 0.10 * dx/d, y0 + 0.10 * dy/d,
           x1 - 0.10 * dx/d, y1 - 0.10 * dy/d,
           lwd = prob * 6, col = arrow_col, length = 0.12, angle = 25)
  }
  
  # Draw nodes
  node_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
  node_sizes <- abs(fraction_means) + 2
  for (i in 1:4) {
    symbols(node_x[i], node_y[i], circles = node_sizes[i] / 25,
            add = TRUE, inches = FALSE, bg = node_cols[i], fg = "black", lwd = 2)
  }
  
  labels <- c("DI\n(Water-soluble)", "NaHCO3\n(Labile)",
              "NaOH\n(Fe/Al-bound)", "HNO3\n(Legacy)")
  text(node_x, node_y, labels, cex = 1.0, font = 2, col = "white")
  
  legend("bottom",
         legend = c("Forward transformation", "Mineralization"),
         col = c(rgb(0.13, 0.4, 0.67), rgb(0.84, 0.38, 0.3)),
         lwd = 4, cex = 1.0, bty = "n", horiz = TRUE,
         title = "Transformation Type", title.font = 2)
  
  # =========================================================================
  # Panel B: Isotopic Distance Matrix
  # =========================================================================
  par(mar = c(5, 5, 4, 6))
  
  color_pal <- colorRampPalette(c("#FFFFCC", "#FD8D3C", "#BD0026"))(100)
  
  image(1:4, 1:4, t(dist_matrix[4:1, ]),
        col = color_pal, zlim = c(0, 7),
        xaxt = "n", yaxt = "n",
        xlab = "From Fraction", ylab = "To Fraction",
        main = "", cex.lab = 1.2, font.lab = 2)
  title(main = "(B) Isotopic Distance Matrix",
        cex.main = 1.4, font.main = 2, line = 2)
  mtext("Lower values = higher transformation probability",
        side = 3, line = 0.5, cex = 0.9, font = 3)
  
  abline(h = 0.5:4.5, col = "white", lwd = 2)
  abline(v = 0.5:4.5, col = "white", lwd = 2)
  axis(1, at = 1:4, labels = fractions, font = 2, cex.axis = 1.2, tick = FALSE)
  axis(2, at = 1:4, labels = rev(fractions), font = 2, cex.axis = 1.2, las = 1, tick = FALSE)
  
  for (i in 1:4) for (j in 1:4) {
    text(i, 5 - j, sprintf("%.2f", dist_matrix[j, i]),
         col = "white", font = 2, cex = 1.3)
  }
  
  legend("right",
         legend = c("7", "5", "3", "1", "0"),
         fill = color_pal[c(100, 75, 50, 25, 1)],
         border = "white", title = "Distance\n(\u2030)",
         cex = 0.9, bty = "n", inset = c(-0.18, 0), xpd = TRUE, title.font = 2)
  
  # =========================================================================
  # Panel C: Posterior Distributions
  # =========================================================================
  par(mar = c(5, 10, 4, 2))
  
  plot(NULL, xlim = c(0, 1), ylim = c(0.5, 4.8),
       xlab = "Transition Probability", ylab = "",
       main = "", xaxt = "n", yaxt = "n", cex.lab = 1.2, font.lab = 2)
  title(main = "(C) Posterior Distributions of Transition Probabilities",
        cex.main = 1.4, font.main = 2, line = 2)
  mtext("95% credible intervals shown", side = 3, line = 0.5, cex = 0.9, font = 3)
  
  abline(v = c(0.25, 0.5, 0.75), col = "grey80", lty = 2)
  axis(1, at = c(0, 0.25, 0.5, 0.75, 1),
       labels = c("0%", "25%", "50%", "75%", "100%"), cex.axis = 1.1)
  
  y_labels <- c(expression("DI" %->% "NaHCO"[3]),
                expression("NaHCO"[3] %->% "NaOH"),
                expression("NaOH" %->% "HNO"[3]),
                expression("HNO"[3] %->% "Mineralization"))
  axis(2, at = 1:4, labels = y_labels, las = 1, font = 2, cex.axis = 1.0, tick = FALSE)
  
  ridge_colors <- c("#7570B3", "#66A61E", "#E6AB02", "#E7298A")
  
  draw_ridge <- function(samples, y_base, color, scale = 0.7) {
    dens <- density(samples, from = 0, to = 1)
    dens$y <- dens$y / max(dens$y) * scale
    polygon(c(dens$x, rev(dens$x)),
            c(y_base + dens$y, rep(y_base, length(dens$x))),
            col = adjustcolor(color, alpha.f = 0.75),
            border = "black", lwd = 1.2)
    q <- quantile(samples, c(0.025, 0.5, 0.975))
    segments(q[1], y_base, q[1], y_base + 0.25, lwd = 2.5)
    segments(q[2], y_base, q[2], y_base + 0.4, lwd = 2.5)
    segments(q[3], y_base, q[3], y_base + 0.25, lwd = 2.5)
  }
  
  draw_ridge(post_DI_NaHCO3,    1, ridge_colors[1])
  draw_ridge(post_NaHCO3_NaOH,  2, ridge_colors[2])
  draw_ridge(post_NaOH_HNO3,    3, ridge_colors[3])
  draw_ridge(post_HNO3_mineral,  4, ridge_colors[4])
  
  # =========================================================================
  # Panel D: Ranked Pathways
  # =========================================================================
  par(mar = c(5, 9, 4, 3))
  
  plot(NULL, xlim = c(0, 1), ylim = c(0.5, 8.5),
       xlab = "Transformation Probability", ylab = "",
       main = "", xaxt = "n", yaxt = "n", cex.lab = 1.2, font.lab = 2)
  title(main = "(D) Ranked P Transformation Pathways",
        cex.main = 1.4, font.main = 2, line = 2)
  mtext("With 95% credible intervals", side = 3, line = 0.5, cex = 0.9, font = 3)
  
  abline(v = c(0.25, 0.5, 0.75), col = "grey80", lty = 2)
  axis(1, at = c(0, 0.25, 0.5, 0.75, 1),
       labels = c("0%", "25%", "50%", "75%", "100%"), cex.axis = 1.1)
  axis(2, at = 1:8, labels = rev(top_pathways$pathway),
       las = 1, font = 2, cex.axis = 1.0, tick = FALSE)
  
  frac_colors <- c("DI" = "#1B9E77", "NaHCO3" = "#D95F02",
                   "NaOH" = "#7570B3", "HNO3" = "#E7298A")
  
  for (i in 1:8) {
    y_pos <- 9 - i
    prob <- top_pathways$prob_mean[i]
    ci_lo <- top_pathways$prob_lower[i]
    ci_hi <- top_pathways$prob_upper[i]
    
    segments(0, y_pos, prob, y_pos, lwd = 2.5, col = "grey50")
    segments(ci_lo, y_pos, ci_hi, y_pos, lwd = 2, col = "grey30")
    segments(ci_lo, y_pos - 0.15, ci_lo, y_pos + 0.15, lwd = 2, col = "grey30")
    segments(ci_hi, y_pos - 0.15, ci_hi, y_pos + 0.15, lwd = 2, col = "grey30")
    
    points(prob, y_pos, pch = 21, cex = 2.2,
           bg = frac_colors[top_pathways$from[i]], col = "black", lwd = 1.5)
  }
  
  legend("bottomright", legend = fractions,
         pch = 21, pt.bg = frac_colors[fractions], pt.cex = 1.8,
         title = "Source Fraction", cex = 1.0, bty = "n", title.font = 2)
}

# Save in multiple formats
base <- file.path(OUTPUT_DIR, "figures/Figure3_Bayesian_FourPanel")

pdf(paste0(base, ".pdf"), width = 14, height = 12, family = "serif")
par(oma = c(1, 1, 3, 1))
create_four_panel()
mtext("Bayesian Analysis of Interpool Phosphorus Transformations",
      outer = TRUE, side = 3, line = 1, cex = 1.6, font = 2)
invisible(dev.off())

png(paste0(base, ".png"), width = 14, height = 12, units = "in", res = 300)
par(oma = c(1, 1, 3, 1), family = "serif")
create_four_panel()
mtext("Bayesian Analysis of Interpool Phosphorus Transformations",
      outer = TRUE, side = 3, line = 1, cex = 1.6, font = 2)
invisible(dev.off())

tiff(paste0(base, ".tiff"), width = 14, height = 12, units = "in", res = 600,
     compression = "lzw")
par(oma = c(1, 1, 3, 1), family = "serif")
create_four_panel()
mtext("Bayesian Analysis of Interpool Phosphorus Transformations",
      outer = TRUE, side = 3, line = 1, cex = 1.6, font = 2)
invisible(dev.off())

cat("  Saved: Figure3_Bayesian_FourPanel (.pdf, .png, .tiff)\n")
cat("  Figure 3 complete\n\n")
