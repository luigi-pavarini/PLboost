# Regenerate the IPL simulation figures from the CURRENT (medium) scenario.
# Two fixes vs the figures previously used in the dissertation:
#   (1) they came from the superseded `low` scenario (median exp(lambda) = 1.47)
#       while the text/Appendix A.2 report the `medium` scenario (2.4242);
#   (2) the lollipop metric normalized each run's selection counts to sum to 1
#       (a share of the boosting update budget), instead of the proportion of
#       replications in which the covariate was selected.
# No refitting: this only replots the stored RData.

save_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must hold PL_simulations/ and IPL_simulations/ with the .RData
fam <- "NO"; M <- 100
beta_alpha <- c(0.70, -0.50)

load(file.path(save_path, "IPL_simulations", "medium_sim_IPL_M100_mstop500_NO.RData"))  # -> sim

img_dir <- file.path(save_path, "IPL_simulations", "imgs", fam)
if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE)
save_png <- function(name, expr, w = 2000, h = 2000) {
  png(file.path(img_dir, paste0(name, ".png")), width = w, height = h, res = 330)
  on.exit(dev.off()); expr
}

save_png("mu_coefs", {
  boxplot(sim$est_mu[2,], sim$est_mu[3,], sim$est_mu[4,], sim$est_mu[5,],
          c(sim$est_mu[c(6:1001),]), ylim = c(-1.3, 0.6), xaxt = "n", las = 1,
          ylab = expression("Coefficients for" ~ eta[mu]), main = "")
  abline(h = 0, lty = 2, col = "gray")
  axis(1, 1, expression(X[1])); axis(1, 2, expression(X[2]))
  axis(1, 3, expression(X[3])); axis(1, 4, expression(X[4])); axis(1, 5, "non-inf.")
  segments(0.5, -0.99, 1.5, -0.99, lty = 2, col = "red4")
  segments(1.5,  0.49, 2.5,  0.49, lty = 2, col = "red4")
  segments(2.5, -1.17, 3.5, -1.17, lty = 2, col = "red4")
  segments(3.5,  0.38, 4.5,  0.38, lty = 2, col = "red4")
  segments(4.5,  0,    5.5,  0,    lty = 2, col = "red4")
})

save_png("sigma_coefs", {
  boxplot(sim$est_sigma[4,], sim$est_sigma[5,], sim$est_sigma[6,],
          c(sim$est_sigma[c(2:3, 7:1001),]), ylim = c(-0.6, 1), xaxt = "n", las = 1,
          ylab = expression("Coefficients for" ~ eta[sigma]), main = "")
  abline(h = 0, lty = 2, col = "gray")
  axis(1, 1, expression(X[3])); axis(1, 2, expression(X[4]))
  axis(1, 3, expression(X[5])); axis(1, 4, "non-inf.")
  segments(0.5, -0.41, 1.5, -0.41, lty = 2, col = "red4")
  segments(1.5, -0.24, 2.5, -0.24, lty = 2, col = "red4")
  segments(2.5,  0.81, 3.5,  0.81, lty = 2, col = "red4")
  segments(3.5,  0,    4.5,  0,    lty = 2, col = "red4")
})

save_png("lambda_coefs", {
  lam <- exp(sim$est_lambda[1,]); r <- range(lam); pad <- 0.6 * diff(r) + 0.05
  boxplot(lam, ylim = c(min(r[1] - pad, 1.1), r[2] + pad), xaxt = "n", las = 1,
          ylab = expression("Coefficients for" ~ lambda), main = "")
  axis(1, 1, "exp(intercept)")
  segments(0.5, 1.2, 1.5, 1.2, lty = 2, col = "red4")
})

save_png("alpha_coefs", {
  boxplot(sim$est_alpha[1000,], sim$est_alpha[1001,], c(sim$est_alpha[c(2:999),]),
          ylim = c(-1, 2), xaxt = "n", las = 1,
          ylab = expression("Coefficients for" ~ eta[alpha]), main = "")
  abline(h = 0, lty = 2, col = "gray")
  axis(1, 1, expression(X[999])); axis(1, 2, expression(X[1000])); axis(1, 3, "non-inf.")
  segments(0.5, beta_alpha[1], 1.5, beta_alpha[1], lty = 2, col = "red4")
  segments(1.5, beta_alpha[2], 2.5, beta_alpha[2], lty = 2, col = "red4")
  segments(2.5, 0,             3.5, 0,             lty = 2, col = "red4")
})

# ---- Lollipops: proportion of replications in which the covariate was selected ----
lollipop_sel <- function(freq_mat, top_n = 8, param_label = "mu") {
  sel_mat  <- (freq_mat > 0)[-1, , drop = FALSE]
  avg_prop <- rowMeans(sel_mat)
  sd_prop  <- sqrt(avg_prop * (1 - avg_prop) / M)     # binomial standard error
  ord <- order(avg_prop, decreasing = TRUE)[seq_len(top_n)]
  p_plot   <- c(avg_prop[ord], mean(avg_prop[-ord]))
  sd_plot  <- c(sd_prop[ord],  sd(avg_prop[-ord]))
  lbl_plot <- c(names(avg_prop)[ord], "others"); n_plot <- top_n + 1
  xlab_expr <- switch(param_label,
    mu    = expression("Proportion of replications selected" ~ eta[mu]),
    sigma = expression("Proportion of replications selected" ~ eta[sigma]),
    alpha = expression("Proportion of replications selected" ~ eta[alpha]))
  op <- par(mar = c(4, 5, 1, 1)); xm <- max(p_plot + sd_plot, na.rm = TRUE) * 1.15
  plot(p_plot, seq_len(n_plot), ylim = c(n_plot, 1), xlim = c(0, xm), yaxt = "n",
       xlab = xlab_expr, ylab = "", main = "", pch = 19, col = "gray30", cex = 1.1)
  segments(0, seq_len(n_plot), p_plot, seq_len(n_plot), col = "gray80", lwd = 0.8, lty = 2)
  xlo <- pmax(p_plot - sd_plot, 0); xhi <- p_plot + sd_plot; cap <- 0.18
  segments(xlo, seq_len(n_plot), xhi, seq_len(n_plot), col = "#4682B4", lwd = 1.4)
  segments(xlo, seq_len(n_plot) - cap, xlo, seq_len(n_plot) + cap, col = "#4682B4", lwd = 1.4)
  segments(xhi, seq_len(n_plot) - cap, xhi, seq_len(n_plot) + cap, col = "#4682B4", lwd = 1.4)
  axis(2, at = seq_len(n_plot), labels = lbl_plot, las = 1, cex.axis = 0.9)
  abline(v = 0, col = "gray60", lty = 2); par(op)
}

save_png("lollipop_mu",    lollipop_sel(sim$freq_sel_mu,    param_label = "mu"),    w = 2000, h = 2000)
save_png("lollipop_sigma", lollipop_sel(sim$freq_sel_sigma, param_label = "sigma"), w = 2000, h = 2000)
save_png("lollipop_alpha", lollipop_sel(sim$freq_sel_alpha, param_label = "alpha"), w = 2000, h = 2000)

cat("median exp(lambda) =", round(median(exp(sim$est_lambda[1,])), 4), "\n")
for (nm in c("mu","sigma","alpha")) {
  fm <- sim[[paste0("freq_sel_", nm)]]
  p <- sort(rowMeans((fm > 0)[-1,,drop=FALSE]), decreasing = TRUE)[1:4]
  cat(nm, ":", paste(names(p), round(p,2), sep="=", collapse="  "), "\n")
}
cat("DONE ->", img_dir, "\n")
