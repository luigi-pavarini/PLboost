# Regenerate the lollipop selection plots at a WIDE/SHORT aspect (2000x1300)
# with fewer entries (top_n = 8) so all six families fit on one condensed page.
# Overwrites imgs/<fam>/lollipop_mu.png and lollipop_sigma.png only.

save_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must hold PL_simulations/ and IPL_simulations/ with the .RData
families  <- c("NO", "TF", "LO", "PE", "Hyp", "SLASH")
M <- 100
TOPN <- 8
LOL_W <- 2000; LOL_H <- 2000

lollipop_sel <- function(freq_mat, top_n = TOPN, param_label = "mu", sd_mult = 1) {
  # Proportion of replications in which each covariate was selected at least once.
  # (Previously this normalized each run's selection COUNTS to sum to 1, which is a
  #  share of the boosting update budget, not a selection frequency.)
  sel_mat  <- (freq_mat > 0)[-1, , drop = FALSE]
  avg_prop <- rowMeans(sel_mat)
  sd_prop  <- sqrt(avg_prop * (1 - avg_prop) / M)   # binomial standard error
  ord       <- order(avg_prop, decreasing = TRUE)[seq_len(top_n)]
  top_avg   <- avg_prop[ord]; top_sd <- sd_prop[ord]; top_names <- names(top_avg)
  others_avg <- mean(avg_prop[-ord]); others_sd <- sd(avg_prop[-ord])
  p_plot <- c(top_avg, others_avg); sd_plot <- c(top_sd, others_sd)
  lbl_plot <- c(top_names, "others"); n_plot <- top_n + 1
  xlab_expr <- switch(param_label,
    mu    = expression("Proportion of replications selected" ~ eta[mu]),
    sigma = expression("Proportion of replications selected" ~ eta[sigma]))
  old_par <- par(mar = c(4, 5, 1, 1))
  x_max <- max(p_plot + sd_mult * sd_plot, na.rm = TRUE) * 1.15
  plot(p_plot, seq_len(n_plot), ylim = c(n_plot, 1), xlim = c(0, x_max),
       yaxt = "n", xlab = xlab_expr, ylab = "", main = "", pch = 19,
       col = "gray30", cex = 1.1)
  segments(x0 = 0, x1 = p_plot, y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "gray80", lwd = 0.8, lty = 2)
  x_lo <- pmax(p_plot - sd_mult * sd_plot, 0); x_hi <- p_plot + sd_mult * sd_plot
  cap <- 0.18
  segments(x0 = x_lo, x1 = x_hi, y0 = seq_len(n_plot), y1 = seq_len(n_plot), col = "#4682B4", lwd = 1.4)
  segments(x0 = x_lo, x1 = x_lo, y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap, col = "#4682B4", lwd = 1.4)
  segments(x0 = x_hi, x1 = x_hi, y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap, col = "#4682B4", lwd = 1.4)
  axis(2, at = seq_len(n_plot), labels = lbl_plot, las = 1, cex.axis = 0.9)
  abline(v = 0, col = "gray60", lty = 2); par(old_par)
}

for (fam in families) {
  rdata <- file.path(save_path, "PL_simulations",
                     paste0("sim_PL_M100_mstop100_", fam, ".RData"))
  if (!file.exists(rdata)) { cat("SKIP:", rdata, "\n"); next }
  load(rdata)  # -> sim
  img_dir <- file.path(save_path, "PL_simulations", "imgs", fam)
  dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)
  save_lol <- function(name, expr) {
    png(file.path(img_dir, paste0(name, ".png")), width = LOL_W, height = LOL_H, res = 330)
    on.exit(dev.off()); expr
  }
  save_lol("lollipop_mu",    lollipop_sel(sim$freq_sel_mu,    param_label = "mu"))
  save_lol("lollipop_sigma", lollipop_sel(sim$freq_sel_sigma, param_label = "sigma"))
  cat("wide lollipops for", fam, "done\n")
}
cat("DONE.\n")
