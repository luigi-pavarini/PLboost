# ==============================================================================
# credit_img_beoi_cv.R
# Visualization for the Credit LGD BEOI application — CV version.
# Reads 5-fold CV outputs (sim_results_credit_VT_cv_BEOI_<batch>.RData) via
# load_cv_results(), which globs and merges EVERY batch tag (_1, _2, _3, ...);
# credit's three batches (repeats 1-3 / 4-6 / 7-10) are pooled automatically, so
# the new _3 needs NO change here.
# Adapts to the number of repeats: each entry is one (repeat, fold) result;
# denominators use n_entries (5 if single pass, k*R if repeated). No lambda/zeta.
# ==============================================================================

fam          <- "Beta"   # used only for the imgs/<fam> output folder name
# --- resolucao de caminho relativa a este arquivo (nao depende de outras pastas) --
.script_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(dirname(normalizePath(o)))
  }
  a  <- commandArgs(trailingOnly = FALSE)
  sp <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(sp)) return(dirname(normalizePath(sp)))
  getwd()
}

app_path     <- Sys.getenv("PLBOOST_OUT", unset = ".")   # set PLBOOST_OUT (must contain essentials/ and data/)
results_path <- file.path(app_path, "cross_val", "results")
base_path    <- file.path(app_path, "cross_val", "imgs", fam)
if (!dir.exists(base_path)) dir.create(base_path, recursive = TRUE)

# ---- Image export settings --------------------------------------------------
img_width  <- 2000
img_height <- 2000
img_res    <- 330
rel_dim    <- 1.8

source(file.path(app_path, "essentials", "cv_load.R"))

cv          <- load_cv_results(results_path, "sim_results_credit_VT", "BEOI")
sim_results <- cv$sim_results
var_dict    <- cv$var_dict

# ── FILTER NULL ENTRIES ───────────────────────────────────────────────────────
success_idx <- which(!sapply(sim_results, is.null))
sim_results <- sim_results[success_idx]
n_success   <- length(sim_results)   # = k * R that succeeded (5 if single pass)
cat("Successful CV entries (folds x repeats):", n_success, "\n")

# ==============================================================================
# HELPERS
# ==============================================================================

map_names <- function(vnames, dict) {
  idx <- match(vnames, dict$new_name)
  ifelse(is.na(idx), vnames, dict$orig_name[idx])
}

# Reconstruct y_all
all_test_df <- do.call(rbind, lapply(sim_results, function(x) {
  data.frame(idx = x$test_idx, y = x$y_test)
}))
all_test_df <- all_test_df[!duplicated(all_test_df$idx), ]
y_all       <- all_test_df$y

# Derived quantities
errors_beoi <- sapply(sim_results, function(x) x$error_beoi) * 100
mae_beoi    <- sapply(sim_results, function(x) x$mae_beoi)
upsilon_vec <- sapply(sim_results, function(x) x$upsilon)
all_ytest   <- unlist(lapply(sim_results, function(x) x$y_test))
all_muhat   <- unlist(lapply(sim_results, function(x) x$mu_hat))
all_resid   <- all_ytest - all_muhat
all_nu_hat  <- unlist(lapply(sim_results, function(x) x$nu_hat))

# Coefficient matrices
build_coef_mat <- function(coef_list, exclude = "(Intercept)") {
  all_vars <- unique(unlist(lapply(coef_list, function(cv) {
    names(cv)[names(cv) != exclude]
  })))
  mat <- matrix(0, nrow = length(all_vars), ncol = length(coef_list),
                dimnames = list(all_vars, NULL))
  for (r in seq_along(coef_list)) {
    cv  <- coef_list[[r]]
    nms <- names(cv)[names(cv) != exclude]
    mat[nms, r] <- cv[nms]
  }
  mat
}

coef_mu_mat    <- build_coef_mat(lapply(sim_results, function(x) x$coef_mu))
coef_sigma_mat <- build_coef_mat(lapply(sim_results, function(x) x$coef_sigma))
coef_nu_mat    <- build_coef_mat(lapply(sim_results, function(x) x$coef_nu))

# ==============================================================================
# PLOT 1: LGD response distribution
# ==============================================================================

png(file.path(base_path, "plot_response_dist.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
prop_one   <- mean(y_all == 1)
dens_obj   <- density(y_all[y_all < 1], bw = "SJ")
dens_x     <- dens_obj$x[dens_obj$x < 1]
dens_y     <- dens_obj$y[dens_obj$x < 1]
pin_x      <- dens_x[which.max(dens_x)]
pin_height <- dens_y[which.max(dens_x)]
hist(y_all, breaks = 25, freq = FALSE,
     col = "gray85", border = "white",
     xlab = "Loss Given Default (LGD)",
     ylab = "Density", main = "", las = 1,
     ylim = c(0, 10), xlim = c(0, 1))
lines(dens_x, dens_y * (1 - prop_one) * (1 / max(dens_obj$x)),
      col = "gray30", lwd = 2)
segments(x0 = pin_x, x1 = pin_x, y0 = 0, y1 = pin_height,
         col = "gray30", lwd = 2)
points(x = pin_x, y = pin_height, pch = 21,
       col = "gray30", bg = "gray30", cex = 1.6)
abline(v = mean(y_all), lty = 3, col = "gray40", lwd = 1.5)
legend("topleft",
       legend = c("Kernel density (lgd < 1)",
                  paste0("lgd = 1 mass: ", round(prop_one * 100, 1), "%"),
                  "Mean"),
       col = c("gray30", "gray30", "gray40"),
       lty = c(1, NA, 3), lwd = c(2, NA, 1.5),
       pch = c(NA, 21, NA), pt.bg = c(NA, "gray30", NA), pt.cex = 1.3,
       bty = "n", cex = 0.85)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 2: MSE distribution
# ==============================================================================

png(file.path(base_path, "plot_mse_dist.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
boxplot(errors_beoi, ylab = expression("MSE" %*% 100),
        xaxt = "n", col = "gray85", border = "gray30",
        las = 1, outline = TRUE, main = "")
axis(1, at = 1, labels = "BEOI")
abline(h = mean(errors_beoi), lty = 2, col = "red4", lwd = 1.5)
legend("topright",
       legend = paste0("Mean = ", round(mean(errors_beoi), 3),
                       "   SD = ", round(sd(errors_beoi), 3)),
       bty = "n", cex = 0.85)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 2b: MAE distribution
# ==============================================================================

png(file.path(base_path, "plot_mae_dist.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
boxplot(mae_beoi, ylab = "MAE",
        xaxt = "n", col = "gray85", border = "gray30",
        las = 1, outline = TRUE, main = "")
axis(1, at = 1, labels = "BEOI")
abline(h = mean(mae_beoi), lty = 2, col = "red4", lwd = 1.5)
legend("topright",
       legend = paste0("Mean = ", round(mean(mae_beoi), 4),
                       "   SD = ", round(sd(mae_beoi), 4)),
       bty = "n", cex = 0.85)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 2c: Upsilon distribution
# ==============================================================================

png(file.path(base_path, "plot_upsilon_dist.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
boxplot(upsilon_vec, ylab = expression(Upsilon),
        xaxt = "n", col = "gray85", border = "gray30",
        las = 1, outline = TRUE, main = "")
axis(1, at = 1, labels = "BEOI")
abline(h = mean(upsilon_vec), lty = 2, col = "red4", lwd = 1.5)
legend("topright",
       legend = paste0("Mean = ", round(mean(upsilon_vec), 4),
                       "   SD = ", round(sd(upsilon_vec), 4)),
       bty = "n", cex = 0.85)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 3: Observed vs predicted — overlapping densities
# ==============================================================================

png(file.path(base_path, "plot_pred_vs_obs.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
d_obs  <- density(all_ytest, bw = 0.02)
d_pred <- density(all_muhat, bw = 0.02)
xlim_r <- range(c(d_obs$x, d_pred$x))
ylim_r <- c(0, max(c(d_obs$y, d_pred$y)) * 1.1)
plot(d_obs, col = "gray30", lwd = 2,
     xlim = xlim_r, ylim = ylim_r,
     xlab = "LGD", ylab = "Density", main = "", las = 1)
polygon(c(d_obs$x, rev(d_obs$x)),
        c(d_obs$y, rep(0, length(d_obs$y))),
        col = adjustcolor("gray30", alpha.f = 0.15), border = NA)
lines(d_pred, col = "red4", lwd = 2, lty = 2)
polygon(c(d_pred$x, rev(d_pred$x)),
        c(d_pred$y, rep(0, length(d_pred$y))),
        col = adjustcolor("red4", alpha.f = 0.15), border = NA)
abline(v = mean(all_ytest), lty = 3, col = "gray30", lwd = 1.2)
abline(v = mean(all_muhat), lty = 3, col = "red4",   lwd = 1.2)
legend("topleft",
       legend = c(expression(y[test]), expression(hat(mu))),
       col = c("gray30", "red4"), lty = c(1, 2), lwd = c(2, 2),
       fill = c(adjustcolor("gray30", 0.15), adjustcolor("red4", 0.15)),
       border = NA, bty = "n", cex = 0.9)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 4a: Residuals — density
# ==============================================================================

png(file.path(base_path, "plot_residuals_density.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
hist(all_resid, breaks = 50, freq = FALSE,
     col = "gray85", border = "white",
     xlab = expression(y[test] - hat(mu)),
     ylab = "Density", main = "", las = 1)
lines(density(all_resid), col = "gray30", lwd = 2)
abline(v = 0, col = "red4", lty = 2, lwd = 1.5)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 4b: Residuals — boxplot
# ==============================================================================

png(file.path(base_path, "plot_residuals_boxplot.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
boxplot(all_resid, ylab = expression(y[test] - hat(mu)),
        xaxt = "n", col = "gray85", border = "gray30",
        outline = TRUE, las = 1, main = "")
abline(h = 0, col = "red4", lty = 2, lwd = 1.5)
axis(1, at = 1, labels = "Residuals")
par(old_par)
dev.off()

# ==============================================================================
# PLOT 5: Nu_hat distribution — clipped KDE
# ==============================================================================

png(file.path(base_path, "plot_nu_dist.png"),
    width = img_width, height = img_height, res = img_res)
old_par  <- par(mar = c(4.5, 4.5, 2, 2))
dens_nu  <- density(all_nu_hat, bw = "SJ")
in_range <- dens_nu$x >= 0 & dens_nu$x <= 1
dx       <- dens_nu$x[in_range]
dy       <- dens_nu$y[in_range]
ylim_top <- max(dy) * 1.10
hist(all_nu_hat, breaks = 20, freq = FALSE,
     col = "gray85", border = "white",
     xlab = expression("Predicted one-inflation probability  " * hat(nu)),
     ylab = "Density", main = "", las = 1,
     xlim = c(0, 1), ylim = c(0, ylim_top))
lines(dx, dy, col = "gray30", lwd = 2)
abline(v = mean(all_nu_hat),   lty = 2, col = "red4",   lwd = 1.5)
abline(v = median(all_nu_hat), lty = 3, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c("Kernel density",
                  paste0("Mean = ",   round(mean(all_nu_hat), 3)),
                  paste0("Median = ", round(median(all_nu_hat), 3))),
       col = c("gray30", "red4", "gray40"),
       lty = c(1, 2, 3), lwd = c(2, 1.5, 1.5),
       bty = "n", cex = 0.85)
par(old_par)
dev.off()

# ==============================================================================
# PLOT 6: Nu_hat vs mu_hat — loess + SE band
# ==============================================================================

png(file.path(base_path, "plot_nu_vs_mu.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(5.5, 4.5, 2, 2))
lo      <- loess(all_nu_hat ~ all_muhat, span = 0.3)
mu_seq  <- seq(min(all_muhat), quantile(all_muhat, 0.98), length.out = 500)
lo_pred <- predict(lo, newdata = data.frame(all_muhat = mu_seq), se = TRUE)
lo_fit  <- lo_pred$fit
lo_hi   <- lo_fit + lo_pred$se.fit
lo_lo   <- lo_fit - lo_pred$se.fit
min_idx <- which.min(lo_fit)
min_mu  <- mu_seq[min_idx]
plot(all_muhat, all_nu_hat,
     pch = 19, cex = 0.3,
     col = adjustcolor("gray30", alpha.f = 0.15),
     xlab = expression("Predicted mean  " * hat(mu)),
     ylab = expression("Predicted one-inflation probability  " * hat(nu)),
     las = 1, main = "")
polygon(c(mu_seq, rev(mu_seq)),
        c(lo_hi,  rev(lo_lo)),
        col = adjustcolor("red4", alpha.f = 0.15), border = NA)
lines(mu_seq, lo_fit, col = "red4", lwd = 2)
abline(v = min_mu, lty = 3, col = "gray40", lwd = 1.2)
text(min_mu, par("usr")[4] * 0.97,
     labels = bquote(hat(mu) == .(round(min_mu, 2))),
     adj = c(-0.1, 1), col = "gray40", cex = 0.75)
rug(all_muhat, col = adjustcolor("gray30", alpha.f = 0.05), ticksize = 0.02)
legend("topleft",
       legend = c("Observations", "Loess trend", "\u00b11 SE band"),
       col    = c(adjustcolor("gray30", 0.5), "red4", adjustcolor("red4", 0.4)),
       pch    = c(19, NA, NA), lty = c(NA, 1, NA), lwd = c(NA, 2, NA),
       fill   = c(NA, NA, adjustcolor("red4", 0.15)),
       border = c(NA, NA, NA),
       bty    = "n", cex = 0.85)
par(old_par)
dev.off()

# ==============================================================================
# PLOTS 7-9: Coefficient forest plots — mu, sigma, nu
# ==============================================================================

plot_coef_forest <- function(coef_mat, var_dict, param_label = "mu",
                             min_freq = 0.10, filename,
                             legend_pos = "bottomright") {
  avg  <- rowMeans(coef_mat)
  sdev <- apply(coef_mat, 1, sd)
  freq <- rowMeans(coef_mat != 0)
  keep <- freq >= min_freq
  avg  <- avg[keep];  sdev <- sdev[keep]
  top10 <- order(abs(avg), decreasing = TRUE)[seq_len(min(10, length(avg)))]
  avg   <- avg[top10];  sdev <- sdev[top10]
  ord  <- order(abs(avg), decreasing = FALSE)
  avg  <- avg[ord];   sdev <- sdev[ord]
  orig <- map_names(names(avg), var_dict)
  n    <- length(avg)
  ci68_lo <- avg - 1.000 * sdev
  ci68_hi <- avg + 1.000 * sdev
  ci90_lo <- avg - 1.645 * sdev
  ci90_hi <- avg + 1.645 * sdev
  ci95_lo <- avg - 1.960 * sdev
  ci95_hi <- avg + 1.960 * sdev
  xlab_expr <- switch(param_label,
    mu    = expression("Average coefficient  " * eta[mu]),
    sigma = expression("Average coefficient  " * eta[sigma]),
    nu    = expression("Average coefficient  " * eta[nu])
  )
  rect_half  <- 0.30
  col_pos    <- adjustcolor("#4682B4", alpha.f = 0.45)
  col_neg    <- adjustcolor("#B44646", alpha.f = 0.45)
  col_pos_90 <- adjustcolor("#4682B4", alpha.f = 0.80)
  col_neg_90 <- adjustcolor("#B44646", alpha.f = 0.80)
  col_bord   <- "#2a2a2a"
  x_range    <- range(c(ci95_lo, ci95_hi, 0)) * 1.15
  x_range[1] <- min(x_range[1], -0.01)
  png(filename, width = img_width * rel_dim,
      height = max(img_height, 55 * n + 150), res = img_res)
  old_par <- par(mar = c(4.5, 11, 2, 3))
  plot(NA, xlim = x_range, ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = xlab_expr, ylab = "", main = "", las = 1)
  abline(v = 0, col = "gray60", lty = 2)
  rect_col    <- ifelse(avg >= 0, col_pos,    col_neg)
  wick_col_90 <- ifelse(avg >= 0, col_pos_90, col_neg_90)
  cap_95      <- rect_half * 0.75
  segments(x0 = ci95_lo, x1 = ci95_hi,
           y0 = seq_len(n), y1 = seq_len(n),
           col = col_bord, lwd = 1.4, lty = 3)
  segments(x0 = ci95_lo, x1 = ci95_lo,
           y0 = seq_len(n) - cap_95, y1 = seq_len(n) + cap_95,
           col = col_bord, lwd = 1.4, lty = 3)
  segments(x0 = ci95_hi, x1 = ci95_hi,
           y0 = seq_len(n) - cap_95, y1 = seq_len(n) + cap_95,
           col = col_bord, lwd = 1.4, lty = 3)
  segments(x0 = ci90_lo, x1 = ci90_hi,
           y0 = seq_len(n), y1 = seq_len(n),
           col = wick_col_90, lwd = 1.4, lty = 1)
  rect(xleft   = ci68_lo, ybottom = seq_len(n) - rect_half,
       xright  = ci68_hi, ytop    = seq_len(n) + rect_half,
       col = rect_col, border = col_bord, lwd = 0.6)
  points(avg, seq_len(n), pch = 22, col = col_bord,
         bg = "white", cex = 1.3, lwd = 1.2)
  axis(2, at = seq_len(n), labels = orig, las = 1, cex.axis = 0.78)
  legend(legend_pos,
         legend = c("\u00b11 SD body (positive)", "\u00b11 SD body (negative)",
                    "\u00b11.645 SD solid (positive)", "\u00b11.645 SD solid (negative)",
                    "\u00b11.960 SD dotted T", "Mean"),
         fill   = c(col_pos, col_neg, NA, NA, NA, NA),
         border = c(col_bord, col_bord, NA, NA, NA, NA),
         lty    = c(NA, NA, 1, 1, 3, NA),
         lwd    = c(NA, NA, 1.4, 1.4, 1.4, NA),
         col    = c(NA, NA, col_pos_90, col_neg_90, col_bord, col_bord),
         pch    = c(NA, NA, NA, NA, NA, 22),
         pt.bg  = c(NA, NA, NA, NA, NA, "white"),
         pt.cex = 1.3, bty = "n", cex = 0.80)
  par(old_par)
  dev.off()
}

# ==============================================================================
# PLOTS 10-12: Selection lollipops — mu, sigma, nu
# ==============================================================================

lollipop_sel_app <- function(coef_mat, var_dict, top_n = 15,
                             param_label = "mu", sd_mult = 1, filename) {
  col_totals <- colSums(coef_mat != 0)
  col_totals[col_totals == 0] <- 1
  prop_mat   <- sweep(coef_mat != 0, 2, col_totals, FUN = "/")
  M          <- ncol(prop_mat)
  avg_prop   <- rowSums(prop_mat) / M
  sd_prop    <- apply(prop_mat, 1, sd)
  top_n_use  <- min(top_n, length(avg_prop))
  ord        <- order(avg_prop, decreasing = TRUE)[seq_len(top_n_use)]
  top_avg    <- avg_prop[ord]
  top_sd     <- sd_prop[ord]
  top_names  <- map_names(names(top_avg), var_dict)
  others_avg <- sum(avg_prop[-ord])
  others_sd  <- sqrt(sum(sd_prop[-ord]^2))
  p_plot     <- c(top_avg, others_avg)
  sd_plot    <- c(top_sd,  others_sd)
  lbl_plot   <- c(top_names, "others")
  n_plot     <- length(p_plot)
  xlab_expr <- switch(param_label,
    mu    = expression("Avg. selection proportion across CV entries" ~ eta[mu]),
    sigma = expression("Avg. selection proportion across CV entries" ~ eta[sigma]),
    nu    = expression("Avg. selection proportion across CV entries" ~ eta[nu])
  )
  png(filename, width = img_width * rel_dim, height = img_height, res = img_res)
  old_par <- par(mar = c(4, 12, 2, 2))
  x_max   <- max(p_plot + sd_mult * sd_plot, na.rm = TRUE) * 1.15
  plot(p_plot, seq_len(n_plot),
       ylim = c(n_plot, 1), xlim = c(0, x_max),
       yaxt = "n", xlab = xlab_expr, ylab = "", main = "",
       pch = 19, col = "gray30", cex = 1.3)
  segments(x0 = 0, x1 = p_plot,
           y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "gray80", lwd = 0.8, lty = 2)
  x_lo <- pmax(p_plot - sd_mult * sd_plot, 0)
  x_hi <- p_plot + sd_mult * sd_plot
  cap  <- 0.18
  segments(x0 = x_lo, x1 = x_hi,
           y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "#4682B4", lwd = 1.4)
  segments(x0 = x_lo, x1 = x_lo,
           y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap,
           col = "#4682B4", lwd = 1.4)
  segments(x0 = x_hi, x1 = x_hi,
           y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap,
           col = "#4682B4", lwd = 1.4)
  axis(2, at = seq_len(n_plot), labels = lbl_plot, las = 1, cex.axis = 0.80)
  abline(v = 0, col = "gray60", lty = 2)
  par(old_par)
  dev.off()
}

# ---- mu ---------------------------------------------------------------------
plot_coef_forest(coef_mu_mat, var_dict, "mu",
                 filename = file.path(base_path, "plot_coef_mu.png"),
                 legend_pos = "bottomleft")
lollipop_sel_app(coef_mu_mat, var_dict, top_n = 15, param_label = "mu",
                 filename = file.path(base_path, "plot_sel_mu.png"))

# ---- sigma ------------------------------------------------------------------
plot_coef_forest(coef_sigma_mat, var_dict, "sigma",
                 filename = file.path(base_path, "plot_coef_sigma.png"))
lollipop_sel_app(coef_sigma_mat, var_dict, top_n = 15, param_label = "sigma",
                 filename = file.path(base_path, "plot_sel_sigma.png"))

# ---- nu ---------------------------------------------------------------------
plot_coef_forest(coef_nu_mat, var_dict, "nu",
                 filename = file.path(base_path, "plot_coef_nu.png"))
lollipop_sel_app(coef_nu_mat, var_dict, top_n = 15, param_label = "nu",
                 filename = file.path(base_path, "plot_sel_nu.png"))

# ==============================================================================
# PRINTED SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  CREDIT LGD APPLICATION (Vermont) — BEOI SUMMARY\n")
cat("================================================================\n\n")

cat("--- Dataset (Vermont only) ---\n")
n_full <- length(y_all)
cat("  n (total observations) :", n_full, "\n")
cat("  p (predictors)         :", nrow(var_dict), "\n")
cat("  lgd range              : [", round(min(y_all), 4), ",",
    round(max(y_all), 4), "]\n")
cat("  lgd mean               :", round(mean(y_all), 4), "\n")
cat("  lgd median             :", round(median(y_all), 4), "\n")
cat("  lgd SD                 :", round(sd(y_all), 4), "\n")
cat("  Proportion lgd = 1     :", round(mean(y_all == 1) * 100, 1), "%\n\n")

cat("--- Prediction performance (", n_success, "CV entries) ---\n")
cat("  Mean MSE x 100         :", round(mean(errors_beoi), 3), "\n")
cat("  SD   MSE x 100         :", round(sd(errors_beoi), 3), "\n")
cat("  Min  MSE x 100         :", round(min(errors_beoi), 3), "\n")
cat("  Max  MSE x 100         :", round(max(errors_beoi), 3), "\n\n")

cat("--- MAE (test set) ---\n")
cat("  Mean MAE               :", round(mean(mae_beoi), 4), "\n")
cat("  SD   MAE               :", round(sd(mae_beoi), 4), "\n")
cat("  Min  MAE               :", round(min(mae_beoi), 4), "\n")
cat("  Max  MAE               :", round(max(mae_beoi), 4), "\n\n")

cat("--- Upsilon (train set, continuous part) ---\n")
cat("  Mean Upsilon           :", round(mean(upsilon_vec), 4), "\n")
cat("  SD   Upsilon           :", round(sd(upsilon_vec), 4), "\n")
cat("  Min  Upsilon           :", round(min(upsilon_vec), 4), "\n")
cat("  Max  Upsilon           :", round(max(upsilon_vec), 4), "\n\n")

cat("--- Residuals ---\n")
cat("  Mean residual          :", round(mean(all_resid), 4), "\n")
cat("  SD   residual          :", round(sd(all_resid), 4), "\n")
cat("  Min  residual          :", round(min(all_resid), 4), "\n")
cat("  Max  residual          :", round(max(all_resid), 4), "\n\n")

cat("--- Predicted nu (test set, all CV entries) ---\n")
cat("  Mean nu_hat            :", round(mean(all_nu_hat), 4), "\n")
cat("  Median nu_hat          :", round(median(all_nu_hat), 4), "\n")
cat("  SD   nu_hat            :", round(sd(all_nu_hat), 4), "\n")
cat("  Range nu_hat           : [", round(min(all_nu_hat), 4), ",",
    round(max(all_nu_hat), 4), "]\n\n")

cat("--- Top selected variables for eta_mu (>= 10% of reps) ---\n")
avg_mu  <- rowMeans(coef_mu_mat)
freq_mu <- rowMeans(coef_mu_mat != 0)
top_mu  <- sort(freq_mu[freq_mu >= 0.10], decreasing = TRUE)
top_mu_names <- map_names(names(top_mu), var_dict)
for (i in seq_along(top_mu)) {
  avg_i <- avg_mu[names(top_mu)[i]]
  sd_i  <- sd(coef_mu_mat[names(top_mu)[i], ])
  ci_lo <- avg_i - 1.96 * sd_i
  ci_hi <- avg_i + 1.96 * sd_i
  cat(sprintf("  %-40s  freq = %.2f  avg = %+.3f  95%% CI [%+.3f, %+.3f]\n",
              top_mu_names[i], top_mu[i], avg_i, ci_lo, ci_hi))
}

cat("\n--- Top selected variables for eta_sigma (>= 10% of reps) ---\n")
avg_sigma  <- rowMeans(coef_sigma_mat)
freq_sigma <- rowMeans(coef_sigma_mat != 0)
top_sigma  <- sort(freq_sigma[freq_sigma >= 0.10], decreasing = TRUE)
top_sigma_names <- map_names(names(top_sigma), var_dict)
for (i in seq_along(top_sigma)) {
  avg_i <- avg_sigma[names(top_sigma)[i]]
  sd_i  <- sd(coef_sigma_mat[names(top_sigma)[i], ])
  ci_lo <- avg_i - 1.96 * sd_i
  ci_hi <- avg_i + 1.96 * sd_i
  cat(sprintf("  %-40s  freq = %.2f  avg = %+.3f  95%% CI [%+.3f, %+.3f]\n",
              top_sigma_names[i], top_sigma[i], avg_i, ci_lo, ci_hi))
}

cat("\n--- Top selected variables for eta_nu (>= 10% of reps) ---\n")
avg_nu  <- rowMeans(coef_nu_mat)
freq_nu <- rowMeans(coef_nu_mat != 0)
top_nu  <- sort(freq_nu[freq_nu >= 0.10], decreasing = TRUE)
top_nu_names <- map_names(names(top_nu), var_dict)
for (i in seq_along(top_nu)) {
  avg_i <- avg_nu[names(top_nu)[i]]
  sd_i  <- sd(coef_nu_mat[names(top_nu)[i], ])
  ci_lo <- avg_i - 1.96 * sd_i
  ci_hi <- avg_i + 1.96 * sd_i
  cat(sprintf("  %-40s  freq = %.2f  avg = %+.3f  95%% CI [%+.3f, %+.3f]\n",
              top_nu_names[i], top_nu[i], avg_i, ci_lo, ci_hi))
}

cat("\n================================================================\n")
cat("Plots saved to:", base_path, "\n")
cat("================================================================\n")
