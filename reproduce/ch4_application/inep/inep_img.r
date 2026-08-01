# ==============================================================================
# inep_img_pl_cv.r
# Visualization for the INEP MEDICINA/SP Conceito Enade (Continuo) application.
# Power-Logit (PL) CV version. STYLE MATCHES credit_img_ipl_cv_10rep.r exactly
# (2000x2000, rel_dim 1.8, gray85/red4 palette, per-family imgs/<fam> folder,
# same plot names) so the two applications are visually consistent in the paper /
# dissertation. PL differences vs the IPL script: fields are error_pl / mae_pl,
# lambda ~ 1 (scalar shape), and there is NO inflation parameter -> the alpha
# plots (plot_alpha_dist, plot_alpha_vs_mu, coef/sel alpha) are dropped.
# Set `fam` to switch family; winner = "PE".
# ==============================================================================

suppressMessages({library(data.table)})

fam          <- "NO"   # ← NO, TF, LO, PE, SN, SLASH, Hyp  (featured family = NO: lowest MSE, beats Beta in every run, simplest generator)
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

img_width  <- 2000
img_height <- 2000
img_res    <- 330
rel_dim    <- 1.8

source(file.path(app_path, "essentials", "cv_load.R"))

prefix   <- "sim_results_medicina_SP"
loaded   <- load_cv_results(results_path, prefix, fam)
sim_results <- loaded$sim_results
sim_results <- sim_results[!vapply(sim_results, is.null, logical(1))]
var_dict <- loaded$var_dict
n_success <- length(sim_results)
cat("Family =", fam, "| successful CV entries:", n_success, "\n")

# ==============================================================================
# HELPERS  (identical to the IPL script; alpha removed)
# ==============================================================================
map_names <- function(vnames, dict) {
  idx <- match(vnames, dict$new_name)
  ifelse(is.na(idx), vnames, dict$orig_name[idx])
}

all_test_df <- do.call(rbind, lapply(sim_results, function(x)
  data.frame(idx = x$test_idx, y = x$y_test)))
all_test_df <- all_test_df[!duplicated(all_test_df$idx), ]
y_all       <- all_test_df$y

errors_pl     <- sapply(sim_results, function(x) x$error_pl) * 100
mae_pl        <- sapply(sim_results, function(x) x$mae_pl)
upsilon_vec   <- sapply(sim_results, function(x) x$upsilon)
zeta_best_vec <- sapply(sim_results, function(x) x$zeta_best)
all_ytest     <- unlist(lapply(sim_results, function(x) x$y_test))
all_muhat     <- unlist(lapply(sim_results, function(x) x$mu_hat))
all_resid     <- all_ytest - all_muhat
lambda_hat_vec <- sapply(sim_results, function(x) exp(x$coef_lambda["(Intercept)"] + attr(x$coef_lambda, "offset")))

has_zeta <- !all(sapply(zeta_best_vec, is.null))

build_coef_mat <- function(coef_list, exclude = "(Intercept)") {
  all_vars <- unique(unlist(lapply(coef_list, function(cv) names(cv)[names(cv) != exclude])))
  mat <- matrix(0, nrow = length(all_vars), ncol = length(coef_list),
                dimnames = list(all_vars, NULL))
  for (r in seq_along(coef_list)) {
    cv <- coef_list[[r]]; nms <- names(cv)[names(cv) != exclude]; mat[nms, r] <- cv[nms]
  }
  mat
}
coef_mu_mat    <- build_coef_mat(lapply(sim_results, function(x) x$coef_mu))
coef_sigma_mat <- build_coef_mat(lapply(sim_results, function(x) x$coef_sigma))

# ==============================================================================
# PLOT 1: Response distribution (tx open on (0,1) — no boundary mass)
# ==============================================================================
png(file.path(base_path, "plot_response_dist.png"),
    width = img_width, height = img_height, res = img_res)
old_par <- par(mar = c(4.5, 4.5, 2, 2))
hist(y_all, breaks = seq(0, 1, length.out = 5), freq = FALSE,  # exactly 4 bars
     col = "gray85", border = "white",
     xlab = "Conceito Enade / 5", ylab = "Density", main = "", las = 1,
     xlim = c(0, 1))
lines(density(y_all, bw = "SJ"), col = "gray30", lwd = 2)
abline(v = mean(y_all), lty = 3, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c("Kernel density", "Mean"),
       col = c("gray30", "gray40"), lty = c(1, 3), lwd = c(2, 1.5),
       bty = "n", cex = 0.85)
par(old_par); dev.off()

# ==============================================================================
# PLOT 2 / 2b / 2c: MSE, MAE, Upsilon distributions
# ==============================================================================
one_box <- function(vec, ylab_expr, fname, dig) {
  png(file.path(base_path, fname), width = img_width, height = img_height, res = img_res)
  op <- par(mar = c(4.5, 4.5, 2, 2))
  boxplot(vec, ylab = ylab_expr, xaxt = "n", col = "gray85", border = "gray30",
          las = 1, outline = TRUE, main = "")
  axis(1, at = 1, labels = "PL")
  abline(h = mean(vec), lty = 2, col = "red4", lwd = 1.5)
  legend("topright", legend = paste0("Mean = ", round(mean(vec), dig),
         "   SD = ", round(sd(vec), dig)), bty = "n", cex = 0.85)
  par(op); dev.off()
}
one_box(errors_pl,   expression("MSE" %*% 100), "plot_mse_dist.png",     3)
one_box(mae_pl,      "MAE",                      "plot_mae_dist.png",     4)
one_box(upsilon_vec, expression(Upsilon),        "plot_upsilon_dist.png", 4)

# ==============================================================================
# PLOT 2d: Zeta selection
# ==============================================================================
if (has_zeta) {
  png(file.path(base_path, "plot_zeta_selection.png"),
      width = img_width, height = img_height, res = img_res)
  op <- par(mar = c(4.5, 4.5, 2, 2))
  zeta_tab <- table(unlist(zeta_best_vec))
  barplot(zeta_tab, xlab = expression(zeta), ylab = "Number of CV entries",
          col = "gray85", border = "gray30", las = 1, main = "")
  abline(h = 0, col = "gray60")
  legend("topright", legend = paste0("Most selected: zeta = ",
         names(zeta_tab)[which.max(zeta_tab)], "  (", max(zeta_tab), "/", n_success, ")"),
         bty = "n", cex = 0.85)
  par(op); dev.off()
}

# ==============================================================================
# PLOT 3: Observed vs predicted — overlapping densities
# ==============================================================================
png(file.path(base_path, "plot_pred_vs_obs.png"),
    width = img_width, height = img_height, res = img_res)
op <- par(mar = c(4.5, 4.5, 2, 2))
d_obs  <- density(all_ytest, bw = 0.02); d_pred <- density(all_muhat, bw = 0.02)
xlim_r <- range(c(d_obs$x, d_pred$x)); ylim_r <- c(0, max(c(d_obs$y, d_pred$y)) * 1.1)
plot(d_obs, col = "gray30", lwd = 2, xlim = xlim_r, ylim = ylim_r,
     xlab = "Conceito Enade / 5", ylab = "Density", main = "", las = 1)
polygon(c(d_obs$x, rev(d_obs$x)), c(d_obs$y, rep(0, length(d_obs$y))),
        col = adjustcolor("gray30", 0.15), border = NA)
lines(d_pred, col = "red4", lwd = 2, lty = 2)
polygon(c(d_pred$x, rev(d_pred$x)), c(d_pred$y, rep(0, length(d_pred$y))),
        col = adjustcolor("red4", 0.15), border = NA)
abline(v = mean(all_ytest), lty = 3, col = "gray30", lwd = 1.2)
abline(v = mean(all_muhat), lty = 3, col = "red4",   lwd = 1.2)
legend("topleft", legend = c(expression(y[test]), expression(hat(mu))),
       col = c("gray30", "red4"), lty = c(1, 2), lwd = c(2, 2),
       fill = c(adjustcolor("gray30", 0.15), adjustcolor("red4", 0.15)),
       border = NA, bty = "n", cex = 0.9)
par(op); dev.off()

# ==============================================================================
# PLOT 4a/4b: Residuals — density + boxplot
# ==============================================================================
png(file.path(base_path, "plot_residuals_density.png"),
    width = img_width, height = img_height, res = img_res)
op <- par(mar = c(4.5, 4.5, 2, 2))
hist(all_resid, breaks = 8, freq = FALSE, col = "gray85", border = "white",
     xlab = expression(y[test] - hat(mu)), ylab = "Density", main = "", las = 1)
lines(density(all_resid), col = "gray30", lwd = 2)
abline(v = 0, col = "red4", lty = 2, lwd = 1.5)
par(op); dev.off()

png(file.path(base_path, "plot_residuals_boxplot.png"),
    width = img_width, height = img_height, res = img_res)
op <- par(mar = c(4.5, 4.5, 2, 2))
boxplot(all_resid, ylab = expression(y[test] - hat(mu)), xaxt = "n",
        col = "gray85", border = "gray30", outline = TRUE, las = 1, main = "")
abline(h = 0, col = "red4", lty = 2, lwd = 1.5)
axis(1, at = 1, labels = "Residuals")
par(op); dev.off()

# ==============================================================================
# PLOT 5: Lambda — boxplot  (lambda ~ 1 scalar shape)
# ==============================================================================
png(file.path(base_path, "plot_lambda_boxplot.png"),
    width = img_width, height = img_height, res = img_res)
op <- par(mar = c(4.5, 4.5, 2, 2))
lam_r <- range(lambda_hat_vec); lam_pad <- 0.6 * diff(lam_r) + 0.05
boxplot(lambda_hat_vec, ylab = expression("Estimated  " * hat(lambda)), xaxt = "n",
        col = "gray85", border = "gray30", las = 1, outline = TRUE, main = "",
        ylim = c(lam_r[1] - lam_pad, lam_r[2] + lam_pad))
axis(1, at = 1, labels = "PL")
abline(h = mean(lambda_hat_vec), lty = 2, col = "red4", lwd = 1.5)
legend("topright", legend = paste0("Mean = ", round(mean(lambda_hat_vec), 3),
       "   SD = ", round(sd(lambda_hat_vec), 3)), bty = "n", cex = 0.85)
par(op); dev.off()

# ==============================================================================
# PLOTS 6-9: Coefficient forests + selection lollipops (mu, sigma)
#            (verbatim styling from the IPL forest/lollipop functions)
# ==============================================================================
plot_coef_forest <- function(coef_mat, var_dict, param_label = "mu",
                             min_freq = 0.10, filename, legend_pos = "bottomright") {
  avg <- rowMeans(coef_mat); sdev <- apply(coef_mat, 1, sd); freq <- rowMeans(coef_mat != 0)
  keep <- freq >= min_freq; avg <- avg[keep]; sdev <- sdev[keep]
  if (length(avg) == 0) return(invisible())
  top10 <- order(abs(avg), decreasing = TRUE)[seq_len(min(10, length(avg)))]
  avg <- avg[top10]; sdev <- sdev[top10]
  ord <- order(abs(avg)); avg <- avg[ord]; sdev <- sdev[ord]
  orig <- map_names(names(avg), var_dict)
  orig <- ifelse(nchar(orig) > 34, paste0(substr(orig, 1, 33), "…"), orig)
  n <- length(avg)
  ci68_lo <- avg - sdev;       ci68_hi <- avg + sdev
  ci90_lo <- avg - 1.645*sdev; ci90_hi <- avg + 1.645*sdev
  ci95_lo <- avg - 1.960*sdev; ci95_hi <- avg + 1.960*sdev
  xlab_expr <- if (param_label == "mu") expression("Average coefficient  " * eta[mu]) else
                                        expression("Average coefficient  " * eta[sigma])
  rect_half <- 0.30
  col_pos <- adjustcolor("#4682B4", 0.45); col_neg <- adjustcolor("#B44646", 0.45)
  col_pos_90 <- adjustcolor("#4682B4", 0.80); col_neg_90 <- adjustcolor("#B44646", 0.80)
  col_bord <- "#2a2a2a"; x_range <- range(c(ci95_lo, ci95_hi, 0)) * 1.15
  x_range[1] <- min(x_range[1], -0.01)
  png(filename, width = img_width * rel_dim, height = max(img_height, 55*n + 150), res = img_res)
  op <- par(mar = c(4.5, 11, 2, 3))
  plot(NA, xlim = x_range, ylim = c(0.5, n + 0.5), yaxt = "n",
       xlab = xlab_expr, ylab = "", main = "", las = 1)
  abline(v = 0, col = "gray60", lty = 2)
  rect_col <- ifelse(avg >= 0, col_pos, col_neg); wick90 <- ifelse(avg >= 0, col_pos_90, col_neg_90)
  cap95 <- rect_half * 0.75
  segments(ci95_lo, seq_len(n), ci95_hi, seq_len(n), col = col_bord, lwd = 1.4, lty = 3)
  segments(ci95_lo, seq_len(n) - cap95, ci95_lo, seq_len(n) + cap95, col = col_bord, lwd = 1.4, lty = 3)
  segments(ci95_hi, seq_len(n) - cap95, ci95_hi, seq_len(n) + cap95, col = col_bord, lwd = 1.4, lty = 3)
  segments(ci90_lo, seq_len(n), ci90_hi, seq_len(n), col = wick90, lwd = 1.4)
  rect(ci68_lo, seq_len(n) - rect_half, ci68_hi, seq_len(n) + rect_half,
       col = rect_col, border = col_bord, lwd = 0.6)
  points(avg, seq_len(n), pch = 22, col = col_bord, bg = "white", cex = 1.3, lwd = 1.2)
  axis(2, at = seq_len(n), labels = orig, las = 1, cex.axis = 0.78)
  par(op); dev.off()
}

lollipop_sel <- function(coef_mat, var_dict, top_n = 15, param_label = "mu", filename) {
  # Proportion of folds in which each variable received a non-zero coefficient.
  # (Previously this normalized each fold's selections to sum to 1, giving a share of
  #  the selection budget rather than a selection frequency, which contradicted the
  #  tables where e.g. prop_enem = 1.00.)
  M <- ncol(coef_mat)
  avg <- rowMeans(coef_mat != 0)
  sdv <- sqrt(avg * (1 - avg) / M)   # binomial standard error
  tn <- min(top_n, length(avg)); ord <- order(avg, decreasing = TRUE)[seq_len(tn)]
  ta <- avg[ord]; ts <- sdv[ord]; tnm <- map_names(names(ta), var_dict)
  tnm <- ifelse(nchar(tnm) > 34, paste0(substr(tnm, 1, 33), "…"), tnm)
  p <- c(ta, mean(avg[-ord])); s <- c(ts, sqrt(sum(sdv[-ord]^2))); lbl <- c(tnm, "others")
  np <- length(p)
  xlab_expr <- if (param_label == "mu") expression("Proportion of folds selected" ~ eta[mu]) else
                                        expression("Proportion of folds selected" ~ eta[sigma])
  png(filename, width = img_width * rel_dim, height = img_height, res = img_res)
  op <- par(mar = c(4, 12, 2, 2)); xm <- max(p + s, na.rm = TRUE) * 1.15
  plot(p, seq_len(np), ylim = c(np, 1), xlim = c(0, xm), yaxt = "n",
       xlab = xlab_expr, ylab = "", main = "", pch = 19, col = "gray30", cex = 1.3)
  segments(0, seq_len(np), p, seq_len(np), col = "gray80", lwd = 0.8, lty = 2)
  xlo <- pmax(p - s, 0); xhi <- p + s; cap <- 0.18
  segments(xlo, seq_len(np), xhi, seq_len(np), col = "#4682B4", lwd = 1.4)
  segments(xlo, seq_len(np) - cap, xlo, seq_len(np) + cap, col = "#4682B4", lwd = 1.4)
  segments(xhi, seq_len(np) - cap, xhi, seq_len(np) + cap, col = "#4682B4", lwd = 1.4)
  axis(2, at = seq_len(np), labels = lbl, las = 1, cex.axis = 0.80)
  par(op); dev.off()
}

plot_coef_forest(coef_mu_mat,    var_dict, "mu",    filename = file.path(base_path, "plot_coef_mu.png"), legend_pos = "bottomleft")
lollipop_sel(   coef_mu_mat,    var_dict, 15, "mu",    filename = file.path(base_path, "plot_sel_mu.png"))
plot_coef_forest(coef_sigma_mat, var_dict, "sigma", filename = file.path(base_path, "plot_coef_sigma.png"))
lollipop_sel(   coef_sigma_mat, var_dict, 15, "sigma", filename = file.path(base_path, "plot_sel_sigma.png"))

cat("\nSaved plots to:", base_path, "\n"); print(list.files(base_path))
