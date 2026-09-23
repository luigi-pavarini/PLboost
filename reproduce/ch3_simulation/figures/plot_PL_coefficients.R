# Regenerate ONLY the coefficient boxplots (mu/sigma/lambda) at a WIDE aspect
# (2000x1450) so all six families fit on a single condensed page in the thesis.
# Lollipop PNGs are left untouched (they stay in the 3+3 layout).

save_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must hold PL_simulations/ and IPL_simulations/ with the .RData
families  <- c("NO", "TF", "LO", "PE", "Hyp", "SLASH")
BOX_W <- 2000; BOX_H <- 2000   # square, matching the rest of the document

for (fam in families) {
  rdata <- file.path(save_path, "PL_simulations",
                     paste0("sim_PL_M100_mstop100_", fam, ".RData"))
  if (!file.exists(rdata)) { cat("SKIP (missing):", rdata, "\n"); next }
  load(rdata)  # -> sim
  img_dir <- file.path(save_path, "PL_simulations", "imgs", fam)
  dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)
  save_box <- function(name, expr) {
    png(file.path(img_dir, paste0(name, ".png")), width = BOX_W, height = BOX_H, res = 330)
    on.exit(dev.off()); expr
  }

  save_box("mu_coefs", {
    par(mar = c(2.6, 4.0, 0.6, 0.6))
    boxplot(sim$est_mu[2,], sim$est_mu[3,], sim$est_mu[4,], sim$est_mu[5,],
            c(sim$est_mu[c(6:1001),]), ylim = c(-1.2, 0.5), xaxt = "n", las = 1,
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

  save_box("sigma_coefs", {
    par(mar = c(2.6, 4.0, 0.6, 0.6))
    boxplot(sim$est_sigma[4,], sim$est_sigma[5,], sim$est_sigma[6,], sim$est_sigma[7,],
            c(sim$est_sigma[c(2:3, 8:1001),]), ylim = c(-0.70, 1.00), xaxt = "n", las = 1,
            ylab = expression("Coefficients for" ~ eta[sigma]), main = "")
    abline(h = 0, lty = 2, col = "gray")
    axis(1, 1, expression(X[3])); axis(1, 2, expression(X[4]))
    axis(1, 3, expression(X[5])); axis(1, 4, expression(X[6])); axis(1, 5, "non-inf.")
    segments(0.5, -0.41, 1.5, -0.41, lty = 2, col = "red4")
    segments(1.5, -0.24, 2.5, -0.24, lty = 2, col = "red4")
    segments(2.5,  0.81, 3.5,  0.81, lty = 2, col = "red4")
    segments(3.5,  0.70, 4.5,  0.70, lty = 2, col = "red4")
    segments(4.5,  0,    5.5,  0,    lty = 2, col = "red4")
  })

  save_box("lambda_coefs", {
    par(mar = c(2.6, 4.0, 0.6, 0.6))
    boxplot(exp(sim$est_lambda[1,]), ylim = c(0, 3), xaxt = "n", las = 1,
            ylab = expression("Coefficients for" ~ lambda), main = "")
    abline(h = 0, lty = 2, col = "gray")
    axis(1, 1, "exp(intercept)")
    segments(0.5, 1.2, 1.5, 1.2, lty = 2, col = "red4")
  })

  cat("wide boxplots for", fam, "done\n")
}
cat("DONE.\n")
