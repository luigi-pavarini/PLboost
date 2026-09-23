# Boosted PL regression — Parameter density plots across families
# ------------------------------------------------------------------------------
# LOGIC: One fixed X is generated once (set.seed(42)) and the TRUE coefficients
# (beta_mu, beta_sigma, lambda) are applied to compute the true mu and sigma.
# These are IDENTICAL across all families — only Y changes, because the family
# controls the shape of the error distribution. This isolates the effect of the
# family assumption on the response, holding everything else fixed.
# ------------------------------------------------------------------------------

rm(list = ls())

# ============================================================================
# SETUP
# ============================================================================
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


base_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # output root
library(PLboost)


require("PLreg")

# ---- Image export settings --------------------------------------------------
img_width  <- 2000
img_height <- 2000
img_res    <- 330

save_path  <- file.path(base_path, "PL_simulations", "0_param_densities")
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# ---- Simulation constants (must match luigi_experimentations_PL.r) ----------
n      <- 800
p      <- 1000
lambda <- 1.2
zeta   <- 2

beta_mu    <- c(-0.99,  0.49, -1.17,  0.38, rep(0, p - 4))
beta_sigma <- c(0, 0, -0.41, -0.24,  0.81,  0.70, rep(0, p - 6))

# True intercepts (as used in the DGP)
intercept_mu    <- 0.5
intercept_sigma <- 0.0

# ---- Families to overlay ----------------------------------------------------
families <- c("NO", "TF", "LO", "SN", "PE", "Hyp", "SLASH")

# ---- Line aesthetics — one entry per family ---------------------------------
fam_lty <- c(
  NO    = 1,
  TF    = 1,
  LO    = 2,
  PE    = 3,
  SN    = 3,
  Hyp   = 2,
  SLASH = 1
)

fam_col <- c(
  NO    = "#9ecae1",
  TF    = "#fc9272",
  LO    = "#fb6a4a",
  PE    = "#de2d26",
  SN    = "#dadaeb",
  Hyp   = "#9e9ac8",
  SLASH = "#6a51a3"
)

fam_lwd <- 1.8

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Boundary-corrected kernel density for data on [lo, hi].
# Uses the reflection method: mirrors data at both boundaries so that the
# kernel does not leak mass outside the support (avoids floating tails).
# Falls back to standard density() for unbounded parameters (lo/hi = NULL).
bounded_density <- function(x, lo = NULL, hi = NULL,
                            bw = "SJ", n_grid = 512) {
  
  x <- x[is.finite(x)]
  if (length(unique(x)) < 2) {
    center <- unique(x)[1]
    grid   <- seq(
      ifelse(is.null(lo), center - 1, lo),
      ifelse(is.null(hi), center + 1, hi),
      length.out = n_grid
    )
    return(list(x = grid, y = rep(0, length(grid))))
  }
  
  bounded <- !is.null(lo) || !is.null(hi)
  
  if (!bounded) {
    return(density(x, bw = bw, n = n_grid))
  }
  
  if (!is.null(lo)) x <- c(x, 2 * lo - x)
  if (!is.null(hi)) x <- c(x, 2 * hi - x)
  
  lo_ext <- if (is.null(lo)) min(x) else lo
  hi_ext <- if (is.null(hi)) max(x) else hi
  
  d_raw <- density(x, bw = bw, from = lo_ext, to = hi_ext, n = n_grid)
  d_raw$y <- d_raw$y * 2
  d_raw
}


# Overlay density lines for all families on a single plot
density_overlay <- function(samp_list, param,
                            xlab_expr,
                            lo = NULL, hi = NULL,
                            xlim = NULL,
                            fam_lty, fam_col, fam_lwd,
                            legend_pos = "topright") {
  
  all_vals <- unlist(lapply(samp_list, `[[`, param))
  all_vals <- all_vals[is.finite(all_vals)]
  
  if (is.null(xlim)) {
    xlim <- range(all_vals)
    rng  <- diff(xlim)
    if (!is.finite(rng) || rng <= 0) rng <- 1
    xlim <- xlim + c(-0.04, 0.04) * rng
  }
  
  if (!is.null(lo)) xlim[1] <- max(xlim[1], lo)
  if (!is.null(hi)) xlim[2] <- min(xlim[2], hi)
  
  dens_list <- lapply(samp_list, function(s) {
    bounded_density(s[[param]], lo = lo, hi = hi)
  })
  
  dens_list <- lapply(dens_list, function(d) {
    keep <- d$x >= xlim[1] & d$x <= xlim[2]
    list(x = d$x[keep], y = d$y[keep])
  })
  
  ylim <- c(0, max(vapply(dens_list, function(d) max(d$y, na.rm = TRUE),
                          numeric(1))) * 1.12)
  if (!is.finite(ylim[2]) || ylim[2] <= 0) ylim[2] <- 1
  
  fam_names <- names(samp_list)
  
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = xlab_expr, ylab = "Density",
       main = "", las = 1)
  
  for (fam in fam_names) {
    d <- dens_list[[fam]]
    lines(d$x, d$y,
          lty = fam_lty[fam],
          col = fam_col[fam],
          lwd = fam_lwd)
  }
  
  legend(legend_pos,
         legend = fam_names,
         lty    = fam_lty[fam_names],
         col    = fam_col[fam_names],
         lwd    = fam_lwd,
         bty    = "n",
         cex    = 0.85)
}

# ============================================================================
# GENERATE ONE FIXED X — TRUE mu AND sigma ARE THE SAME FOR ALL FAMILIES
# Only Y differs across families (family controls the error shape, not mu/sigma)
# ============================================================================

set.seed(42)
X_mat <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)

# True linear predictors using true coefficients
eta_mu_true    <- intercept_mu    + X_mat %*% beta_mu
eta_sigma_true <- intercept_sigma + X_mat %*% beta_sigma

# True mu and sigma — identical for all families
mu_true    <- plogis(eta_mu_true)
mu_true    <- pmax(pmin(as.vector(mu_true),    1 - 1e-6), 1e-6)
sigma_true <- exp(eta_sigma_true)
sigma_true <- as.vector(sigma_true)

samp <- list()

for (fam in families) {
  
  # Generate Y under this family, using the same true mu, sigma, lambda
  y_fam <- PLreg::rPL(n      = n,
                      mu     = mu_true,
                      sigma  = sigma_true,
                      lambda = lambda,
                      family = fam,
                      zeta   = zeta)
  y_fam <- pmax(pmin(as.vector(y_fam), 1 - 1e-6), 1e-6)
  
  # mu and sigma are the same for every family — stored for reference plots
  samp[[fam]] <- list(
    Y     = y_fam,
    mu    = mu_true,
    sigma = sigma_true
  )
  
  message("Generated Y for family: ", fam)
}

# Pre-compute mu and sigma densities once (they are family-invariant)
d_mu    <- bounded_density(mu_true,    lo = 0,    hi = 1)
d_sigma <- bounded_density(sigma_true, lo = NULL, hi = NULL)

# ============================================================================
# PLOT 1 — Density of Y by family  (bounded [0, 1])
# This is the key plot: shows how the family assumption shapes the response
# ============================================================================

png(file.path(save_path, "param_density_Y.png"),
    width = img_width, height = img_height, res = img_res)

density_overlay(samp, param = "Y",
                xlab_expr  = expression(Y),
                lo = 0, hi = 1,
                fam_lty = fam_lty, fam_col = fam_col, fam_lwd = fam_lwd,
                legend_pos = "topright")

dev.off()

# ============================================================================
# PLOT 2 — Density of true mu  (same for all families — single curve)
# ============================================================================

png(file.path(save_path, "param_density_mu.png"),
    width = img_width, height = img_height, res = img_res)

hist(mu_true, breaks = 15, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(mu), ylab = "Density", las = 1)
lines(d_mu$x, d_mu$y, col = "gray30", lwd = 2)

dev.off()

# ============================================================================
# PLOT 3 — Density of true sigma  (same for all families — single curve)
# ============================================================================

png(file.path(save_path, "param_density_sigma.png"),
    width = img_width, height = img_height, res = img_res)
hist(sigma_true, breaks = 15, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(sigma), ylab = "Density", las = 1)
lines(d_sigma$x, d_sigma$y, col = "gray30", lwd = 2)

dev.off()


# ============================================================================
# COMBINED — Y by family | true mu | true sigma
# ============================================================================

png(file.path(save_path, "param_density_ALL.png"),
    width = img_width * 1.5, height = img_height, res = img_res)

old_par <- par(mfrow = c(1, 3), mar = c(4.5, 4, 2, 1), oma = c(0, 0, 2, 0))

# Panel 1: Y — differs by family
density_overlay(samp, param = "Y",
                xlab_expr  = expression(Y),
                lo = 0, hi = 1,
                fam_lty = fam_lty, fam_col = fam_col, fam_lwd = fam_lwd,
                legend_pos = "topright")
title(sub = expression(Y), line = 0.3, cex.sub = 0.9)

# Panel 2: mu — single curve (same for all families)
hist(mu_true, breaks = 20, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(mu), ylab = "Density", las = 1)
lines(d_mu$x, d_mu$y, col = "gray30", lwd = 2)
title(sub = expression(mu), line = 0.3, cex.sub = 0.9)

# Panel 3: sigma — histogram + density (same for all families)
hist(sigma_true, breaks = 20, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(sigma), ylab = "Density", las = 1)
lines(d_sigma$x, d_sigma$y, col = "gray30", lwd = 2)
title(sub = expression(sigma), line = 0.3, cex.sub = 0.9)

mtext("True parameter densities and Y by family", outer = TRUE, cex = 1.1)

par(old_par)
dev.off()

message("All density plots saved to: ", save_path)