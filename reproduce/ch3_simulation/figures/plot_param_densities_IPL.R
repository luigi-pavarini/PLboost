# Boosted IPL regression — Parameter density plots (single scenario: low)

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

# ---- Selectors --------------------------------------------------------------
fam      <- "NO"
scenario <- "low"

# ---- Image export settings --------------------------------------------------
img_width  <- 2000
img_height <- 2000
img_res    <- 330

save_path <- file.path(base_path, "IPL_simulations", "0_param_densities")
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# ---- Simulation constants (must match luigi_experimentations_IPL_n.R) -------
n     <- 800
p     <- 1000   # V1...V1000; alpha uses V999 and V1000
zeta  <- 2
c_par <- 0

# ============================================================================
# HELPER FUNCTION
# ============================================================================

# Boundary-corrected kernel density for data on [lo, hi].
# Uses the reflection method so the kernel does not leak mass outside support.
# Falls back to standard density() when lo/hi = NULL.
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
  
  if (!bounded) return(density(x, bw = bw, n = n_grid))
  
  if (!is.null(lo)) x <- c(x, 2 * lo - x)
  if (!is.null(hi)) x <- c(x, 2 * hi - x)
  
  lo_ext <- if (is.null(lo)) min(x) else lo
  hi_ext <- if (is.null(hi)) max(x) else hi
  
  d_raw   <- density(x, bw = bw, from = lo_ext, to = hi_ext, n = n_grid)
  d_raw$y <- d_raw$y * 2   # reflection doubles mass
  d_raw
}

# ============================================================================
# GENERATE TRUE PARAMETERS  (single scenario: low)
# ============================================================================

set.seed(42)
X_mat <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
colnames(X_mat) <- paste0("V", 1:p)

beta_mu_true    <- c(-0.99, 0.49, -1.17, 0.38, rep(0, p - 4))
beta_sigma_true <- c(0, 0, -0.41, -0.24, 0.81, rep(0, p - 5))
lambda_true     <- 1.2

# Low scenario: alpha driven by V999 and V1000
alpha_intercept <- -1.5
beta_alpha_true <- c(0.70, -0.50)

eta_mu_true    <- 0.5 + X_mat %*% beta_mu_true
eta_sigma_true <- 0.0 + X_mat %*% beta_sigma_true
eta_alpha_true <- alpha_intercept + X_mat[, c("V999", "V1000")] %*% beta_alpha_true

mu_true    <- pmax(pmin(plogis(eta_mu_true),    1 - 1e-6), 1e-6)
sigma_true <- exp(eta_sigma_true)
alpha_true <- pmax(pmin(plogis(eta_alpha_true), 1 - 1e-6), 1e-6)

set.seed(101)
y_true <- rIPL(
  n      = n,
  mu     = as.vector(mu_true),
  sigma  = as.vector(sigma_true),
  alpha  = as.vector(alpha_true),
  lambda = lambda_true,
  family = fam,
  zeta   = zeta,
  c      = c_par
)
y_true <- pmax(pmin(as.vector(y_true), 1), 0)

mu_vec    <- as.vector(mu_true)
sigma_vec <- as.vector(sigma_true)
alpha_vec <- as.vector(alpha_true)

message("alpha range = [", round(min(alpha_vec), 3), ", ", round(max(alpha_vec), 3), "]",
        " | mean = ", round(mean(alpha_vec), 3))

# Pre-compute densities
d_mu    <- bounded_density(mu_vec,    lo = 0,    hi = 1)
d_sigma <- bounded_density(sigma_vec, lo = NULL, hi = NULL)
d_alpha <- bounded_density(alpha_vec, lo = 0,    hi = 1)

# ============================================================================
# PLOT 1 — Y: histogram + density curve + point mass at 0
# ============================================================================

prop_zero <- mean(y_true == 0)
y_cont    <- y_true[y_true > 0]

# Density on open interval (0, 1) — no boundary reflection at 0 or 1
d_cont    <- bounded_density(y_cont, lo = NULL, hi = NULL)
keep      <- d_cont$x > 0 & d_cont$x < 1
dens_x    <- d_cont$x[keep]

# ---- Summary of Y (printed to console — use for table in text) --------------
cat("\n--- Summary of Y ---\n")
cat(sprintf("  Point mass at 0  : %.1f%%  (%d out of %d observations)\n",
            100 * prop_zero, sum(y_true == 0), n))
cat(sprintf("  Continuous part  : %.1f%%  (%d observations, y > 0)\n",
            100 * (1 - prop_zero), sum(y_true > 0), n))
cat("\n  Summary of continuous component (y > 0):\n")
s <- summary(y_cont)
cat(sprintf("    Min.   : %.4f\n", s["Min."]))
cat(sprintf("    1st Qu.: %.4f\n", s["1st Qu."]))
cat(sprintf("    Median : %.4f\n", s["Median"]))
cat(sprintf("    Mean   : %.4f\n", s["Mean"]))
cat(sprintf("    3rd Qu.: %.4f\n", s["3rd Qu."]))
cat(sprintf("    Max.   : %.4f\n", s["Max."]))
cat(sprintf("    SD     : %.4f\n", sd(y_cont)))
cat("--------------------\n\n")

png(file.path(save_path, paste0("param_density_Y_", fam, ".png")),
    width = img_width, height = img_height, res = img_res)

old_par <- par(mar = c(4.5, 4.5, 3, 2), family = "serif")

prop_one  <- mean(y_true == 1)
y_mid     <- y_true[y_true > 0 & y_true < 1]

# Point-mass counts
n_zero <- sum(y_true == 0)
n_one  <- sum(y_true == 1)

# Histogram of interior (0, 1) in frequency
h        <- hist(y_mid, breaks = 10, plot = FALSE)
ylim_top <- max(h$counts, n_zero, n_one, na.rm = TRUE) * 1.15

hist(y_mid, breaks = 10,
     freq = TRUE,
     col  = "white", border = "black",
     xlim = c(0, 1), ylim = c(0, ylim_top),
     xlab = "y", ylab = "Frequência",
     main = "",#paste0("n = ", n),
     las = 1)

# Point-mass stems and dots at 0 and 1
segments(x0 = 0, x1 = 0, y0 = 0, y1 = n_zero, lwd = 2, col = "black")
points(x = 0, y = n_zero, pch = 19, cex = 1.6, col = "black")

par(old_par)
dev.off()

# ============================================================================
# PLOT 2 — mu: histogram + density curve  (bounded [0, 1])
# ============================================================================

png(file.path(save_path, paste0("param_density_mu_", fam, ".png")),
    width = img_width, height = img_height, res = img_res)

par(family = "serif")
hist(mu_vec, breaks = 20, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(mu), ylab = "Density", las = 1)
lines(d_mu$x, d_mu$y, col = "gray30", lwd = 2)

dev.off()

# ============================================================================
# PLOT 3 — sigma: histogram + density curve  (unbounded > 0)
# ============================================================================

png(file.path(save_path, paste0("param_density_sigma_", fam, ".png")),
    width = img_width, height = img_height, res = img_res)

par(family = "serif")
hist(sigma_vec, breaks = 25, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(sigma), ylab = "Density", las = 1)
lines(d_sigma$x, d_sigma$y, col = "gray30", lwd = 2)

dev.off()

# ============================================================================
# PLOT 4 — alpha: histogram + density curve  (bounded [0, 1])
# ============================================================================

png(file.path(save_path, paste0("param_density_alpha_", fam, ".png")),
    width = img_width, height = img_height, res = img_res)

par(family = "serif")
hist(alpha_vec, breaks = 20, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(alpha), ylab = "Density", las = 1)
lines(d_alpha$x, d_alpha$y, col = "gray30", lwd = 2)

dev.off()

# ============================================================================
# BONUS — Combined 2x2 panel
# ============================================================================

png(file.path(save_path, paste0("param_density_ALL_", fam, ".png")),
    width = img_width * 1.6, height = img_height * 1.6, res = img_res)

old_par <- par(mfrow = c(2, 2), mar = c(4.5, 4.2, 2, 1.5), family = "serif")

# Panel 1: Y
ylim_top_y <- max(h$counts, n_zero, n_one, na.rm = TRUE) * 1.15

hist(y_mid, breaks = 10,
     freq = TRUE, col = "white", border = "black",
     xlim = c(0, 1), ylim = c(0, ylim_top_y),
     xlab = "y", ylab = "Frequência",
     main = paste0("n = ", n), las = 1)
segments(x0 = 0, x1 = 0, y0 = 0, y1 = n_zero, lwd = 2, col = "black")
points(x = 0, y = n_zero, pch = 19, cex = 1.6, col = "black")
# Panel 2: mu
hist(mu_vec, breaks = 20, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(mu), ylab = "Density", las = 1)
lines(d_mu$x, d_mu$y, col = "gray30", lwd = 2)

# Panel 3: sigma
hist(sigma_vec, breaks = 25, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(sigma), ylab = "Density", las = 1)
lines(d_sigma$x, d_sigma$y, col = "gray30", lwd = 2)

# Panel 4: alpha
hist(alpha_vec, breaks = 20, freq = FALSE, col = "gray85", border = "white",
     main = "", xlab = expression(alpha), ylab = "Density", las = 1)
lines(d_alpha$x, d_alpha$y, col = "gray30", lwd = 2)

par(old_par)
dev.off()

message("All IPL parameter-density plots saved to: ", save_path)