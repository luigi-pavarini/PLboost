# ==============================================================================
# testing_PL_diss.R
# Dissertation validation script: Boosting (PowerLogit) vs PLreg MLE
#
# For each scenario, fits both estimators and produces:
#   - Side-by-side coefficient comparison plot (PNG)
#   - Printed comparison table in the console
#
# Requires: gamboostLSS, mboost, faux, PLreg + your local source files
# ==============================================================================

library(PLboost)
library(gamboostLSS)
library(mboost)
library(faux)
library(PLreg)
base_path <- NULL  # not needed: the implementation now ships as the PLboost package


# output directory (same folder as this script, or change as needed)
out_dir <- Sys.getenv("PLBOOST_OUT", unset = ".")

# ==============================================================================
# SETTINGS
# ==============================================================================

set.seed(123)
n           <- 100
p           <- 5
fam         <- "NO"
zeta        <- 2
lambda_true <- 1.2
nu          <- 0.01
mstop       <- 500

img_width  <- 2000*1.8
img_height <- 2000
img_res    <- 330

family_pl <- PowerLogit(family = fam, stabilization = "MAD", zeta = zeta)

# ==============================================================================
# HELPERS
# ==============================================================================

# gamboostLSS can label the intercept differently — normalise to "(Intercept)"
normalise_intercept <- function(x) {
  names(x) <- gsub("^.*[Ii]ntercept.*$", "(Intercept)", names(x))
  x
}

# ==============================================================================
# HELPER: build comparison data.frame and plot
# ==============================================================================

make_comparison_plot <- function(boost_mu, boost_sigma, boost_lambda,
                                 mle_mu,   mle_sigma,   mle_lambda,
                                 true_mu,  true_sigma,  true_lambda,
                                 scenario_label, out_path) {
  
  # colours matching application plot palette (same lightness as #4682B4 / #B44646)
  col_boost      <- adjustcolor("#E07B39", alpha.f = 0.80)  # orange  — Boosting
  col_boost_fill <- adjustcolor("#E07B39", alpha.f = 0.45)
  col_mle        <- adjustcolor("#3DA85A", alpha.f = 0.80)  # green   — MLE
  col_mle_fill   <- adjustcolor("#3DA85A", alpha.f = 0.45)
  col_true       <- "red4"
  col_bord       <- "#2a2a2a"
  bar_half       <- 0.28   # half-height of each horizontal bar
  
  submodels <- list(
    list(name  = "mu",
         vars  = names(true_mu),
         true  = true_mu,  boost = boost_mu,  mle = mle_mu,
         xlab  = expression("Coefficient value  " * eta[mu])),
    list(name  = "sigma",
         vars  = names(true_sigma),
         true  = true_sigma, boost = boost_sigma, mle = mle_sigma,
         xlab  = expression("Coefficient value  " * eta[sigma])),
    list(name  = "lambda",
         vars  = "lambda",
         true  = c(lambda = unname(true_lambda)),
         boost = c(lambda = unname(boost_lambda)),
         mle   = c(lambda = unname(mle_lambda)),
         xlab  = expression(log(lambda)))
  )
  
  base <- sub("\\.png$", "", out_path)
  
  for (sm in submodels) {
    panel_path <- paste0(base, "_", sm$name, ".png")
    
    vars   <- sm$vars
    n_vars <- length(vars)
    
    true_vals  <- as.numeric(sm$true[vars])
    boost_vals <- as.numeric(sm$boost[vars])
    mle_vals   <- as.numeric(sm$mle[vars])
    
    # two horizontal rows per variable: Boosting (top), MLE (bottom)
    y_boost <- seq(n_vars, 1, by = -1) * 3 - 1   # top track
    y_mle   <- y_boost - 1                         # bottom track
    y_label <- (y_boost + y_mle) / 2              # centre for axis label
    
    all_vals <- c(true_vals, boost_vals, mle_vals, 0)
    xrange   <- range(all_vals, na.rm = TRUE) * 1.15
    xrange[1] <- min(xrange[1], -0.05)
    xrange[2] <- xrange[2] * 1.20
    
    ylim <- c(min(y_mle) - bar_half - 0.5, max(y_boost) + bar_half + 0.5)
    
    img_height_px <- max(img_height, 55 * n_vars * 2 + 150)
    
    png(panel_path, width = img_width, height = img_height_px, res = img_res)
    par(mar = c(4.5, 9, 2, 2))
    
    plot(NA, xlim = xrange, ylim = ylim,
         yaxt = "n", xlab = sm$xlab, ylab = "", main = "", las = 1,
         cex.lab = 0.95)
    
    # vertical zero reference line
    abline(v = 0, col = "gray60", lty = 2, lwd = 1)
    
    # true value: vertical dashed line per variable spanning both tracks
    for (k in seq_along(vars)) {
      segments(true_vals[k], y_mle[k]   - bar_half,
               true_vals[k], y_boost[k] + bar_half,
               col = col_true, lty = 2, lwd = 1.8)
    }
    
    # Boosting bars (filled rectangle from 0 to estimate)
    rect(xleft   = pmin(0, boost_vals),
         ybottom = y_boost - bar_half,
         xright  = pmax(0, boost_vals),
         ytop    = y_boost + bar_half,
         col     = col_boost_fill,
         border  = col_boost,
         lwd     = 0.8)
    
    # MLE bars
    rect(xleft   = pmin(0, mle_vals),
         ybottom = y_mle - bar_half,
         xright  = pmax(0, mle_vals),
         ytop    = y_mle + bar_half,
         col     = col_mle_fill,
         border  = col_mle,
         lwd     = 0.8)
    
    # mean dot on each bar (square marker, white fill)
    points(boost_vals, y_boost, pch = 22, col = col_bord,
           bg = "white", cex = 1.1, lwd = 1.0)
    points(mle_vals,   y_mle,   pch = 22, col = col_bord,
           bg = "white", cex = 1.1, lwd = 1.0)
    
    axis(2, at = y_label, labels = vars, las = 1, cex.axis = 0.80)
    
    legend("bottomright",
           legend = c("True value", "Boosting", "PLreg MLE"),
           fill   = c(NA,           col_boost_fill, col_mle_fill),
           border = c(NA,           col_boost,      col_mle),
           lty    = c(2,            NA,             NA),
           lwd    = c(1.8,          NA,             NA),
           col    = c(col_true,     NA,             NA),
           pch    = c(NA,           22,             22),
           pt.bg  = c(NA,           "white",        "white"),
           pt.cex = 1.1,
           bty    = "n",
           cex    = 0.80)
    
    dev.off()
    message("Plot saved: ", panel_path)
  }
}

# ==============================================================================
# HELPER: print console table
# ==============================================================================

print_table <- function(boost_mu, boost_sigma, boost_lambda,
                        mle_mu,   mle_sigma,   mle_lambda,
                        true_mu,  true_sigma,  true_lambda,
                        label) {
  
  cat("\n", strrep("=", 62), "\n")
  cat(" Coefficient comparison:", label, "\n")
  cat(strrep("=", 62), "\n")
  
  cat(sprintf("%-12s  %8s  %8s  %8s\n", "Parameter", "True", "Boosting", "PLreg"))
  cat(strrep("-", 44), "\n")
  
  # mu
  for (v in names(true_mu)) {
    bv <- if (v %in% names(boost_mu)) boost_mu[v] else NA
    mv <- if (v %in% names(mle_mu))   mle_mu[v]   else NA
    cat(sprintf("mu: %-8s  %8.4f  %8.4f  %8.4f\n", v, true_mu[v], bv, mv))
  }
  
  # sigma
  for (v in names(true_sigma)) {
    bv <- if (v %in% names(boost_sigma)) boost_sigma[v] else NA
    mv <- if (v %in% names(mle_sigma))   mle_sigma[v]   else NA
    cat(sprintf("sigma: %-5s  %8.4f  %8.4f  %8.4f\n", v, true_sigma[v], bv, mv))
  }
  
  # lambda
  cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f\n",
              "lambda", true_lambda, boost_lambda, mle_lambda))
  cat(strrep("=", 62), "\n\n")
}

# ==============================================================================
# SCENARIO 1
# ==============================================================================

cat("\n>>> Fitting Scenario 1...\n")

beta_mu_1    <- c(Intercept = 0.5,  V1 = -0.99, V2 =  0.49, V3 = -1.17, V4 =  0.38, V5 = 0)
beta_sigma_1 <- c(Intercept = 0.0,  V1 =  0.00, V2 =  0.00, V3 = -0.41, V4 = -0.24, V5 = 0.81)

X <- rnorm_multi(n, vars = p, mu = rep(0, p), sd = rep(1, p), r = 0.5,
                 varnames = paste0("V", 1:p), empirical = FALSE)
X_mat <- as.matrix(X)

eta_mu_1    <- beta_mu_1["Intercept"]    + X_mat %*% beta_mu_1[paste0("V", 1:p)]
eta_sigma_1 <- beta_sigma_1["Intercept"] + X_mat %*% beta_sigma_1[paste0("V", 1:p)]

mu_1    <- pmax(pmin(plogis(eta_mu_1),    1 - 1e-6), 1e-6)
sigma_1 <- exp(eta_sigma_1)

y1 <- PLreg::rPL(n = n, mu = mu_1, sigma = sigma_1,
                 lambda = lambda_true, family = fam, zeta = zeta)
y1 <- pmax(pmin(y1, 1 - 1e-6), 1e-6)

data_pl1 <- data.frame(y = y1, X)

# --- Boosting ---
model_boost1 <- glmboostLSS(
  formula  = list(mu     = y ~ V1 + V2 + V3 + V4,
                  sigma  = y ~ V3 + V4 + V5,
                  lambda = y ~ 1),
  data     = data_pl1,
  families = family_pl,
  control  = boost_control(mstop = mstop, nu = nu, trace = FALSE)
)

b1_mu     <- normalise_intercept(coef(model_boost1, parameter = "mu"))
b1_sigma  <- normalise_intercept(coef(model_boost1, parameter = "sigma"))
b1_lambda <- unlist(coef(model_boost1, parameter = "lambda"))  # already log scale

# --- PLreg MLE (lambda fixed at true value) ---
model_reg1 <- PLreg::PLreg(
  y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5,
  data   = data_pl1,
  family = fam,
  lambda = lambda_true,
  zeta   = zeta
)
all1      <- coef(model_reg1)
mle1_mu   <- all1[!grepl("^\\(sigma\\)", names(all1))]
mle1_sig  <- setNames(all1[grepl("^\\(sigma\\)", names(all1))],
                      gsub("^\\(sigma\\)_", "", names(all1)[grepl("^\\(sigma\\)", names(all1))]))

# true values including intercepts
true_mu_1    <- c("(Intercept)" =  0.50, V1 = -0.99, V2 =  0.49, V3 = -1.17, V4 =  0.38)
true_sigma_1 <- c("(Intercept)" =  0.00, V3 = -0.41, V4 = -0.24, V5 =  0.81)

print_table(
  boost_mu     = b1_mu,     boost_sigma  = b1_sigma,  boost_lambda = b1_lambda,
  mle_mu       = mle1_mu,   mle_sigma    = mle1_sig,  mle_lambda   = log(lambda_true),
  true_mu      = true_mu_1, true_sigma   = true_sigma_1, true_lambda = log(lambda_true),
  label        = "Scenario 1 (PL-NO, n=100, p=5)"
)

make_comparison_plot(
  boost_mu     = b1_mu,     boost_sigma  = b1_sigma,  boost_lambda = b1_lambda,
  mle_mu       = mle1_mu,   mle_sigma    = mle1_sig,  mle_lambda   = log(lambda_true),
  true_mu      = true_mu_1, true_sigma   = true_sigma_1, true_lambda = log(lambda_true),
  scenario_label = "Scenario 1 — PL-NO  (n = 100, p = 5)",
  out_path       = file.path(out_dir, "diss_PL_scenario1.png")
)

# ==============================================================================
# SCENARIO 2
# ==============================================================================

cat("\n>>> Fitting Scenario 2...\n")

beta_mu_2    <- c(Intercept = 0.5,  V1 = -0.50, V2 =  1.48, V3 = -1.45, V4 =  0.53, V5 = 0)
beta_sigma_2 <- c(Intercept = 0.0,  V1 =  0.00, V2 =  0.00, V3 = -0.69, V4 = -1.13, V5 = 0.91)

eta_mu_2    <- beta_mu_2["Intercept"]    + X_mat %*% beta_mu_2[paste0("V", 1:p)]
eta_sigma_2 <- beta_sigma_2["Intercept"] + X_mat %*% beta_sigma_2[paste0("V", 1:p)]

mu_2    <- pmax(pmin(plogis(eta_mu_2),    1 - 1e-6), 1e-6)
sigma_2 <- exp(eta_sigma_2)

y2 <- PLreg::rPL(n = n, mu = mu_2, sigma = sigma_2,
                 lambda = lambda_true, family = fam, zeta = zeta)
y2 <- pmax(pmin(y2, 1 - 1e-6), 1e-6)

data_pl2 <- data.frame(y = y2, X)

# --- Boosting ---
model_boost2 <- glmboostLSS(
  formula  = list(mu     = y ~ V1 + V2 + V3 + V4,
                  sigma  = y ~ V3 + V4 + V5,
                  lambda = y ~ 1),
  data     = data_pl2,
  families = family_pl,
  control  = boost_control(mstop = mstop, nu = nu, trace = FALSE)
)

b2_mu     <- normalise_intercept(coef(model_boost2, parameter = "mu"))
b2_sigma  <- normalise_intercept(coef(model_boost2, parameter = "sigma"))
b2_lambda <- unlist(coef(model_boost2, parameter = "lambda"))  # already log scale

# --- PLreg MLE ---
model_reg2 <- PLreg::PLreg(
  y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5,
  data   = data_pl2,
  family = fam,
  lambda = lambda_true,
  zeta   = zeta
)
all2      <- coef(model_reg2)
mle2_mu   <- all2[!grepl("^\\(sigma\\)", names(all2))]
mle2_sig  <- setNames(all2[grepl("^\\(sigma\\)", names(all2))],
                      gsub("^\\(sigma\\)_", "", names(all2)[grepl("^\\(sigma\\)", names(all2))]))

true_mu_2    <- c("(Intercept)" =  0.50, V1 = -0.50, V2 =  1.48, V3 = -1.45, V4 =  0.53)
true_sigma_2 <- c("(Intercept)" =  0.00, V3 = -0.69, V4 = -1.13, V5 =  0.91)

print_table(
  boost_mu     = b2_mu,     boost_sigma  = b2_sigma,  boost_lambda = b2_lambda,
  mle_mu       = mle2_mu,   mle_sigma    = mle2_sig,  mle_lambda   = log(lambda_true),
  true_mu      = true_mu_2, true_sigma   = true_sigma_2, true_lambda = log(lambda_true),
  label        = "Scenario 2 (PL-NO, n=100, p=5)"
)

make_comparison_plot(
  boost_mu     = b2_mu,     boost_sigma  = b2_sigma,  boost_lambda = b2_lambda,
  mle_mu       = mle2_mu,   mle_sigma    = mle2_sig,  mle_lambda   = log(lambda_true),
  true_mu      = true_mu_2, true_sigma   = true_sigma_2, true_lambda = log(lambda_true),
  scenario_label = "Scenario 2 — PL-NO  (n = 100, p = 5)",
  out_path       = file.path(out_dir, "diss_PL_scenario2.png")
)

cat("\n>>> Done. Plots saved to:", out_dir, "\n")