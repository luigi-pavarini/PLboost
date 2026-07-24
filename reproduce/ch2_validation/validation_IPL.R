# ==============================================================================
# testing_IPL_diss.R
# Dissertation validation script: Boosting (IPL) vs IPLreg (PLreg MLE + GLM)
#
# For each of the three inflation scenarios, fits both estimators and produces:
#   - Side-by-side coefficient comparison plot across all 4 submodels (PNG)
#   - Printed comparison table in the console
#
# Requires: gamboostLSS, mboost, faux, PLreg, Formula + your local source files
# ==============================================================================

library(PLboost)
library(gamboostLSS)
library(mboost)
library(faux)
library(PLreg)
library(Formula)
base_path <- NULL  # not needed: the implementation now ships as the PLboost package
out_dir <- Sys.getenv("PLBOOST_OUT", unset = ".")

# ==============================================================================
# SETTINGS
# ==============================================================================

set.seed(123)
n           <- 3000
p           <- 5
p_alpha     <- 2
fam         <- "NO"
zeta        <- 2
lambda_true <- 1.2
c_infl      <- 1
nu          <- 0.01
mstop       <- 1000

img_width  <- 2000*1.8
img_height <- 2000
img_res    <- 330

family_ipl <- IPL(family = fam, stabilization = "MAD", c = c_infl, zeta = zeta)

# ==============================================================================
# IPLreg: PLreg (continuous component) + binomial GLM (inflation component)
# ==============================================================================

# PLreg stores coefficients with submodel prefix e.g. "(mu)_V1", "(sigma)_V3"
# Strip prefix so names match variable names used in boosting
normalise_intercept <- function(x) {
  names(x) <- gsub("^.*[Ii]ntercept.*$", "(Intercept)", names(x))
  x
}

IPLreg <- function(formula, data, family, zeta = 2, type = "ML", c) {
  oformula <- Formula(formula)
  mf       <- model.frame(oformula, data = data)
  y        <- model.response(mf)
  ind      <- ifelse(y == c, 1, 0)
  
  data.PL  <- data[ind == 0, ]
  data.glm <- cbind(data, ind)
  
  formula.PL  <- formula(oformula, lhs = 1, rhs = -3)
  formula.glm <- update(formula(oformula, lhs = 1, rhs = 3), ind ~ .)
  
  fit.PL  <- PLreg(formula.PL, data = data.PL, family = family,
                   zeta = zeta, type = type)
  fit.glm <- glm(formula.glm, family = binomial(), data = data.glm)
  
  return(list(fitPL = fit.PL, fitglm = fit.glm))
}

# ==============================================================================
# COVARIATES (shared across all three scenarios)
# ==============================================================================

X <- rnorm_multi(n, vars = p, mu = rep(0, p), sd = rep(1, p), r = 0.5,
                 varnames = paste0("V", 1:p), empirical = FALSE)
X_mat   <- as.matrix(X)
X_alpha <- matrix(runif(n * p_alpha), nrow = n)
colnames(X_alpha) <- paste0("U", 1:p_alpha)

beta_mu    <- c(-0.99,  0.49, -1.17,  0.38, rep(0, p - 4))
beta_sigma <- c( 0.00,  0.00, -0.41, -0.24, 0.81)

eta_mu    <-  0.5 + X_mat %*% beta_mu
eta_sigma <-  0.0 + X_mat %*% beta_sigma

mu    <- pmax(pmin(plogis(eta_mu),  1 - 1e-6), 1e-6)
sigma <- exp(eta_sigma)

# true values (with intercepts) for display
true_mu    <- c("(Intercept)" =  0.50, V1 = -0.99, V2 =  0.49, V3 = -1.17, V4 =  0.38)
true_sigma <- c("(Intercept)" =  0.00, V3 = -0.41, V4 = -0.24, V5 =  0.81)

# ==============================================================================
# HELPERS
# ==============================================================================

make_ipl_plot <- function(boost_mu, boost_sigma, boost_lambda, boost_alpha,
                          mle_mu,   mle_sigma,   mle_lambda,   mle_alpha,
                          true_mu,  true_sigma,  true_lambda,  true_alpha,
                          scenario_label, out_path) {
  
  # colours matching application plot palette
  col_boost      <- adjustcolor("#E07B39", alpha.f = 0.80)  # orange  — Boosting
  col_boost_fill <- adjustcolor("#E07B39", alpha.f = 0.45)
  col_mle        <- adjustcolor("#3DA85A", alpha.f = 0.80)  # green   — MLE
  col_mle_fill   <- adjustcolor("#3DA85A", alpha.f = 0.45)
  col_true       <- "red4"
  col_bord       <- "#2a2a2a"
  bar_half       <- 0.28
  
  submodels <- list(
    list(name  = "mu",
         vars  = names(true_mu),
         true  = true_mu,    boost = boost_mu,    mle = mle_mu,
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
         xlab  = expression(log(lambda))),
    list(name  = "alpha",
         vars  = names(true_alpha),
         true  = true_alpha, boost = boost_alpha, mle = mle_alpha,
         xlab  = expression("Coefficient value  " * eta[alpha]))
  )
  
  base <- sub("\\.png$", "", out_path)
  
  for (sm in submodels) {
    panel_path <- paste0(base, "_", sm$name, ".png")
    
    vars   <- sm$vars
    n_vars <- length(vars)
    
    true_vals  <- as.numeric(sm$true[vars])
    boost_vals <- as.numeric(sm$boost[vars])
    mle_vals   <- as.numeric(sm$mle[vars])
    
    y_boost <- seq(n_vars, 1, by = -1) * 3 - 1
    y_mle   <- y_boost - 1
    y_label <- (y_boost + y_mle) / 2
    
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
    
    abline(v = 0, col = "gray60", lty = 2, lwd = 1)
    
    for (k in seq_along(vars)) {
      segments(true_vals[k], y_mle[k]   - bar_half,
               true_vals[k], y_boost[k] + bar_half,
               col = col_true, lty = 2, lwd = 1.8)
    }
    
    rect(xleft   = pmin(0, boost_vals),
         ybottom = y_boost - bar_half,
         xright  = pmax(0, boost_vals),
         ytop    = y_boost + bar_half,
         col     = col_boost_fill,
         border  = col_boost,
         lwd     = 0.8)
    
    rect(xleft   = pmin(0, mle_vals),
         ybottom = y_mle - bar_half,
         xright  = pmax(0, mle_vals),
         ytop    = y_mle + bar_half,
         col     = col_mle_fill,
         border  = col_mle,
         lwd     = 0.8)
    
    points(boost_vals, y_boost, pch = 22, col = col_bord,
           bg = "white", cex = 1.1, lwd = 1.0)
    points(mle_vals,   y_mle,   pch = 22, col = col_bord,
           bg = "white", cex = 1.1, lwd = 1.0)
    
    axis(2, at = y_label, labels = vars, las = 1, cex.axis = 0.80)
    
    legend("bottomright",
           legend = c("True value", "Boosting", "IPLreg MLE"),
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

print_ipl_table <- function(boost_mu, boost_sigma, boost_lambda, boost_alpha,
                            mle_mu,   mle_sigma,   mle_lambda,   mle_alpha,
                            true_mu,  true_sigma,  true_lambda,  true_alpha,
                            label) {
  
  cat("\n", strrep("=", 66), "\n")
  cat(" IPL Coefficient comparison:", label, "\n")
  cat(strrep("=", 66), "\n")
  cat(sprintf("%-14s  %8s  %8s  %8s\n", "Parameter", "True", "Boosting", "IPLreg"))
  cat(strrep("-", 46), "\n")
  
  for (v in names(true_mu)) {
    bv <- if (v %in% names(boost_mu)) boost_mu[v] else NA
    mv <- if (v %in% names(mle_mu))   mle_mu[v]   else NA
    cat(sprintf("mu: %-10s  %8.4f  %8.4f  %8.4f\n", v, true_mu[v], bv, mv))
  }
  for (v in names(true_sigma)) {
    bv <- if (v %in% names(boost_sigma)) boost_sigma[v] else NA
    mv <- if (v %in% names(mle_sigma))   mle_sigma[v]   else NA
    cat(sprintf("sigma: %-7s  %8.4f  %8.4f  %8.4f\n", v, true_sigma[v], bv, mv))
  }
  cat(sprintf("%-14s  %8.4f  %8.4f  %8.4f\n",
              "lambda", true_lambda, boost_lambda, mle_lambda))
  for (v in names(true_alpha)) {
    bv <- if (v %in% names(boost_alpha)) boost_alpha[v] else NA
    mv <- if (v %in% names(mle_alpha))   mle_alpha[v]   else NA
    cat(sprintf("alpha: %-7s  %8.4f  %8.4f  %8.4f\n", v, true_alpha[v], bv, mv))
  }
  cat(strrep("=", 66), "\n\n")
}

# ==============================================================================
# SCENARIO 1 — Low inflation: alpha in [0.047, 0.142]
# ==============================================================================

cat("\n>>> Fitting Scenario 1 — Low inflation...\n")

beta_alpha_1 <- c(0.7, -0.5)
eta_alpha_1  <- -2.5 + X_alpha %*% beta_alpha_1
alpha_1      <- plogis(eta_alpha_1)
true_alpha_1 <- c("(Intercept)" = -2.5, U1 =  0.7, U2 = -0.5)

y1 <- rIPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true,
           alpha = alpha_1, family = fam, zeta = zeta, c = c_infl)
data_ipl1 <- data.frame(y = y1, X_alpha, X)

# Boosting
model_boost1 <- glmboostLSS(
  formula = list(mu     = y ~ V1 + V2 + V3 + V4,
                 sigma  = y ~ V3 + V4 + V5,
                 lambda = y ~ 1,
                 alpha  = y ~ U1 + U2),
  data     = data_ipl1,
  families = family_ipl,
  control  = boost_control(mstop = mstop, nu = nu, trace = FALSE)
)

b1_mu     <- normalise_intercept(coef(model_boost1, parameter = "mu"))
b1_sigma  <- normalise_intercept(coef(model_boost1, parameter = "sigma"))
b1_lambda <- unlist(coef(model_boost1, parameter = "lambda"))  # already log scale
b1_alpha  <- normalise_intercept(coef(model_boost1, parameter = "alpha"))

# IPLreg MLE
fit1      <- IPLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
                    data = data_ipl1, family = fam, zeta = zeta,
                    type = "ML", c = c_infl)
all1      <- coef(fit1$fitPL)
mle1_mu   <- all1[!grepl("^\\(sigma\\)", names(all1))]
mle1_sig  <- setNames(all1[grepl("^\\(sigma\\)", names(all1))],
                      gsub("^\\(sigma\\)_", "", names(all1)[grepl("^\\(sigma\\)", names(all1))]))
mle1_lam  <- log(fit1$fitPL$coefficients$skewness)
mle1_alph <- coef(fit1$fitglm)

print_ipl_table(
  boost_mu = b1_mu, boost_sigma = b1_sigma, boost_lambda = b1_lambda, boost_alpha = b1_alpha,
  mle_mu   = mle1_mu, mle_sigma = mle1_sig, mle_lambda  = mle1_lam,  mle_alpha   = mle1_alph,
  true_mu  = true_mu, true_sigma = true_sigma, true_lambda = log(lambda_true), true_alpha = true_alpha_1,
  label = "Scenario 1 — Low inflation  [0.047, 0.142]"
)

make_ipl_plot(
  boost_mu = b1_mu, boost_sigma = b1_sigma, boost_lambda = b1_lambda, boost_alpha = b1_alpha,
  mle_mu   = mle1_mu, mle_sigma = mle1_sig, mle_lambda  = mle1_lam,  mle_alpha   = mle1_alph,
  true_mu  = true_mu, true_sigma = true_sigma, true_lambda = log(lambda_true), true_alpha = true_alpha_1,
  scenario_label = "Scenario 1 — IPL-NO  |  Low inflation  alpha in [0.047, 0.142]  (n = 3000, p = 5)",
  out_path = file.path(out_dir, "diss_IPL_scenario1_low.png")
)

# ==============================================================================
# SCENARIO 2 — Medium inflation: alpha in [0.142, 0.310]
# ==============================================================================

cat("\n>>> Fitting Scenario 2 — Medium inflation...\n")

beta_alpha_2 <- c(0.70, -0.3)
eta_alpha_2  <- -1.5 + X_alpha %*% beta_alpha_2
alpha_2      <- plogis(eta_alpha_2)
true_alpha_2 <- c("(Intercept)" = -1.5, U1 =  0.70, U2 = -0.3)

y2 <- rIPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true,
           alpha = alpha_2, family = fam, zeta = zeta, c = c_infl)
data_ipl2 <- data.frame(y = y2, X_alpha, X)

model_boost2 <- glmboostLSS(
  formula = list(mu     = y ~ V1 + V2 + V3 + V4,
                 sigma  = y ~ V3 + V4 + V5,
                 lambda = y ~ 1,
                 alpha  = y ~ U1 + U2),
  data     = data_ipl2,
  families = family_ipl,
  control  = boost_control(mstop = mstop, nu = nu, trace = FALSE)
)

b2_mu     <- normalise_intercept(coef(model_boost2, parameter = "mu"))
b2_sigma  <- normalise_intercept(coef(model_boost2, parameter = "sigma"))
b2_lambda <- unlist(coef(model_boost2, parameter = "lambda"))  # already log scale
b2_alpha  <- normalise_intercept(coef(model_boost2, parameter = "alpha"))

fit2      <- IPLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
                    data = data_ipl2, family = fam, zeta = zeta,
                    type = "ML", c = c_infl)
all2      <- coef(fit2$fitPL)
mle2_mu   <- all2[!grepl("^\\(sigma\\)", names(all2))]
mle2_sig  <- setNames(all2[grepl("^\\(sigma\\)", names(all2))],
                      gsub("^\\(sigma\\)_", "", names(all2)[grepl("^\\(sigma\\)", names(all2))]))
mle2_lam  <- log(fit2$fitPL$coefficients$skewness)
mle2_alph <- coef(fit2$fitglm)

print_ipl_table(
  boost_mu = b2_mu, boost_sigma = b2_sigma, boost_lambda = b2_lambda, boost_alpha = b2_alpha,
  mle_mu   = mle2_mu, mle_sigma = mle2_sig, mle_lambda  = mle2_lam,  mle_alpha   = mle2_alph,
  true_mu  = true_mu, true_sigma = true_sigma, true_lambda = log(lambda_true), true_alpha = true_alpha_2,
  label = "Scenario 2 — Medium inflation  [0.142, 0.310]"
)

make_ipl_plot(
  boost_mu = b2_mu, boost_sigma = b2_sigma, boost_lambda = b2_lambda, boost_alpha = b2_alpha,
  mle_mu   = mle2_mu, mle_sigma = mle2_sig, mle_lambda  = mle2_lam,  mle_alpha   = mle2_alph,
  true_mu  = true_mu, true_sigma = true_sigma, true_lambda = log(lambda_true), true_alpha = true_alpha_2,
  scenario_label = "Scenario 2 — IPL-NO  |  Medium inflation  alpha in [0.142, 0.310]  (n = 3000, p = 5)",
  out_path = file.path(out_dir, "diss_IPL_scenario2_medium.png")
)

# ==============================================================================
# SCENARIO 3 — High inflation: alpha in [0.231, 0.401]
# ==============================================================================

cat("\n>>> Fitting Scenario 3 — High inflation...\n")

beta_alpha_3 <- c(0.4, -0.4)
eta_alpha_3  <- -0.8 + X_alpha %*% beta_alpha_3
alpha_3      <- plogis(eta_alpha_3)
true_alpha_3 <- c("(Intercept)" = -0.8, U1 =  0.4, U2 = -0.4)

y3 <- rIPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true,
           alpha = alpha_3, family = fam, zeta = zeta, c = c_infl)
data_ipl3 <- data.frame(y = y3, X_alpha, X)

model_boost3 <- glmboostLSS(
  formula = list(mu     = y ~ V1 + V2 + V3 + V4,
                 sigma  = y ~ V3 + V4 + V5,
                 lambda = y ~ 1,
                 alpha  = y ~ U1 + U2),
  data     = data_ipl3,
  families = family_ipl,
  control  = boost_control(mstop = mstop, nu = nu, trace = FALSE)
)

b3_mu     <- normalise_intercept(coef(model_boost3, parameter = "mu"))
b3_sigma  <- normalise_intercept(coef(model_boost3, parameter = "sigma"))
b3_lambda <- unlist(coef(model_boost3, parameter = "lambda"))  # already log scale
b3_alpha  <- normalise_intercept(coef(model_boost3, parameter = "alpha"))

fit3      <- IPLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
                    data = data_ipl3, family = fam, zeta = zeta,
                    type = "ML", c = c_infl)
all3      <- coef(fit3$fitPL)
mle3_mu   <- all3[!grepl("^\\(sigma\\)", names(all3))]
mle3_sig  <- setNames(all3[grepl("^\\(sigma\\)", names(all3))],
                      gsub("^\\(sigma\\)_", "", names(all3)[grepl("^\\(sigma\\)", names(all3))]))
mle3_lam  <- log(fit3$fitPL$coefficients$skewness)
mle3_alph <- coef(fit3$fitglm)

print_ipl_table(
  boost_mu = b3_mu, boost_sigma = b3_sigma, boost_lambda = b3_lambda, boost_alpha = b3_alpha,
  mle_mu   = mle3_mu, mle_sigma = mle3_sig, mle_lambda  = mle3_lam,  mle_alpha   = mle3_alph,
  true_mu  = true_mu, true_sigma = true_sigma, true_lambda = log(lambda_true), true_alpha = true_alpha_3,
  label = "Scenario 3 — High inflation  [0.231, 0.401]"
)

make_ipl_plot(
  boost_mu = b3_mu, boost_sigma = b3_sigma, boost_lambda = b3_lambda, boost_alpha = b3_alpha,
  mle_mu   = mle3_mu, mle_sigma = mle3_sig, mle_lambda  = mle3_lam,  mle_alpha   = mle3_alph,
  true_mu  = true_mu, true_sigma = true_sigma, true_lambda = log(lambda_true), true_alpha = true_alpha_3,
  scenario_label = "Scenario 3 — IPL-NO  |  High inflation  alpha in [0.231, 0.401]  (n = 3000, p = 5)",
  out_path = file.path(out_dir, "diss_IPL_scenario3_high.png")
)

cat("\n>>> Done. Plots saved to:", out_dir, "\n")