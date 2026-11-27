# ==============================================================================
# example_residuals.R
#
# Quantile residuals with simulated envelopes for PL and APL fits, by maximum
# likelihood and by boosting, and worm plots of the boosting fits, on data
# simulated as in the validation of Sections 3.1 and 3.2. Experimental, not
# part of the master thesis.
# Run:  Rscript example_residuals.R   (draws to Rplots.pdf when not interactive)
# ==============================================================================

library(PLboost)
library(gamboostLSS)
library(mboost)
library(faux)
library(PLreg)
library(Formula)

source(system.file("diagnostics", "residuals_boosting.R", package = "PLboost"))
source(system.file("diagnostics", "wp_PLboost.R", package = "PLboost"))

# ==============================================================================
# SETTINGS
# ==============================================================================
# Same configuration as validation_multirep.R, with one replication instead of
# M = 100: a residual plot belongs to one fit.

semente     <- 1206
set.seed(semente)
fam         <- "NO"
zeta        <- 2
lambda_true <- 1.2
nu          <- 0.01
p           <- 5
rep_env     <- 40     # envelope replications; each one is a full refit
conf        <- 0.95

beta_mu    <- c(-0.99, 0.49, -1.17, 0.38, 0)
beta_sigma <- c(0, 0, -0.41, -0.24, 0.81)

summary_env <- function(e)
  cat(sprintf("%-16s mean %7.4f  sd %6.4f  outside the band %5.1f%%\n",
              e$main, mean(e$res), sd(e$res), 100 * e$prop_outside))

# ==============================================================================
# ==============================================================================
# PL
# ==============================================================================
# ==============================================================================

n_pl     <- 500
mstop_pl <- 800

# ==============================================================================
# PL DATA
# ==============================================================================

X <- rnorm_multi(n_pl, vars = p, mu = rep(0, p), sd = rep(1, p), r = 0.5,
                 varnames = paste0("V", 1:p), empirical = FALSE)
X_mat <- as.matrix(X)

#mu    <- pmax(pmin(plogis(0.5 + X_mat %*% beta_mu), 1 - 1e-6), 1e-6)
mu    <- plogis(0.5 + X_mat %*% beta_mu)

sigma <- exp(0.0 + X_mat %*% beta_sigma)
y     <- rPL(n_pl, mu, sigma, lambda_true, family = fam, zeta = zeta)
#y     <- pmax(pmin(y, 1 - 1e-6), 1e-6)
data_pl <- data.frame(y, X)

# ==============================================================================
# PL MAXIMUM LIKELIHOOD
# ==============================================================================

fit_pl_mle <- PLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5, data = data_pl,
                    family = fam, zeta = zeta)

env_pl_mle <- envelope_PLboost(fit_pl_mle, data_pl, rep = rep_env, conf = conf,
                               main = "PL - MLE", plot = FALSE)

# ==============================================================================
# PL BOOSTING
# ==============================================================================

fit_pl_boost <- glmboostLSS(
  formula  = list(mu = y ~ V1 + V2 + V3 + V4, sigma = y ~ V3 + V4 + V5,
                  lambda = y ~ 1),
  data     = data_pl,
  families = PowerLogit(family = fam, stabilization = "MAD", zeta = zeta),
  control  = boost_control(mstop = mstop_pl, nu = nu, trace = FALSE))

env_pl_boost <- envelope_PLboost(fit_pl_boost, data_pl, family = fam, zeta = zeta,
                                 rep = rep_env, conf = conf,
                                 main = "PL - boosting", plot = FALSE)

# ==============================================================================
# PL PLOTS AND RESULTS
# ==============================================================================

par(mfrow = c(1, 2))
plot_envelope(env_pl_mle)
plot_envelope(env_pl_boost)
par(mfrow = c(1, 1))

summary_env(env_pl_mle)
summary_env(env_pl_boost)

# ==============================================================================
# PL WORM PLOT (boosting)
# ==============================================================================

wp_pl_boost <- plot_worm(env_pl_boost$res, main = "PL - boosting")
cat(sprintf("PL - boosting    outside the worm plot band %5.1f%%\n", 100 * wp_pl_boost))

# ==============================================================================
# ==============================================================================
# APL
# ==============================================================================
# ==============================================================================

n_apl      <- 800
mstop_apl  <- 1000
p_alpha    <- 2
c_infl     <- 1
beta_alpha <- c(0.70, -0.3)          # medium boundary mass, alpha in [0.142, 0.310]

# ==============================================================================
# IPLreg: APL by maximum likelihood, PLreg (continuous) + binomial GLM (discrete)
# ==============================================================================

IPLreg <- function(formula, data, family, zeta = 2, type = "ML", c) {
  oformula <- Formula(formula)
  mf       <- model.frame(oformula, data = data)
  y        <- model.response(mf)
  ind      <- ifelse(y == c, 1, 0)

  data.PL  <- data[ind == 0, ]
  data.glm <- cbind(data, ind)

  formula.PL  <- formula(oformula, lhs = 1, rhs = -3)
  formula.glm <- update(formula(oformula, lhs = 1, rhs = 3), ind ~ .)

  fit.PL  <- PLreg(formula.PL, data = data.PL, family = family, zeta = zeta, type = type)
  fit.glm <- glm(formula.glm, family = binomial(), data = data.glm)

  out <- list(fitPL = fit.PL, fitglm = fit.glm, formula = formula, data = data,
              c = c, f_mu = formula(oformula, lhs = 0, rhs = 1),
              f_sigma = formula(oformula, lhs = 0, rhs = 2))
  class(out) <- "IPLreg"
  out
}

# ==============================================================================
# APL DATA
# ==============================================================================

X <- rnorm_multi(n_apl, vars = p, mu = rep(0, p), sd = rep(1, p), r = 0.5,
                 varnames = paste0("V", 1:p), empirical = FALSE)
X_mat   <- as.matrix(X)
X_alpha <- matrix(runif(n_apl * p_alpha), nrow = n_apl)
colnames(X_alpha) <- paste0("U", 1:p_alpha)

mu    <- pmax(pmin(plogis(0.5 + X_mat %*% beta_mu), 1 - 1e-6), 1e-6)
sigma <- exp(0.0 + X_mat %*% beta_sigma)
alpha <- plogis(-1.5 + X_alpha %*% beta_alpha)
y     <- rAPL(n_apl, alpha = alpha, mu = mu, sigma = sigma, lambda = lambda_true,
              family = fam, zeta = zeta, c = c_infl)
data_apl <- data.frame(y, X_alpha, X)

# ==============================================================================
# APL MAXIMUM LIKELIHOOD
# ==============================================================================

fit_apl_mle <- IPLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
                      data = data_apl, family = fam, zeta = zeta,
                      type = "ML", c = c_infl)

env_apl_mle <- envelope_PLboost(fit_apl_mle, data_apl, c = c_infl,
                                rep = rep_env, conf = conf,
                                main = "APL - MLE", plot = FALSE)

# ==============================================================================
# APL BOOSTING
# ==============================================================================

fit_apl_boost <- glmboostLSS(
  formula  = list(mu = y ~ V1 + V2 + V3 + V4, sigma = y ~ V3 + V4 + V5,
                  lambda = y ~ 1, alpha = y ~ U1 + U2),
  data     = data_apl,
  families = APL(family = fam, stabilization = "MAD", c = c_infl, zeta = zeta),
  control  = boost_control(mstop = mstop_apl, nu = nu, trace = FALSE))

env_apl_boost <- envelope_PLboost(fit_apl_boost, data_apl, family = fam, zeta = zeta,
                                  c = c_infl, rep = rep_env, conf = conf,
                                  main = "APL - boosting", plot = FALSE)

# ==============================================================================
# APL PLOTS AND RESULTS
# ==============================================================================

par(mfrow = c(1, 2))
plot_envelope(env_apl_mle)
plot_envelope(env_apl_boost)
par(mfrow = c(1, 1))

summary_env(env_apl_mle)
summary_env(env_apl_boost)

# ==============================================================================
# APL WORM PLOT (boosting)
# ==============================================================================

wp_apl_boost <- plot_worm(env_apl_boost$res, main = "APL - boosting")
cat(sprintf("APL - boosting   outside the worm plot band %5.1f%%\n", 100 * wp_apl_boost))
