library(gamboostLSS)
library(mboost)
library(faux)
library(PLreg)

# Needs to be settled!!!
base_path <- "..."
source(file.path(base_path, "\essentials\v_function.r"))
source(file.path(base_path, "\essentials\families_PL.R"))

# ==============================================================================
# SETTINGS
# ==============================================================================

set.seed(123)
n           <- 100
p           <- 5
fam         <- "NO" #"NO","LO","TF","PE","Hyp","SN","SLASH"
zeta        <- 2
lambda_true <- 1.2

# ==============================================================================
# COVARIATES
# ==============================================================================

X <- rnorm_multi(n, vars = p,
                 mu = rep(0, p), sd = rep(1, p), r = 0.5,
                 varnames = paste0("V", 1:p), empirical = FALSE)
X_mat <- as.matrix(X)

# ==============================================================================
# SCENARIO 1
# ==============================================================================

beta_mu    <- c(-0.99,  0.49, -1.17,  0.38, 0)
beta_sigma <- c(0, 0, -0.41, -0.24,  0.81)

eta_mu    <-  0.5 + X_mat %*% beta_mu
eta_sigma <-  0.0 + X_mat %*% beta_sigma

mu    <- pmax(pmin(plogis(eta_mu), 1 - 1e-6), 1e-6)
sigma <- exp(eta_sigma)

y <- PLreg::rPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true, family = fam, zeta = zeta)
hist(y)
data_pl <- data.frame(y, X)

# PL Boosting
family_pl   <- PowerLogit(family = fam, stabilization = "MAD", zeta = zeta)
model_boost <- glmboostLSS(
  formula  = list(mu     = y ~ V1 + V2 + V3 + V4,
                  sigma  = y ~ V3 + V4 + V5,
                  lambda = y ~ 1),
  data     = data_pl,
  families = family_pl,
  control  = boost_control(mstop = c(1000, 1000, 1000), nu = 0.01, trace = T)
)

coef(model_boost, parameter = "mu")
coef(model_boost, parameter = "sigma")
exp(coef(model_boost, parameter = "lambda"))

# PL Regression
model_reg <- PLreg::PLreg(
  y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5,
  data = data_pl, family = fam, lambda = lambda_true, zeta = zeta
)
summary(model_reg)

# ==============================================================================
# SCENARIO 2
# ==============================================================================

beta_mu    <- c(-0.5,  1.48, -1.45,  0.53, 0)
beta_sigma <- c(0, 0, -0.69, -1.13,  0.91)

eta_mu    <- 0.5 + X_mat %*% beta_mu
eta_sigma <- 0.0 + X_mat %*% beta_sigma

mu    <- pmax(pmin(plogis(eta_mu), 1 - 1e-6), 1e-6)
sigma <- exp(eta_sigma)

y <- PLreg::rPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true, family = fam, zeta = zeta)
y <- pmax(pmin(y, 1 - 1e-6), 1e-6)
hist(y)
summary(y)
data_pl <- data.frame(y, X)

# PL Boosting
model_boost <- glmboostLSS(
  formula  = list(mu     = y ~ V1 + V2 + V3 + V4,
                  sigma  = y ~ V3 + V4 + V5,
                  lambda = y ~ 1),
  data     = data_pl,
  families = family_pl,
  control  = boost_control(mstop = c(1000, 1000, 1000), nu = 0.01, trace = T)
)

coef(model_boost, parameter = "mu")
coef(model_boost, parameter = "sigma")
exp(unlist(coef(model_boost, parameter = "lambda")))

# PL Regression
model_reg <- PLreg::PLreg(
  y ~ V1 + V2 + V3 + V4| V3 + V4 + V5,
  data = data_pl, family = fam, lambda = lambda_true, zeta = zeta
)
summary(model_reg)