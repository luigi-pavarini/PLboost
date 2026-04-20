library(gamboostLSS)
library(mboost)
library(faux)
library(PLreg)
library(Formula)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else NULL
script_dir <- if (!is.null(script_path)) dirname(normalizePath(script_path)) else getwd()
essentials_dir <- file.path(script_dir, "essentials")
source(file.path(essentials_dir, "v_function.r"))
source(file.path(essentials_dir, "IPL_distributions.R"))
source(file.path(essentials_dir, "families_IPL.R"))
source(file.path(essentials_dir, "families_PL.R"))

# ==============================================================================
# SETTINGS
# ==============================================================================

set.seed(123)
n           <- 3000
p           <- 5
p_alpha     <- 2
fam         <- "NO" #"NO","LO","TF","PE","Hyp","SN","SLASH"
zeta        <- 2
lambda_true <- 1.2
c_infl      <- 1   # 0 for zero-inflation, 1 for one-inflation

# ==============================================================================
# IPLreg: PLreg (continuous) + binomial GLM (discrete)
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
  
  return(list(fitPL = fit.PL, fitglm = fit.glm))
}

# ==============================================================================
# COVARIATES
# ==============================================================================

X <- rnorm_multi(n, vars = p,
                 mu = rep(0, p), sd = rep(1, p), r = 0.5,
                 varnames = paste0("V", 1:p), empirical = FALSE)
X_mat   <- as.matrix(X)
X_alpha <- matrix(runif(n * p_alpha), nrow = n)
colnames(X_alpha) <- paste0("U", 1:p_alpha)

# ==============================================================================
# BETA COEFFICIENTS
# ==============================================================================

beta_mu    <- c(-0.99, 0.49, -1.17, 0.38, rep(0, p - 4))
beta_sigma <- c(0, 0, -0.41, -0.24, 0.81, rep(0, p - 5))

eta_mu    <-  0.5 + X_mat %*% beta_mu
eta_sigma <-  0.0 + X_mat %*% beta_sigma

mu    <- pmax(pmin(plogis(eta_mu), 1 - 1e-6), 1e-6)
sigma <- exp(eta_sigma)

# ==============================================================================
# SCENARIO 1 (Low): alpha in [0.047, 0.142]
# ==============================================================================

beta_alpha <- c(0.7, -0.5)
eta_alpha  <- -2.5 + X_alpha %*% beta_alpha
alpha      <- plogis(eta_alpha)
hist(alpha)
y <- rIPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true,
          alpha = alpha, family = fam, zeta = zeta, c = c_infl)
hist(y)
data_ipl <- data.frame(y, X_alpha, X)

# IPL Boosting
family_ipl  <- IPL(family = fam, stabilization = "MAD", c = c_infl, zeta = zeta)
model_boost <- glmboostLSS(
  formula = list(mu     = y ~ V1 + V2 + V3 + V4,
                 sigma  = y ~ V3 + V4 + V5,
                 lambda = y ~ 1,
                 alpha  = y ~ U1 + U2),
  data     = data_ipl,
  families = family_ipl,
  control  = boost_control(mstop = c(1000, 1000, 1000, 1000), nu = 0.01, trace = T)
)

coef(model_boost, parameter = "mu")
coef(model_boost, parameter = "sigma")
exp(unlist(coef(model_boost, parameter = "lambda")))
coef(model_boost, parameter = "alpha")

# IPL Regression
fit1 <- IPLreg(
  formula = y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
  data = data_ipl, family = fam, zeta = zeta, type = "ML", c = c_infl
)
fit1$fitPL
fit1$fitglm

# ==============================================================================
# SCENARIO 2 (Medium): alpha in [0.142, 0.310]
# ==============================================================================

beta_alpha <- c(0.70, -0.3)
eta_alpha  <- -1.5 + X_alpha %*% beta_alpha
alpha      <- plogis(eta_alpha)

y <- rIPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true,
          alpha = alpha, family = fam, zeta = zeta, c = c_infl)

data_ipl <- data.frame(y, X_alpha, X)

# IPL Boosting
model_boost <- glmboostLSS(
  formula = list(mu     = y ~ V1 + V2 + V3 + V4,
                 sigma  = y ~ V3 + V4 + V5,
                 lambda = y ~ 1,
                 alpha  = y ~ U1 + U2),
  data     = data_ipl,
  families = family_ipl,
  control  = boost_control(mstop = c(1000, 1000, 1000, 1000), nu = 0.01, trace = T)
)

coef(model_boost, parameter = "mu")
coef(model_boost, parameter = "sigma")
exp(unlist(coef(model_boost, parameter = "lambda")))
coef(model_boost, parameter = "alpha")

# IPL Regression
fit2 <- IPLreg(
  formula = y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
  data = data_ipl, family = fam, zeta = zeta, type = "ML", c = c_infl
)
fit2$fitPL
fit2$fitglm

# ==============================================================================
# SCENARIO 3 (High): alpha in [0.231, 0.401]
# ==============================================================================

beta_alpha <- c(0.4, -0.4)
eta_alpha  <- -0.8 + X_alpha %*% beta_alpha
alpha      <- plogis(eta_alpha)

y <- rIPL(n = n, mu = mu, sigma = sigma, lambda = lambda_true,
          alpha = alpha, family = fam, zeta = zeta, c = c_infl)

data_ipl <- data.frame(y, X_alpha, X)

# IPL Boosting
model_boost <- glmboostLSS(
  formula = list(mu     = y ~ V1 + V2 + V3 + V4,
                 sigma  = y ~ V3 + V4 + V5,
                 lambda = y ~ 1,
                 alpha  = y ~ U1 + U2),
  data     = data_ipl,
  families = family_ipl,
  control  = boost_control(mstop = c(1000, 1000, 1000, 1000), nu = 0.01, trace = T)
)

coef(model_boost, parameter = "mu")
coef(model_boost, parameter = "sigma")
exp(unlist(coef(model_boost, parameter = "lambda")))
coef(model_boost, parameter = "alpha")

# IPL Regression
fit3 <- IPLreg(
  formula = y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2,
  data = data_ipl, family = fam, zeta = zeta, type = "ML", c = c_infl
)
fit3$fitPL
fit3$fitglm

