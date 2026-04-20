# Template: boosted Inflated Power-Logit (IPL) regression
#
# Run from the repository root, e.g.
#   Rscript examples/iplboost_template.R

rm(list = ls())

suppressPackageStartupMessages({
  library(mboost)
  library(gamboostLSS)
})

source(file.path("essentials", "v_function.r"))
source(file.path("essentials", "IPL_distributions.R"))
source(file.path("essentials", "families_IPL.R"))

# --- data ---

# df <- ... (a data.frame with response + covariates)
# IPL is for y in [0,1], with mass at c = 0 or c = 1.
# y <- ...  # must be your response variable (e.g. "lgd")

# --- settings ---

fam  <- "NO"  # "NO","LO","TF","PE","Hyp","SN","SLASH" # must be a valid family name
zeta <- 2 # must be > 0
c_infl <- 1 # 0 for mass at 0, 1 for mass at 1

ipl_family <- IPL(
  family = fam,
  zeta = zeta,
  c = c_infl,
  stabilization = "MAD" # "none", "MAD", "L2"
)

# --- formulas ---

# You can use all columns in df except the response:
# mu_formula    <- y ~ .
# sigma_formula <- y ~ .
# alpha_formula <- y ~ .
#
# Or specify subsets for each parameter:
mu_formula    <- NULL # e.g. "y ~ V1 + V2 + V3 + V4"
sigma_formula <- NULL # e.g. "y ~ V3 + V4 + V5"
alpha_formula <- NULL # e.g. "y ~ U1 + U2"

# --- fit ---

# There are many options for the control parameters. Here are some common ones:
# - mstop: vector of length 4 (mu, sigma, lambda, alpha)
# - nu: between 0 and 1
# - trace: TRUE or FALSE (follow up on the progress of the boosting algorithm)
# - ...: other options look at the help of mboost and/or gamboostLSS

fit <- glmboostLSS(
  formula = list(
    mu     = mu_formula,
    sigma  = sigma_formula,
    lambda = y ~ 1, # keep lambda as intercept-only
    alpha  = alpha_formula
  ),
  data     = df,
  families = ipl_family,
  control  = boost_control(
    mstop = c(mu = 500, sigma = 500, lambda = 50, alpha = 500), # must be a vector of length 4
    nu = 0.01, # must be between 0 and 1
    trace = TRUE # TRUE or FALSE
  )
)

print(mstop(fit))
print(coef(fit, parameter = "mu"))
print(coef(fit, parameter = "sigma"))
print(exp(unlist(coef(fit, parameter = "lambda"))))
print(coef(fit, parameter = "alpha"))

