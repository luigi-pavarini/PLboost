# Template: boosted Power-Logit (PL) regression
#
# Run from the repository root, e.g.
#   Rscript examples/plboost_template.R

rm(list = ls())

suppressPackageStartupMessages({
  library(mboost)
  library(gamboostLSS)
})

# Needs to be settled!!!
source(".../essentials/v_function.r")
source(".../essentials/families_PL.R")

# --- data ---

# df <- ... (a data.frame with response + covariates)
# PL is for y in (0,1). If your data has 0/1, use IPL.
# y <- ...  # must be your response variable (e.g. "lgd")

# --- settings ---

fam  <- "NO"  # "NO","LO","TF","PE","Hyp","SN","SLASH" # must be a valid family name
zeta <- 2 # must be > 0

pl_family <- PowerLogit(
  family = fam,
  zeta = zeta,
  stabilization = "MAD" # "none", "MAD", "L2"
)

# --- formulas ---

# Option A (simple): use all columns in df except the response
# mu_formula    <- as.formula(paste(y_name, "~ ."))
# sigma_formula <- as.formula(paste(y_name, "~ ."))
#
# Option B: specify subsets / different covariates for each parameter
stop("TODO: define `mu_formula` and `sigma_formula` (and keep `lambda` as intercept-only).")

mu_formula    <- NULL # e.g. "y ~ V1 + V2 + V3 + V4"
sigma_formula <- NULL # e.g. "y ~ V3 + V4 + V5"

# --- fit ---

# There are many options for the control parameters. Here are some common ones:
# - mstop: vector of length 3 (mu, sigma, lambda)
# - nu: between 0 and 1
# - trace: TRUE or FALSE (follow up on the progress of the boosting algorithm)
# - ...: other options look at the help of mboost and/or gamboostLSS

fit <- glmboostLSS(
  formula = list(
    mu     = mu_formula,
    sigma  = sigma_formula,
    lambda = y ~ 1 # keep lambda as intercept-only
  ),
  data     = df,
  families = pl_family,
  control  = boost_control(
    mstop = c(mu = 500, sigma = 500, lambda = 50), # must be a vector of length 3
    nu = 0.01, # must be between 0 and 1
    trace = TRUE # TRUE or FALSE
  )
)

print(mstop(fit))
print(coef(fit, parameter = "mu"))
print(coef(fit, parameter = "sigma"))
print(exp(unlist(coef(fit, parameter = "lambda"))))

