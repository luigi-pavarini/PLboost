# Boosted IPL regression - Adapted with X_alpha uniform control approach

rm(list = ls())

# ==============================================================================
# SETUP
# ==============================================================================
# ---- setup -------------------------------------------------------------------
library(PLboost)
save_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # results go to $PLBOOST_OUT/IPL_simulations/
dir.create(file.path(save_path, "IPL_simulations"), showWarnings = FALSE, recursive = TRUE)

require("PLreg")
require("gamboostLSS")
require("mboost")
require("faux")
require("parallel")

# Parallel cross-validation: use 8 cores (mclapply reads options("mc.cores")).
options(mc.cores = 8)

# Canonical dissertation scenario: alpha intercept -1.5 (Section 3.2), labelled
# "medium". Keep this in sync with `example` in the summary and plot scripts.
example    = "medium"

M          = 100
n          = 800

p          = 1000             # Main covariates (high-dimensional)
fam        = 'NO' #"NO", "TF", "LO", "SN", "PE", "Hyp", "SLASH"
# Grid as reported in the dissertation (Section 3.2: "upper bound of mstop = 800").
mstop_init = 500
m_stop_cv  = 800
#mstop_init = 900   # experimental, does NOT match the dissertation
#m_stop_cv  = 1200
nu         = 0.01

# Grid for cross-validation (adapted for IPL with 4 parameters)
grid <- make.grid(max = c(mu = m_stop_cv, sigma = m_stop_cv, lambda = m_stop_cv, alpha = m_stop_cv),
                  min = mstop_init,
                  length.out = 2, dense_mu_grid = TRUE)

# ==============================================================================
# COVARIATES GENERATION - x1000 and x999 used for alpha control
# ==============================================================================

# Main covariates (normal) for all parameters — V1...V1000
# Alpha is controlled by V999 and V1000 (last two columns)
names_v <- paste0(rep("V", p), 1:p)
X <- rnorm_multi(n, vars = p, 
                 mu = rep(0, p),
                 sd = rep(1, p),
                 r = 0.5,
                 varnames = names_v,
                 empirical = FALSE)

# ==============================================================================
# TRUE COEFFICIENT VALUES - Adapted structure
# ==============================================================================

# Coefficients for main parameters (using 998 covariates)
beta_mu    <- c(-0.99, 0.49, -1.17, 0.38, rep(0, p - 4))
beta_sigma <- c(0, 0, -0.41, -0.24, 0.81, rep(0, p - 5))

# SCENARIO 1 (Low):
beta_alpha <- c(0.7, -0.5)
# eta_alpha = -2.5 + 0.7*V999 + (-0.5)*V1000
# Range: [-2.5+0+(-0.5), -2.5+0.7+0] = [-3.0, -1.8] -> alpha in [0.047, 0.142]

# SCENARIO 2 (Medium):
#beta_alpha <- c(0.70, -0.3)  # Extends proven second scenario
# eta_alpha = -1.5 + 0.70*V999 + (-0.3)*V1000
# Range: [-1.5+0-0.3, -1.5+0.7+0] = [-1.8, -0.8] -> alpha in [0.142, 0.310]

# SCENARIO 3 (High):
#beta_alpha <- c(0.4, -0.4)
# eta_alpha = -0.8 + 0.4*V999 + (-0.4)*V1000
# Range: [-0.8+0+(-0.4), -0.8+0.4+0] = [-1.2, -0.4] -> alpha in [0.231, 0.401]

# ==============================================================================
# LINEAR PREDICTORS
# ==============================================================================

X_mat     <- as.matrix(X)
eta_mu    <- 0.5 + X_mat %*% beta_mu
eta_sigma <- 0.0 + X_mat %*% beta_sigma

# Canonical scenario for the dissertation (Section 3.2): kappa_alpha = (-1.5, 0.70, -0.50),
# so alpha ~ logit^-1(-1.5) ~ 0.18 on average, ~20% of observations at the boundary.
eta_alpha <- -1.5 + X_mat[, c("V999", "V1000")] %*% beta_alpha

# Parameter transformations
mu     <- plogis(eta_mu)
mu     <- pmax(pmin(mu, 1 - 1e-6), 1e-6)
sigma  <- exp(eta_sigma)
lambda <- 1.2
alpha  <- plogis(eta_alpha)

# ==============================================================================
# STORAGE MATRICES
# Column order in data_ipl: y, V1...V1000
# => coef() vector order: (Intercept)[1], V1[2]...V999[1000], V1000[1001]
# ==============================================================================

var_names_full <- c("(Intercept)", names_v)
n_coef <- p + 1  # 1001

estimativas_mu     <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
estimativas_sigma  <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
estimativas_lambda <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
estimativas_alpha  <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))

# Selection frequency matrices (rows = intercept + all vars, cols = M runs)
freq_mu_mat    <- matrix(0L, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
freq_sigma_mat <- matrix(0L, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
freq_alpha_mat <- matrix(0L, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))

m_stop   <- matrix(0, nrow = 4, ncol = M)  # 4 parameters: mu, sigma, lambda, alpha
nvar_mu     <- NULL
nvar_sigma  <- NULL
nvar_lambda <- NULL
nvar_alpha  <- NULL

# ==============================================================================
# SIMULATION LOOP
# ==============================================================================

i = 1
k = 0
while(i <= M){
  tryCatch({
    k = k + 1  
    
    y <- rIPL(n = n, mu = mu, sigma = sigma, alpha = alpha,
              lambda = lambda, family = fam, zeta = 2, c = 0)
    
    # Combine covariates: V1...V1000 (alpha uses V999 and V1000)
    # Column order: y | V1...V1000
    data_ipl <- data.frame(y, X)
    
    # IPL family
    family_ipl <- IPL(family = fam, stabilization = "MAD", c = 0, zeta = 2)
    
    # Initial model for cross-validation
    model <- glmboostLSS(
      formula = list(
        mu     = y ~ .,  # V1-V1000
        sigma  = y ~ .,  # V1-V1000
        lambda = y ~ 1,  # intercept only
        alpha  = y ~ .   # V1-V1000
      ),
      data     = data_ipl,
      families = family_ipl,
      control  = boost_control(mstop = c(mstop_init, mstop_init, mstop_init, mstop_init),
                               nu = nu, trace = TRUE)
    )
    print("check point 1")
    
    # Cross-validation
    cvr <- cvrisk(model, grid = grid, folds = cv(model.weights(model)),
                  papply = mclapply)
    print("check point 2")
    
    m_stop[, i] <- mstop(cvr)
    print("check point 3")
    
    # Final model with optimal stopping iterations
    model <- glmboostLSS(
      formula = list(
        mu     = y ~ .,
        sigma  = y ~ .,
        lambda = y ~ 1,
        alpha  = y ~ .
      ),
      data     = data_ipl,
      families = family_ipl,
      control  = boost_control(mstop = c(m_stop[1, i], m_stop[2, i],
                                         m_stop[3, i], m_stop[4, i]),
                               nu = nu, trace = FALSE)
    )
    
    # Extract coefficients. off2int = TRUE folds each submodel's offset into its
    # intercept (mboost coef.glmboost, cf[[1]] <- cf[[1]] + offset), so the stored
    # intercept is the fitted value predict() uses, not just the boosted increment.
    # Slopes are unchanged. This matters most for lambda, whose offset is log(0.1).
    coefi <- coef(model, off2int = TRUE)
    estimativas_mu[, i]     <- coef(model, which = "", off2int = TRUE)$mu
    estimativas_sigma[, i]  <- coef(model, which = "", off2int = TRUE)$sigma
    estimativas_lambda[, i] <- coef(model, which = "", off2int = TRUE)$lambda
    estimativas_alpha[, i]  <- coef(model, which = "", off2int = TRUE)$alpha
    
    # ------------------------------------------------------------------
    # Selection frequency: count boosting iterations each var was picked
    # selected() returns base-learner indices; map to var_names_full
    # Index 1 = intercept, 2 = V1, ..., 1000 = V999, 1001 = V1000
    # ------------------------------------------------------------------
    sel_mu    <- table(selected(model$mu))
    sel_sigma <- table(selected(model$sigma))
    sel_alpha <- table(selected(model$alpha))
    
    # Convert numeric base-learner index to variable name
    idx_to_name <- function(sel_tbl) {
      idx <- as.integer(names(sel_tbl))
      setNames(as.integer(sel_tbl), var_names_full[idx])
    }
    
    named_mu    <- idx_to_name(sel_mu)
    named_sigma <- idx_to_name(sel_sigma)
    named_alpha <- idx_to_name(sel_alpha)
    
    # Accumulate into matrices
    matched_mu    <- intersect(names(named_mu),    var_names_full)
    matched_sigma <- intersect(names(named_sigma), var_names_full)
    matched_alpha <- intersect(names(named_alpha), var_names_full)
    
    freq_mu_mat[matched_mu,       i] <- named_mu[matched_mu]
    freq_sigma_mat[matched_sigma, i] <- named_sigma[matched_sigma]
    freq_alpha_mat[matched_alpha, i] <- named_alpha[matched_alpha]
    
    nvar_mu[i]     <- length(coefi$mu)
    nvar_sigma[i]  <- length(coefi$sigma)
    nvar_lambda[i] <- length(coefi$lambda)
    nvar_alpha[i]  <- length(coefi$alpha)
    
    print(paste("i =", i))
    i = i + 1
    print(paste("k =", k))
    
  }, error = function(e){ cat("ERROR :", conditionMessage(e), "\n") })
}

# ==============================================================================
# STORE RESULTS + AUTO-SAVE
# ==============================================================================

sim <- list(
  est_mu          = estimativas_mu,
  est_sigma       = estimativas_sigma,
  est_lambda      = estimativas_lambda,
  est_alpha       = estimativas_alpha,
  freq_sel_mu     = freq_mu_mat,     # n_coef x M: selection counts per var per run
  freq_sel_sigma  = freq_sigma_mat,
  freq_sel_alpha  = freq_alpha_mat,
  var_mu          = nvar_mu,
  var_sigma       = nvar_sigma,
  var_lambda      = nvar_lambda,
  var_alpha       = nvar_alpha,
  m_stop          = m_stop,
  true_beta_mu    = beta_mu,
  true_beta_sigma = beta_sigma,
  true_beta_alpha = beta_alpha,
  true_lambda     = lambda
)

sim_name <- paste0(example, "_sim_IPL_M", M, "_mstop", mstop_init, "_", fam, ".RData")
# Save into the IPL_simulations/ subfolder, where summary_table_IPL.r and
# experimentations_plots_IPL.R read from.
save(sim, file = file.path(save_path, "IPL_simulations", sim_name))
cat("\nSaved:", file.path(save_path, "IPL_simulations", sim_name), "\n")
