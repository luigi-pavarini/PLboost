# Boosted PL regression

rm(list = ls())

# ==============================================================================
# SETUP
# ==============================================================================
# ---- setup -------------------------------------------------------------------
library(PLboost)
save_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # results go to $PLBOOST_OUT/PL_simulations/
dir.create(file.path(save_path, "PL_simulations"), showWarnings = FALSE, recursive = TRUE)

require("PLreg")
require("gamboostLSS")
require("gamlss")
require("faux")
require("parallel")

# Parallel cross-validation (mclapply reads options("mc.cores")). Each core forks
# a full copy of the p=1000 model, so 8 cores is memory-hungry: run with other
# apps closed. If R still gets killed mid-loop, drop to 4 or 2 — the checkpoint
# below lets a rerun continue from where it stopped either way.
options(mc.cores = 10)

# M          = 100
# n          = 800
# p          = 1000
# fam        = 'SN'
# mstop_init = 500
# m_stop_cv  = 800
# nu         = 0.01

M          = 100
n          = 800
p          = 1000
fam        = 'SLASH' #"NO", "TF", "LO", "PE", "Hyp", "SLASH"
mstop_init = 100
m_stop_cv  = 500
nu         = 0.1

# Grid for cross-validation
grid <- make.grid(max = c(mu = m_stop_cv, sigma = m_stop_cv, lambda = m_stop_cv),
                  min = mstop_init,
                  length.out = 2, dense_mu_grid = TRUE)

# ==============================================================================
# COVARIATES GENERATION
# ==============================================================================

names_v <- paste0(rep("V", p), 1:p)
X <- rnorm_multi(n, vars = p, 
                 mu = c(rep(0, p)),
                 sd = c(rep(1, p)),
                 r = 0.5,
                 varnames = names_v,
                 empirical = FALSE)

# ==============================================================================
# TRUE COEFFICIENT VALUES
# ==============================================================================

beta_mu    <- c(-0.99, 0.49, -1.17, 0.38, rep(0, p - 4))
beta_sigma <- c(0, 0, -0.41, -0.24, 0.81, 0.70, rep(0, p - 6))

# ==============================================================================
# LINEAR PREDICTORS
# ==============================================================================

X_mat     <- as.matrix(X)
eta_mu    <- 0.5 + X_mat %*% beta_mu
eta_sigma <- 0.0 + X_mat %*% beta_sigma

# Parameter transformations
mu <- plogis(eta_mu)
mu[mu >= 1] <- 1 - 1e-6
mu[mu <= 0] <- 1e-6

sigma  <- exp(eta_sigma)
lambda <- 1.2

# ==============================================================================
# STORAGE MATRICES
# ==============================================================================

var_names_full <- c("(Intercept)", names_v)
n_coef <- p + 1  # 1001

estimativas_mu     <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
estimativas_sigma  <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
estimativas_lambda <- matrix(0, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))

# Matrices to store selection counts: rows = (intercept + p vars), cols = M iterations
freq_mu_mat    <- matrix(0L, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))
freq_sigma_mat <- matrix(0L, nrow = n_coef, ncol = M, dimnames = list(var_names_full, NULL))

m_stop     <- matrix(0, nrow = 3, ncol = M)  # 3: mu, sigma, lambda
nvar_mu     <- NULL
nvar_sigma  <- NULL
nvar_lambda <- NULL

# ==============================================================================
# SIMULATION LOOP
# ==============================================================================

# ---- checkpoint: resume a run that was killed mid-loop (e.g. out of memory) ----
# After every successful replication the accumulators, the fixed design X, the
# derived mu/sigma, the counters i/k and the RNG state are saved. If the script
# is rerun and a checkpoint for this family exists, it picks up where it stopped
# instead of starting over. The file is deleted once all M replications finish.
ckpt_file <- file.path(save_path, "PL_simulations", paste0(".ckpt_", fam, ".RData"))
if (file.exists(ckpt_file)) {
  message("Resuming from checkpoint: ", ckpt_file)
  load(ckpt_file)
} else {
  i <- 1
  k <- 0
}
while(i <= M){
  tryCatch({
    k = k + 1  
    
    yc <- PLreg::rPL(n = n, mu = mu, sigma = sigma,
                     lambda = lambda, family = fam, zeta = 2)
    
    yc[yc >= 1] <- 1 - 1e-6
    yc[yc <= 0] <- 1e-6
    
    data <- cbind(yc, X)
    
    # Initial model for cross-validation
    model <- glmboostLSS(formula = list(mu     = yc ~ .,
                                        sigma  = yc ~ .,
                                        lambda = yc ~ 1),
                         families = PowerLogit(zeta = 2, family = fam),
                         data = data,
                         control = boost_control(mstop = c(mstop_init, mstop_init, mstop_init),
                                                 nu = nu, trace = TRUE),
                         center = TRUE)
    
    # Cross-validation
    cvr <- cvrisk(model, grid = grid, folds = cv(model.weights(model)),
                  papply = mclapply)
    m_stop[, i] <- mstop(cvr)
    
    # Final model with optimal stopping iterations
    model <- glmboostLSS(formula = list(mu     = yc ~ .,
                                        sigma  = yc ~ .,
                                        lambda = yc ~ 1),
                         families = PowerLogit(zeta = 2, family = fam),
                         data = data,
                         control = boost_control(mstop = c(m_stop[1, i],
                                                           m_stop[2, i],
                                                           m_stop[3, i]),
                                                 nu = nu, trace = FALSE),
                         center = TRUE)
    
    # off2int = TRUE folds each submodel's offset into its intercept (mboost
    # coef.glmboost), so stored intercepts are the fitted values predict() uses.
    # Slopes are unchanged; the lambda intercept becomes exp(offset + coef) below.
    coefi <- coef(model, off2int = TRUE)
    estimativas_mu[, i]     <- coef(model, which = "", off2int = TRUE)$mu
    estimativas_sigma[, i]  <- coef(model, which = "", off2int = TRUE)$sigma
    estimativas_lambda[, i] <- coef(model, which = "", off2int = TRUE)$lambda
    
    # Count how many boosting iterations each variable was selected in this run
    sel_mu    <- table(selected(model$mu))
    sel_sigma <- table(selected(model$sigma))
    
    names(sel_mu)    <- paste0(rep("V", length(names(sel_mu))),    as.numeric(names(sel_mu))    - 1)
    names(sel_sigma) <- paste0(rep("V", length(names(sel_sigma))), as.numeric(names(sel_sigma)) - 1)
    
    # Accumulate counts into matrix column i (matched by row name)
    matched_mu    <- intersect(names(sel_mu),    var_names_full)
    matched_sigma <- intersect(names(sel_sigma), var_names_full)
    freq_mu_mat[matched_mu,       i] <- as.integer(sel_mu[matched_mu])
    freq_sigma_mat[matched_sigma, i] <- as.integer(sel_sigma[matched_sigma])
    
    nvar_mu[i]     <- length(coefi$mu)
    nvar_sigma[i]  <- length(coefi$sigma)
    nvar_lambda[i] <- length(coefi$lambda)
    
    print(paste("i =", i))
    i = i + 1
    print(paste("k =", k))

    # checkpoint after each success so an OOM kill loses at most one replication
    save(estimativas_mu, estimativas_sigma, estimativas_lambda,
         freq_mu_mat, freq_sigma_mat, m_stop,
         nvar_mu, nvar_sigma, nvar_lambda,
         X, mu, sigma, lambda, i, k, .Random.seed,
         file = ckpt_file)

  }, error = function(e){ cat("ERROR :", conditionMessage(e), "\n") })
  gc()
}

# ==============================================================================
# STORE RESULTS + AUTO-SAVE
# ==============================================================================

sim <- list(
  est_mu          = estimativas_mu, 
  est_sigma       = estimativas_sigma,
  est_lambda      = estimativas_lambda,
  freq_sel_mu     = freq_mu_mat,    # (p+1) x M matrix: selection counts per variable per run
  freq_sel_sigma  = freq_sigma_mat,
  var_mu          = nvar_mu,
  var_sigma       = nvar_sigma,
  var_lambda      = nvar_lambda,
  m_stop          = m_stop,
  true_beta_mu    = beta_mu,
  true_beta_sigma = beta_sigma,
  true_lambda     = lambda
)

sim_name <- paste0("sim_PL_M", M, "_mstop", mstop_init, "_", fam, ".RData")
save(sim, file = file.path(save_path, "PL_simulations", sim_name))
cat("\nSaved:", file.path(save_path, "PL_simulations", sim_name), "\n")
