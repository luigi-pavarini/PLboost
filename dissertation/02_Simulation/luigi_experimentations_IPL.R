# Boosted IPL regression - Adapted with X_alpha uniform control approach

rm(list = ls())

# ==============================================================================
# SETUP
# ==============================================================================

# Needs to be settled!!!
base_path  <- "..."
save_path  <- "..."

source(file.path(base_path, "essentials", "IPL_distributions.R"))
source(file.path(base_path, "essentials", "families_IPL.R"))

require("PLreg")
require("gamboostLSS")
require("mboost")
require("faux")

#example    = "high"
example    = "medium"
#example    = "low"

M          = 100
n          = 800
p_alpha    = 2       # Uniform covariates for alpha control
p          = 1000 - p_alpha   # Main covariates (high-dimensional)
fam        = 'NO' #"NO", "TF", "LO", "SN", "PE", "Hyp", "SLASH"
mstop_init = 1000
m_stop_cv  = 1500
nu         = 0.01

# Grid for cross-validation (adapted for IPL with 4 parameters)
grid <- make.grid(max = c(mu = m_stop_cv, sigma = m_stop_cv, lambda = m_stop_cv, alpha = m_stop_cv),
                  min = mstop_init,
                  length.out = 2, dense_mu_grid = TRUE)

# ==============================================================================
# COVARIATES GENERATION - New approach with X_alpha
# ==============================================================================

# Main covariates (normal) for all parameters
names_v <- paste0(rep("V", p), 1:p)
X <- rnorm_multi(n, vars = p, 
                 mu = rep(0, p),
                 sd = rep(1, p),
                 r = 0.5,
                 varnames = names_v,
                 empirical = FALSE)

# SEPARATE uniform covariates for alpha control
X_alpha <- matrix(runif(n * p_alpha), nrow = n)
colnames(X_alpha) <- paste0("U", 1:p_alpha)

# ==============================================================================
# TRUE COEFFICIENT VALUES - Adapted structure
# ==============================================================================

# Coefficients for main parameters (using 998 covariates)
beta_mu    <- c(-0.99, 0.49, -1.17, 0.38, rep(0, p - 4))
beta_sigma <- c(0, 0, -0.41, -0.24, 0.81, rep(0, p - 5))

# SCENARIO 1 (Low):
# beta_alpha <- c(0.7, -0.5)
# eta_alpha = -2.5 + 0.7*U1 + (-0.5)*U2
# Range: [-2.5+0+(-0.5), -2.5+0.7+0] = [-3.0, -1.8] -> alpha in [0.047, 0.142]

# SCENARIO 2 (Medium):
beta_alpha <- c(0.70, -0.3)  # Extends proven second scenario
# eta_alpha = -1.5 + 0.70*U1 + (-0.3)*U2  
# Range: [-1.5+0-0.3, -1.5+0.7+0] = [-1.8, -0.8] -> alpha in [0.142, 0.310]

# SCENARIO 2.1 (Medium):
#beta_alpha <- c(0.6, -0.4)  # Extends proven second scenario
# eta_alpha = -1.4 + 0.6*U1 + (-0.4)*U2  
# Range: [-1.4+0-0.4, -1.4+0.6+0] = [-1.8, -0.8] -> alpha in [0.142, 0.310]

# SCENARIO 3 (High):
# beta_alpha <- c(0.4, -0.4)
# eta_alpha = -0.8 + 0.4*U1 + (-0.4)*U2
# Range: [-0.8+0+(-0.4), -0.8+0.4+0] = [-1.2, -0.4] -> alpha in [0.231, 0.401]

# SCENARIO 3.1 (High):
# beta_alpha <- c(0.4, -0.6)
# eta_alpha = -0.6 + 0.4*U1 + (-0.6)*U2
# Range: [-0.6+0+(-0.6), -0.6+0.4+0] = [-1.2, -0.2] -> alpha in [0.231, 0.450]

# ==============================================================================
# LINEAR PREDICTORS
# ==============================================================================

X_mat     <- as.matrix(X)
eta_mu    <- 0.5 + X_mat %*% beta_mu
eta_sigma <- 0.0 + X_mat %*% beta_sigma

#eta_alpha <- -2.5 + X_alpha %*% beta_alpha  
eta_alpha <- -1.5 + X_alpha %*% beta_alpha
#eta_alpha <- -0.8 + X_alpha %*% beta_alpha  # Proven intercept + uniform control

# Parameter transformations
mu     <- plogis(eta_mu)
mu     <- pmax(pmin(mu, 1 - 1e-6), 1e-6)
sigma  <- exp(eta_sigma)
lambda <- 1.2
alpha  <- plogis(eta_alpha)

# ==============================================================================
# STORAGE MATRICES
# Column order in data_ipl: y, U1, U2, V1...V998
# => coef() vector order: (Intercept)[1], U1[2], U2[3], V1[4]...V998[1001]
# ==============================================================================

var_names_full <- c("(Intercept)", paste0("U", 1:p_alpha), names_v)
n_coef <- p + p_alpha + 1  # 1001

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
    
    # Combine covariates: uniform first, then normal
    # Column order: y | U1 U2 | V1...V998
    data_ipl <- data.frame(y, X_alpha, X)
    
    # IPL family
    family_ipl <- IPL(family = fam, stabilization = "MAD", c = 0, zeta = 2)
    
    # Initial model for cross-validation
    model <- glmboostLSS(
      formula = list(
        mu     = y ~ .,  # V1-V998 + U1-U2
        sigma  = y ~ .,  # V1-V998 + U1-U2
        lambda = y ~ 1,  # intercept only
        alpha  = y ~ .   # V1-V998 + U1-U2
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
    
    # Extract coefficients
    coefi <- coef(model, off2int = TRUE)
    estimativas_mu[, i]     <- coef(model, which = "")$mu
    estimativas_sigma[, i]  <- coef(model, which = "")$sigma
    estimativas_lambda[, i] <- coef(model, which = "")$lambda
    estimativas_alpha[, i]  <- coef(model, which = "")$alpha
    
    # ------------------------------------------------------------------
    # Selection frequency: count boosting iterations each var was picked
    # selected() returns base-learner indices; map to var_names_full
    # Index 1 = intercept, 2 = U1, 3 = U2, 4 = V1, ..., 1001 = V998
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
save(sim, file = file.path(save_path, "IPL_simulations", sim_name))
load(file.path(save_path, "IPL_simulations", sim_name))

#load(file.path(save_path, "IPL_simulations", "high_sim_IPL_M100_mstop500_NO.RData"))
#load(file.path(save_path, "IPL_simulations", "low_sim_IPL_M100_mstop500_NO.RData"))
load(file.path(save_path, "IPL_simulations", "medium_sim_IPL_M100_mstop500_NO.RData"))

# ==============================================================================
# VISUALIZATION
# Column order reminder: (Intercept)[1] | U1[2] U2[3] | V1[4]...V998[1001]
# mu    true: V1[4], V2[5], V3[6], V4[7]         non-inf: [c(2:3, 8:1001),]
# sigma true: V3[6], V4[7], V5[8]                non-inf: [c(2:5, 9:1001),]
# alpha true: U1[2], U2[3]                        non-inf: [c(4:1001),]
# ==============================================================================

# Stopping iterations
boxplot(sim$m_stop[1,], sim$m_stop[2,], sim$m_stop[3,], sim$m_stop[4,],
        names = c("mu", "sigma", "lambda", "alpha"),
        ylab = "Stopping iterations")

# ---- Mu coefficients ----
boxplot(sim$est_mu[4,], sim$est_mu[5,], sim$est_mu[6,], sim$est_mu[7,],
        c(sim$est_mu[c(2:3, 8:1001),]),
        ylim = c(-1.3, 0.6), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[mu]),
        main = "" #"IPL Model - Mu Parameter"
        )
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, expression(X[1]))
axis(1, 2, expression(X[2]))
axis(1, 3, expression(X[3]))
axis(1, 4, expression(X[4]))
axis(1, 5, "non-inf.")
segments(x0 = 0.5, x1 = 1.5, y0 = -0.99, y1 = -0.99, lty = 2, col = "red4")
segments(x0 = 1.5, x1 = 2.5, y0 =  0.49, y1 =  0.49, lty = 2, col = "red4")
segments(x0 = 2.5, x1 = 3.5, y0 = -1.17, y1 = -1.17, lty = 2, col = "red4")
segments(x0 = 3.5, x1 = 4.5, y0 =  0.38, y1 =  0.38, lty = 2, col = "red4")
segments(x0 = 4.5, x1 = 5.5, y0 =  0,    y1 =  0,    lty = 2, col = "red4")

# ---- Sigma coefficients ----
# beta_sigma = c(0, 0, -0.41, -0.24, 0.81, 0...) -> V3[6], V4[7], V5[8] are true
boxplot(sim$est_sigma[6,], sim$est_sigma[7,], sim$est_sigma[8,],
        c(sim$est_sigma[c(2:5, 9:1001),]),
        ylim = c(-0.6, 1), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[sigma]),
        main = "" #"IPL Model - Sigma Parameter"
        )
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, expression(X[3]))
axis(1, 2, expression(X[4]))
axis(1, 3, expression(X[5]))
axis(1, 4, "non-inf.")
segments(x0 = 0.5, x1 = 1.5, y0 = -0.41, y1 = -0.41, lty = 2, col = "red4")
segments(x0 = 1.5, x1 = 2.5, y0 = -0.24, y1 = -0.24, lty = 2, col = "red4")
segments(x0 = 2.5, x1 = 3.5, y0 =  0.81, y1 =  0.81, lty = 2, col = "red4")
segments(x0 = 3.5, x1 = 4.5, y0 =  0,    y1 =  0,    lty = 2, col = "red4")

# ---- Lambda ----
boxplot(exp(sim$est_lambda[1,]), ylim = c(0, 3), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ lambda),
        main = "" #"IPL Model - Lambda Parameter"
        )
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, "exp(intercept)")
segments(x0 = 0.5, x1 = 1.5, y0 = 1.2, y1 = 1.2, lty = 2, col = "red4")

# ---- Alpha coefficients ----
# beta_alpha = c(0.70, -0.3) -> U1[2], U2[3] are true; non-inf = [c(4:1001),]
boxplot(sim$est_alpha[2,], sim$est_alpha[3,],
        c(sim$est_alpha[c(4:1001),]),
        ylim = c(-1, 2), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[alpha]),
        main = "" #"IPL Model - Alpha Parameter"
        )
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, expression(U[1]))
axis(1, 2, expression(U[2]))
axis(1, 3, "non-inf.")
segments(x0 = 0.5, x1 = 1.5, y0 = beta_alpha[1], y1 = beta_alpha[1], lty = 2, col = "red4")
segments(x0 = 1.5, x1 = 2.5, y0 = beta_alpha[2], y1 = beta_alpha[2], lty = 2, col = "red4")
segments(x0 = 2.5, x1 = 3.5, y0 = 0,             y1 = 0,             lty = 2, col = "red4")

# ==============================================================================
# LOLLIPOP PLOTS
# - Each column (run) is normalized to proportions first (don't mix samples)
# - Then average proportion across M runs is computed
# - Error bars show +/- sd_mult * sd across runs  (default: 1)
# - Stems: solid gray; error bars: dashed steel blue
# - Top 15 variables shown individually; rest collapsed into "others"
# ==============================================================================

lollipop_sel <- function(freq_mat, top_n = 15, param_label = "mu", sd_mult = 1) {
  
  # Step 1: normalize each run (column) to proportions — avoids mixing samples
  col_totals <- colSums(freq_mat)
  col_totals[col_totals == 0] <- 1  # guard against zero-sum columns
  prop_mat <- sweep(freq_mat, 2, col_totals, FUN = "/")
  
  # Drop intercept row
  prop_mat <- prop_mat[-1, ]
  
  # Step 2: average proportion and sd across M runs
  avg_prop <- rowSums(prop_mat) / M
  sd_prop  <- apply(prop_mat, 1, sd)
  
  # Step 3: select top_n by average proportion
  ord       <- order(avg_prop, decreasing = TRUE)[seq_len(top_n)]
  top_avg   <- avg_prop[ord]
  top_sd    <- sd_prop[ord]
  top_names <- names(top_avg)
  
  # Step 4: collapse remaining into "others"
  others_avg <- sum(avg_prop[-ord])
  others_sd  <- sqrt(sum(sd_prop[-ord]^2))  # pooled sd for "others"
  
  p_plot   <- c(top_avg, others_avg)
  sd_plot  <- c(top_sd,  others_sd)
  lbl_plot <- c(top_names, "others")
  n_plot   <- top_n + 1
  
  # Step 5: axis label
  xlab_expr <- switch(param_label,
                      mu    = expression("Avg. selection proportion across M runs" ~ eta[mu]),
                      sigma = expression("Avg. selection proportion across M runs" ~ eta[sigma]),
                      alpha = expression("Avg. selection proportion across M runs" ~ eta[alpha])
  )
  
  # Step 6: plot
  old_par <- par(mar = c(4, 6, 3, 2))
  
  x_max <- max(p_plot + sd_mult * sd_plot, na.rm = TRUE) * 1.15
  
  plot(p_plot, seq_len(n_plot),
       ylim = c(n_plot, 1),
       xlim = c(0, x_max),
       yaxt = "n",
       xlab = xlab_expr,
       ylab = "",
       main = "", #paste0("IPL Model - Top ", top_n, " selected (", param_label, ")"),
       pch  = 19, col = "gray30", cex = 1.3)
  
  # Lollipop stems: light dashed gray
  segments(x0 = 0, x1 = p_plot,
           y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "gray80", lwd = 0.8, lty = 2)
  
  # Error bars: sd_mult * sd, solid steel blue
  x_lo <- pmax(p_plot - sd_mult * sd_plot, 0)
  x_hi <- p_plot + sd_mult * sd_plot
  cap  <- 0.18  # vertical cap half-height in plot units
  
  # Horizontal bar
  segments(x0 = x_lo, x1 = x_hi,
           y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "#4682B4", lwd = 1.4)
  # Left cap
  segments(x0 = x_lo, x1 = x_lo,
           y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap,
           col = "#4682B4", lwd = 1.4)
  # Right cap
  segments(x0 = x_hi, x1 = x_hi,
           y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap,
           col = "#4682B4", lwd = 1.4)
  
  axis(2, at = seq_len(n_plot), labels = lbl_plot, las = 1, cex.axis = 0.85)
  abline(v = 0, col = "gray60", lty = 2)
  par(old_par)
}

lollipop_sel(sim$freq_sel_mu,    top_n = 15, param_label = "mu",    sd_mult = 1)
lollipop_sel(sim$freq_sel_sigma, top_n = 15, param_label = "sigma", sd_mult = 1)
lollipop_sel(sim$freq_sel_alpha, top_n = 15, param_label = "alpha", sd_mult = 1)