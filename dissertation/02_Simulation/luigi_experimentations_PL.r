# Boosted PL regression

rm(list = ls())

# ==============================================================================
# SETUP
# ==============================================================================

# Needs to be settled!!!
base_path  <- "..."
save_path  <- "..."

source(file.path(base_path, "essentials", "families_PL.R"))

require("PLreg")
require("gamboostLSS")
require("gamlss")
require("faux")

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

i = 1
k = 0
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
    
    coefi <- coef(model, off2int = TRUE)
    estimativas_mu[, i]     <- coef(model, which = "")$mu
    estimativas_sigma[, i]  <- coef(model, which = "")$sigma
    estimativas_lambda[, i] <- coef(model, which = "")$lambda
    
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
    
  }, error = function(e){ cat("ERROR :", conditionMessage(e), "\n") })
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
load(file.path(save_path, "PL_simulations", sim_name))

#load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop100_NO.RData"))
#load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop100_TF.RData"))
#load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop100_LO.RData"))
#load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop500_SN.RData"))
#load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop100_PE.RData"))
#load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop100_Hyp.RData"))
load(file.path(save_path, "PL_simulations", "sim_PL_M100_mstop100_SLASH.RData"))

# ==============================================================================
# VISUALIZATION
# Column order: (Intercept)[1] | V1[2]...V1000[1001]
# mu    true: V1[2], V2[3], V3[4], V4[5]           non-inf: [c(6:1001),]
# sigma true: V3[4], V4[5], V5[6], V6[7]           non-inf: [c(2:3, 8:1001),]
# ==============================================================================

# Stopping iterations
boxplot(sim$m_stop[1,], sim$m_stop[2,], sim$m_stop[3,],
        names = c("mu", "sigma", "lambda"),
        ylab = "Stopping iterations")

# ---- Mu coefficients ----
boxplot(sim$est_mu[2,], sim$est_mu[3,], sim$est_mu[4,],
        sim$est_mu[5,], c(sim$est_mu[c(6:1001),]),
        ylim = c(-1.2, 0.5), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[mu]),
        main = "" #PL Model - Mu Parameter"
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
# beta_sigma = c(0, 0, -0.41, -0.24, 0.81, 0.70, 0...) -> V3[4], V4[5], V5[6], V6[7] are true
boxplot(sim$est_sigma[4,], sim$est_sigma[5,], sim$est_sigma[6,],
        sim$est_sigma[7,], c(sim$est_sigma[c(2:3, 8:1001),]),
        ylim = c(-0.70, 1.00), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[sigma]),
        main = "" #"PL Model - Sigma Parameter"
        )
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, expression(X[3]))
axis(1, 2, expression(X[4]))
axis(1, 3, expression(X[5]))
axis(1, 4, expression(X[6]))
axis(1, 5, "non-inf.")
segments(x0 = 0.5, x1 = 1.5, y0 = -0.41, y1 = -0.41, lty = 2, col = "red4")
segments(x0 = 1.5, x1 = 2.5, y0 = -0.24, y1 = -0.24, lty = 2, col = "red4")
segments(x0 = 2.5, x1 = 3.5, y0 =  0.81, y1 =  0.81, lty = 2, col = "red4")
segments(x0 = 3.5, x1 = 4.5, y0 =  0.70, y1 =  0.70, lty = 2, col = "red4")
segments(x0 = 4.5, x1 = 5.5, y0 =  0,    y1 =  0,    lty = 2, col = "red4")

# ---- Lambda ----
boxplot(exp(sim$est_lambda[1,]), ylim = c(0, 10), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ lambda),
        main = "" #"PL Model - Lambda Parameter"
        )
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, "exp(intercept)")
segments(x0 = 0.5, x1 = 1.5, y0 = 1.2, y1 = 1.2, lty = 2, col = "red4")

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
                      sigma = expression("Avg. selection proportion across M runs" ~ eta[sigma])
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
       main = "",#paste0("PL Model - Top ", top_n, " selected (", param_label, ")"),
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
