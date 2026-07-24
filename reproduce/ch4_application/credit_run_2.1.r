# ==============================================================================
# credit_application_mod_2.1.R
# IPL boosted regression on LGD — replications 1-50 (80/20)
# Vermont filter applied — all IPL families + BEOI benchmark
# tryCatch per family: failures recorded as NULL, loop always reaches 100
# ==============================================================================

library(data.table)
library(mboost)
library(PLboost)
library(gamboostLSS)
library(gamlss.dist)
library(EnvStats)
app_path <- Sys.getenv("PLBOOST_OUT", unset = ".")  # set PLBOOST_OUT or run from the repo root
results_path <- file.path(app_path, "01_replications", "results")

source(file.path(app_path, "essentials", "families_IPL.R"))
load(file.path(app_path, "data", "df_lgd_all.RData"))

# ==============================================================================
# STATE FILTER — Vermont only
# ==============================================================================

addr_state_cols <- grep("^addr_state_", names(df_lgd), value = TRUE)
chosen_state    <- "VT"
dummy_col       <- paste0("addr_state_", chosen_state)

if (!dummy_col %in% names(df_lgd))
  stop("State '", chosen_state, "' not found.")

df_lgd <- df_lgd[get(dummy_col) == 1]
df_lgd[, (addr_state_cols) := NULL]
cat("Filtered to", chosen_state, "— rows remaining:", nrow(df_lgd), "\n")

# ==============================================================================
# NA CHECK
# ==============================================================================

na_counts <- sort(colSums(is.na(df_lgd)), decreasing = TRUE)
na_counts <- na_counts[na_counts > 0]
print(na_counts)
cat("Rows:", nrow(df_lgd), "  Cols:", ncol(df_lgd), "\n")

# ==============================================================================
# RENAME PREDICTORS TO V1, V2, ... + BUILD DICTIONARY
# ==============================================================================

pred_cols <- setdiff(names(df_lgd), "lgd")
new_names <- paste0("V", seq_along(pred_cols))

var_dict <- data.table(
  new_name  = new_names,
  orig_name = pred_cols
)

setnames(df_lgd, old = pred_cols, new = new_names)

# ==============================================================================
# BOOSTING SETTINGS
# ==============================================================================

mstop_init <- 500 #1000
m_stop_cv  <- 800 #1200
nu         <- 0.001

# ==============================================================================
# IPL FAMILY LIST + ZETA GRIDS
# ==============================================================================

ipl_families <- list(
  NO    = 2,
  TF    = c(3, 4, 5, 6, 8, 10, 15, 20),
  LO    = 2,
  PE    = seq(0.5, 5, by = 0.5),
  SN    = seq(0.5, 5, by = 0.5),
  SLASH = seq(0.5, 5, by = 0.5),
  Hyp   = seq(0.5, 5, by = 0.5)
)

grid_ipl <- make.grid(
  max           = c(mu = m_stop_cv, sigma = m_stop_cv,
                    lambda = m_stop_cv, alpha = m_stop_cv),
  min           = mstop_init,
  length.out    = 2,
  dense_mu_grid = TRUE
)

grid_beoi <- make.grid(
  max           = c(mu = m_stop_cv, sigma = m_stop_cv, nu = m_stop_cv),
  min           = mstop_init,
  length.out    = 2,
  dense_mu_grid = TRUE
)

# ==============================================================================
# UPSILON FUNCTIONS
# ==============================================================================

# IPL: Upsilon computed on CONTINUOUS part only (lgd > 0 & lgd < 1)
# Uses pPL on the continuous observations — ignores the inflated point (lgd = 1)
compute_upsilon_ipl <- function(y, mu, sigma, lambda, zeta, family) {
  cont  <- y > 0 & y < 1
  y_c   <- y[cont]
  n     <- length(y_c)
  fda   <- sort(pPL(y_c,
                    mu     = mu[cont],
                    sigma  = sigma[cont],
                    lambda = lambda,
                    zeta   = zeta,
                    family = family))
  mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
}

# BEOI: Upsilon computed on CONTINUOUS part only (lgd > 0 & lgd < 1)
# Uses gamlss.dist::pBEINF0 or pBE depending on which parameterisation fits
compute_upsilon_beoi <- function(y, mu, sigma) {
  cont <- y > 0 & y < 1
  y_c  <- y[cont]
  n    <- length(y_c)
  fda <- sort(gamlss.dist::pBE(y_c,
                               mu    = mu[cont],
                               sigma = 1 / sqrt(sigma[cont] + 1)))
  mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
}

# ==============================================================================
# RESULTS CONTAINERS — one list per family, 100 entries each
# NULL entries mean that family failed on that replication
# ==============================================================================

options(mc.cores = 8)

results_containers <- list(
  NO    = vector("list", 50),
  TF    = vector("list", 50),
  LO    = vector("list", 50),
  PE    = vector("list", 50),
  SN    = vector("list", 50),
  SLASH = vector("list", 50),
  Hyp   = vector("list", 50),
  BEOI  = vector("list", 50)
)

# ==============================================================================
# 100 REPLICATIONS (for loop — always reaches 100)
# ==============================================================================

n      <- nrow(df_lgd)
n_test <- round(n * 0.20)

set.seed(123)
for (i in seq_len(50)) {
  
  test_idx  <- sample(n, n_test)
  train_idx <- setdiff(seq_len(n), test_idx)
  
  df_train <- df_lgd[train_idx]
  df_test  <- df_lgd[test_idx]
  
  cat("\n========== replication", i, "/ 100 ==========\n")
  
  # ── IPL FAMILIES ─────────────────────────────────────────────────────────────
  for (fam_name in names(ipl_families)) {
    
    zeta_grid            <- ipl_families[[fam_name]]
    zeta_results         <- vector("list", length(zeta_grid))
    upsilon_zeta_profile <- numeric(length(zeta_grid))
    
    for (z in seq_along(zeta_grid)) {
      
      zeta_z <- zeta_grid[z]
      
      tryCatch({
        
        # Initial fit
        fam_z <- IPL(family = fam_name, stabilization = "MAD", c = 1, zeta = zeta_z)
        
        model_z <- glmboostLSS(
          formula = list(mu     = lgd ~ .,
                         sigma  = lgd ~ .,
                         lambda = lgd ~ 1,
                         alpha  = lgd ~ .),
          data     = df_train,
          families = fam_z,
          control  = boost_control(
            mstop = c(mstop_init, mstop_init, mstop_init, mstop_init),
            nu = nu, trace = T)
        )
        
        # CV to find best mstop
        cvr_z <- cvrisk(model_z, grid = grid_ipl,
                        folds = mboost::cv(model.weights(model_z)),
                        papply = mclapply)
        opt_z <- mstop(cvr_z)
        
        # Refit at optimal mstop
        fam_z <- IPL(family = fam_name, stabilization = "MAD", c = 1, zeta = zeta_z)
        
        model_z <- glmboostLSS(
          formula = list(mu     = lgd ~ .,
                         sigma  = lgd ~ .,
                         lambda = lgd ~ 1,
                         alpha  = lgd ~ .),
          data     = df_train,
          families = fam_z,
          control  = boost_control(
            mstop = c(opt_z["mu"], opt_z["sigma"],
                      opt_z["lambda"], opt_z["alpha"]),
            nu = nu, trace = FALSE)
        )
        
        # Predict on TRAIN — Upsilon (continuous part only)
        mu_train_z     <- as.numeric(predict(model_z$mu,     newdata = df_train, type = "response"))
        sigma_train_z  <- as.numeric(predict(model_z$sigma,  newdata = df_train, type = "response"))
        lambda_train_z <- as.numeric(predict(model_z$lambda, newdata = df_train, type = "response"))[1]
        
        upsilon_z <- compute_upsilon_ipl(
          y      = df_train$lgd,
          mu     = mu_train_z,
          sigma  = sigma_train_z,
          lambda = lambda_train_z,
          zeta   = zeta_z,
          family = fam_name
        )
        upsilon_zeta_profile[z] <- upsilon_z
        
        # Predict on TEST — MSE, MAE
        mu_hat_z     <- as.numeric(predict(model_z$mu,     newdata = df_test, type = "response"))
        sigma_hat_z  <- as.numeric(predict(model_z$sigma,  newdata = df_test, type = "response"))
        lambda_hat_z <- as.numeric(predict(model_z$lambda, newdata = df_test, type = "response"))
        alpha_hat_z  <- as.numeric(predict(model_z$alpha,  newdata = df_test, type = "response"))
        
        mse_z <- mean((df_test$lgd - mu_hat_z)^2)
        mae_z <- mean(abs(df_test$lgd - mu_hat_z))
        
        zeta_results[[z]] <- list(
          zeta         = zeta_z,
          upsilon      = upsilon_z,
          mse          = mse_z,
          mae          = mae_z,
          mstop_mu     = opt_z["mu"],
          mstop_sigma  = opt_z["sigma"],
          mstop_lambda = opt_z["lambda"],
          mstop_alpha  = opt_z["alpha"],
          mu_hat       = mu_hat_z,
          sigma_hat    = sigma_hat_z,
          lambda_hat   = lambda_hat_z,
          alpha_hat    = alpha_hat_z,
          coef_mu      = coef(model_z$mu,     which = ""),
          coef_sigma   = coef(model_z$sigma,  which = ""),
          coef_lambda  = coef(model_z$lambda, which = ""),
          coef_alpha   = coef(model_z$alpha,  which = "")
        )
        
        cat("  rep", i, "| fam =", fam_name, "| zeta =", zeta_z,
            "| Upsilon =", round(upsilon_z, 5),
            "| MSE x100 =", round(mse_z * 100, 3),
            "| MAE =", round(mae_z, 4), "\n")

        rm(model_z, cvr_z, fam_z)
        gc()

      }, error = function(e) {
        cat("  rep", i, "| fam =", fam_name, "| zeta =", zeta_z,
            "| FAILED:", conditionMessage(e), "\n")
        # zeta_results[[z]] stays NULL
      })
    }
    
    # Select best zeta among successful ones
    valid <- which(!sapply(zeta_results, is.null))
    
    if (length(valid) == 0) {
      cat("  rep", i, "| fam =", fam_name, "| ALL ZETAS FAILED — recorded as NULL\n")
      results_containers[[fam_name]][[i]] <- NULL
      next
    }
    
    best_z <- valid[which.min(upsilon_zeta_profile[valid])]
    best   <- zeta_results[[best_z]]
    
    cat(">>> rep", i, "| fam =", fam_name, "| zeta_best =", best$zeta,
        "| Upsilon_best =", round(best$upsilon, 5), "\n")
    
    results_containers[[fam_name]][[i]] <- list(
      replicate            = i,
      zeta_best            = best$zeta,
      upsilon              = best$upsilon,
      upsilon_zeta_profile = setNames(upsilon_zeta_profile, zeta_grid),
      error_ipl            = best$mse,
      mae_ipl              = best$mae,
      mstop_mu             = best$mstop_mu,
      mstop_sigma          = best$mstop_sigma,
      mstop_lambda         = best$mstop_lambda,
      mstop_alpha          = best$mstop_alpha,
      y_test               = df_test$lgd,
      mu_hat               = best$mu_hat,
      sigma_hat            = best$sigma_hat,
      lambda_hat           = best$lambda_hat,
      alpha_hat            = best$alpha_hat,
      coef_mu              = best$coef_mu,
      coef_sigma           = best$coef_sigma,
      coef_lambda          = best$coef_lambda,
      coef_alpha           = best$coef_alpha,
      test_idx             = test_idx
    )

    rm(zeta_results)
    gc()
  }
  
  # ── BEOI ─────────────────────────────────────────────────────────────────────
  tryCatch({
    
    fam_beoi <- as.families(gamlss.dist::BEOI(), stabilization = "MAD")
    
    model_beoi <- glmboostLSS(
      formula = list(mu    = lgd ~ .,
                     sigma = lgd ~ .,
                     nu    = lgd ~ .),
      data     = df_train,
      families = fam_beoi,
      control  = boost_control(
        mstop = c(mstop_init, mstop_init, mstop_init),
        nu = nu, trace = T)
    )
    
    cvr_beoi <- cvrisk(model_beoi, grid = grid_beoi,
                       folds = mboost::cv(model.weights(model_beoi)),
                       papply = mclapply)
    opt_beoi <- mstop(cvr_beoi)
    
    fam_beoi <- as.families(gamlss.dist::BEOI(), stabilization = "MAD")
    
    model_beoi <- glmboostLSS(
      formula = list(mu    = lgd ~ .,
                     sigma = lgd ~ .,
                     nu    = lgd ~ .),
      data     = df_train,
      families = fam_beoi,
      control  = boost_control(
        mstop = c(opt_beoi["mu"], opt_beoi["sigma"], opt_beoi["nu"]),
        nu = nu, trace = FALSE)
    )
    
    mu_hat_beoi  <- as.numeric(predict(model_beoi$mu,    newdata = df_test, type = "response"))
    sigma_hat_beoi <- as.numeric(predict(model_beoi$sigma, newdata = df_test, type = "response"))
    nu_hat_beoi  <- as.numeric(predict(model_beoi$nu,    newdata = df_test, type = "response"))
    
    mse_beoi <- mean((df_test$lgd - mu_hat_beoi)^2)
    mae_beoi <- mean(abs(df_test$lgd - mu_hat_beoi))
    
    # Upsilon on TRAIN — continuous part only
    mu_train_beoi    <- as.numeric(predict(model_beoi$mu,    newdata = df_train, type = "response"))
    sigma_train_beoi <- as.numeric(predict(model_beoi$sigma, newdata = df_train, type = "response"))
    
    upsilon_beoi <- compute_upsilon_beoi(
      y     = df_train$lgd,
      mu    = mu_train_beoi,
      sigma = sigma_train_beoi
    )
    
    cat("  rep", i, "| fam = BEOI",
        "| Upsilon(train) =", round(upsilon_beoi, 5),
        "| MSE x100 =", round(mse_beoi * 100, 3),
        "| MAE =", round(mae_beoi, 4), "\n")
    
    results_containers[["BEOI"]][[i]] <- list(
      replicate   = i,
      upsilon     = upsilon_beoi,
      error_beoi  = mse_beoi,
      mae_beoi    = mae_beoi,
      mstop_mu    = opt_beoi["mu"],
      mstop_sigma = opt_beoi["sigma"],
      mstop_nu    = opt_beoi["nu"],
      y_test      = df_test$lgd,
      mu_hat      = mu_hat_beoi,
      sigma_hat   = sigma_hat_beoi,
      nu_hat      = nu_hat_beoi,
      coef_mu     = coef(model_beoi$mu,    which = ""),
      coef_sigma  = coef(model_beoi$sigma, which = ""),
      coef_nu     = coef(model_beoi$nu,    which = ""),
      test_idx    = test_idx
    )

    rm(model_beoi, cvr_beoi, fam_beoi)
    gc()

  }, error = function(e) {
    cat("  rep", i, "| fam = BEOI | FAILED:", conditionMessage(e),
        "— recorded as NULL\n")
  })

  print(paste("replication", i, "/ 50 done"))
  gc()
}

# ==============================================================================
# RESULTS SUMMARY
# ==============================================================================

cat("\n================================================================\n")
cat("  SUMMARY ACROSS 100 REPLICATIONS\n")
cat("================================================================\n\n")

for (fam_name in names(ipl_families)) {
  res     <- results_containers[[fam_name]]
  success <- which(!sapply(res, is.null))
  n_ok    <- length(success)
  cat("---", fam_name, "--- successes:", n_ok, "/ 50  (failures:", 50 - n_ok, ")\n")
  if (n_ok > 0) {
    errors <- sapply(res[success], function(x) x$error_ipl) * 100
    mae    <- sapply(res[success], function(x) x$mae_ipl)
    ups    <- sapply(res[success], function(x) x$upsilon)
    cat("  MSE x100  Mean:", round(mean(errors), 3), "  SD:", round(sd(errors), 3), "\n")
    cat("  MAE       Mean:", round(mean(mae),    4), "  SD:", round(sd(mae),    4), "\n")
    cat("  Upsilon   Mean:", round(mean(ups),    4), "  SD:", round(sd(ups),    4), "\n")
  }
}

res_beoi  <- results_containers[["BEOI"]]
success_b <- which(!sapply(res_beoi, is.null))
n_ok_b    <- length(success_b)
cat("--- BEOI --- successes:", n_ok_b, "/ 50  (failures:", 50 - n_ok_b, ")\n")
if (n_ok_b > 0) {
  errors_b <- sapply(res_beoi[success_b], function(x) x$error_beoi) * 100
  mae_b    <- sapply(res_beoi[success_b], function(x) x$mae_beoi)
  ups_b    <- sapply(res_beoi[success_b], function(x) x$upsilon)
  cat("  MSE x100  Mean:", round(mean(errors_b), 3), "  SD:", round(sd(errors_b), 3), "\n")
  cat("  MAE       Mean:", round(mean(mae_b),    4), "  SD:", round(sd(mae_b),    4), "\n")
  cat("  Upsilon   Mean:", round(mean(ups_b),    4), "  SD:", round(sd(ups_b),    4), "\n")
}

# ==============================================================================
# SAVE — one .RData per family
# ==============================================================================

for (fam_name in names(ipl_families)) {
  sim_results <- results_containers[[fam_name]]
  save(sim_results, var_dict,
       file = file.path(results_path,
                        paste0("sim_results_credit_VT_", fam_name, "_1.RData")))
  cat("Saved:", paste0("sim_results_credit_VT_", fam_name, "_1.RData"), "\n")
}

sim_results <- results_containers[["BEOI"]]
save(sim_results, var_dict,
     file = file.path(results_path, "sim_results_credit_VT_BEOI_1.RData"))
cat("Saved: sim_results_credit_VT_BEOI_1.RData\n")