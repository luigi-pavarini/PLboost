# ==============================================================================
# credit_application_mod_2.R
# IPL boosted regression on LGD — 100 train/test replications (80/20)
# Vermont filter applied, IPL-NO family with c = 1 (one-inflation at lgd = 1)
# tryCatch ensures the loop never breaks on a bad split
# ==============================================================================

library(data.table)
library(mboost)
library(gamboostLSS)

# Needs to be settled!!!
base_path <- "..."
source(file.path(base_path, "essentials", "families_IPL.R"))
load(file.path(base_path, "df_lgd_all.RData"))

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

mstop_init <- 1000
m_stop_cv  <- 1200
nu         <- 0.001

grid <- make.grid(
  max           = c(mu = m_stop_cv, sigma = m_stop_cv,
                    lambda = m_stop_cv, alpha = m_stop_cv),
  min           = mstop_init,
  length.out    = 2,
  dense_mu_grid = TRUE
)

# ==============================================================================
# 100 TRAIN/TEST REPLICATIONS (80/20)
# while loop + tryCatch: skips failed splits, always reaches M successes
# ==============================================================================

n      <- nrow(df_lgd)
n_test <- round(n * 0.20)
M      <- 100

sim_results <- vector("list", M)

set.seed(42)
i <- 1  # successful replication counter
k <- 0  # total attempt counter

while (i <= M) {
  tryCatch({
    k <- k + 1
    
    test_idx  <- sample(n, n_test)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    df_train <- df_lgd[train_idx]
    df_test  <- df_lgd[test_idx]
    
    # ── Initial model ─────────────────────────────────────────────────────────
    fam_r <- IPL(family = "NO", stabilization = "MAD", c = 1, zeta = 2)
    
    model_r <- glmboostLSS(
      formula = list(
        mu     = lgd ~ .,
        sigma  = lgd ~ .,
        lambda = lgd ~ 1,
        alpha  = lgd ~ .
      ),
      data     = df_train,
      families = fam_r,
      control  = boost_control(
        mstop = c(mstop_init, mstop_init, mstop_init, mstop_init),
        nu = nu, trace = TRUE
      )
    )
    
    # ── Cross-validation ──────────────────────────────────────────────────────
    cvr_r <- cvrisk(model_r, grid = grid,
                    folds = cv(model.weights(model_r)),
                    papply = mclapply)
    opt_r <- mstop(cvr_r)
    
    # ── Final model at optimal stopping ───────────────────────────────────────
    fam_r <- IPL(family = "NO", stabilization = "MAD", c = 1, zeta = 2)
    
    model_r <- glmboostLSS(
      formula = list(
        mu     = lgd ~ .,
        sigma  = lgd ~ .,
        lambda = lgd ~ 1,
        alpha  = lgd ~ .
      ),
      data     = df_train,
      families = fam_r,
      control  = boost_control(
        mstop = c(opt_r["mu"], opt_r["sigma"],
                  opt_r["lambda"], opt_r["alpha"]),
        nu = nu, trace = FALSE
      )
    )
    
    # ── Predictions on test set ───────────────────────────────────────────────
    mu_hat     <- predict(model_r$mu,     newdata = df_test, type = "response")
    sigma_hat  <- predict(model_r$sigma,  newdata = df_test, type = "response")
    lambda_hat <- predict(model_r$lambda, newdata = df_test, type = "response")
    alpha_hat  <- predict(model_r$alpha,  newdata = df_test, type = "response")
    
    error_r <- mean((df_test$lgd - mu_hat)^2)
    
    # ── Predicted alpha on training set (for inflation distribution plot) ─────
    alpha_hat_train <- predict(model_r$alpha, newdata = df_train,
                               type = "response")
    
    # ── Coefficients ──────────────────────────────────────────────────────────
    coef_mu     <- coef(model_r$mu,     which = "")
    coef_sigma  <- coef(model_r$sigma,  which = "")
    coef_lambda <- coef(model_r$lambda, which = "")
    coef_alpha  <- coef(model_r$alpha,  which = "")
    
    sim_results[[i]] <- list(
      replicate       = i,
      error_ipl       = error_r,
      mstop_mu        = opt_r["mu"],
      mstop_sigma     = opt_r["sigma"],
      mstop_lambda    = opt_r["lambda"],
      mstop_alpha     = opt_r["alpha"],
      y_test          = df_test$lgd,
      mu_hat          = as.numeric(mu_hat),
      sigma_hat       = as.numeric(sigma_hat),
      lambda_hat      = as.numeric(lambda_hat),
      alpha_hat       = as.numeric(alpha_hat),
      alpha_hat_train = as.numeric(alpha_hat_train),
      coef_mu         = coef_mu,
      coef_sigma      = coef_sigma,
      coef_lambda     = coef_lambda,
      coef_alpha      = coef_alpha,
      test_idx        = test_idx
    )
    
    print(paste("i =", i, "/ 100 done  (attempt k =", k, ")"))
    i <- i + 1
    
  }, error = function(e) {
    cat("ERROR on attempt k =", k, ":", conditionMessage(e), "\n")
  })
}

cat("\nCompleted", i - 1, "successful replications in", k, "total attempts.\n")
cat("Failed attempts:", k - (i - 1), "\n")

# ==============================================================================
# QUICK RESULTS SUMMARY
# ==============================================================================

errors_ipl <- sapply(sim_results, function(x) x$error_ipl)
cat("\nMean MSE x 100 :", round(mean(errors_ipl) * 100, 3), "\n")
cat("SD   MSE x 100 :", round(sd(errors_ipl)   * 100, 3), "\n")

# ==============================================================================
# SAVE
# ==============================================================================

save(sim_results, var_dict,
     file = file.path(base_path, "sim_results_credit_VT.RData"))

cat("\nResults saved to sim_results_credit_VT.RData\n")