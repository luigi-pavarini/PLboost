library(data.table)
library(mboost)
library(gamboostLSS)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else NULL
base_path <- if (!is.null(script_path)) dirname(normalizePath(script_path)) else getwd()
source(file.path(base_path, "essentials", "families_PL.R"))
load(file.path(base_path, "df_itagrade.RData"))

# ── NA COUNT PER COLUMN ───────────────────────────────────────────────────────
na_counts <- sort(colSums(is.na(df)), decreasing = TRUE)
na_counts <- na_counts[na_counts > 0]
print(na_counts)
nrow(df)
ncol(df)

# ── RENAME PREDICTORS TO V1, V2, V3 ... + BUILD DICTIONARY ───────────────────
pred_cols <- setdiff(names(df), "y")
new_names <- paste0("V", seq_along(pred_cols))

var_dict <- data.table(
  new_name  = new_names,
  orig_name = pred_cols
)

setnames(df, old = pred_cols, new = new_names)

mstop_init <- 200
m_stop_cv  <- 1500
nu         <- 0.001

grid <- make.grid(max    = c(mu = m_stop_cv, sigma = m_stop_cv,
                             lambda = m_stop_cv),
                  min    = mstop_init,
                  length.out = 2, dense_mu_grid = TRUE)

# ── PREDICTION — 100 TRAIN/TEST SPLITS (80/20) ───────────────────────────────

n      <- nrow(df)
n_test <- 35
n_rep  <- 100

sim_results <- vector("list", n_rep)

set.seed(42)
for (r in seq_len(n_rep)) {
  
  test_idx  <- sample(n, n_test)
  train_idx <- setdiff(seq_len(n), test_idx)
  
  df_train <- df[train_idx]
  df_test  <- df[test_idx]
  
  fam_r <- PowerLogit(family = "NO", stabilization = "MAD", zeta = 2)
  
  model_r <- glmboostLSS(
    formula = list(
      mu     = y ~ .,
      sigma  = y ~ .,
      lambda = y ~ 1
    ),
    data     = df_train,
    families = fam_r,
    control  = boost_control(mstop = c(mstop_init, mstop_init,
                                       mstop_init),
                             nu = nu, trace = T)
  )
  
  cvr_r <- cvrisk(model_r, grid = grid, folds = cv(model.weights(model_r)),
                  papply = mclapply)
  opt_r <- mstop(cvr_r)
  
  fam_r <- PowerLogit(family = "NO", stabilization = "MAD", zeta = 2)
  
  model_r <- glmboostLSS(
    formula = list(
      mu     = y ~ .,
      sigma  = y ~ .,
      lambda = y ~ 1
    ),
    data     = df_train,
    families = fam_r,
    control  = boost_control(mstop = c(opt_r["mu"], opt_r["sigma"],
                                       opt_r["lambda"]),
                             nu = nu, trace = F)
  )
  
  mu_hat    <- predict(model_r$mu,     newdata = df_test, type = "response")
  sigma_hat <- predict(model_r$sigma,  newdata = df_test, type = "response")
  lambda_hat<- predict(model_r$lambda, newdata = df_test, type = "response")
  
  error_r   <- mean((df_test$y - mu_hat)^2)
  
  # coefficients for each submodel (named vectors, intercept + selected vars)
  coef_mu     <- coef(model_r$mu,     which = "")
  coef_sigma  <- coef(model_r$sigma,  which = "")
  coef_lambda <- coef(model_r$lambda, which = "")
  
  sim_results[[r]] <- list(
    replicate   = r,
    error_pl    = error_r,
    mstop_mu    = opt_r["mu"],
    mstop_sigma = opt_r["sigma"],
    mstop_lambda= opt_r["lambda"],
    y_test      = df_test$y,
    mu_hat      = as.numeric(mu_hat),
    sigma_hat   = as.numeric(sigma_hat),
    lambda_hat  = as.numeric(lambda_hat),
    coef_mu     = coef_mu,
    coef_sigma  = coef_sigma,
    coef_lambda = coef_lambda,
    test_idx    = test_idx
  )
  
  print(paste("replicate", r, "/ 100 done"))
}

# ── RESULTS (Table 4 equivalent) ─────────────────────────────────────────────
errors_pl <- sapply(sim_results, function(x) x$error_pl)

mean(errors_pl) * 100
sd(errors_pl)   * 100

# ── SAVE ─────────────────────────────────────────────────────────────────────
save(sim_results, var_dict,
     file = file.path(base_path, "sim_results_itagrade.RData"))