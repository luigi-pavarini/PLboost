# ==============================================================================
# inep_run_cv_10rep.r
# Power-Logit (PL) boosted regression on the INEP Conceito Enade (Continuo),
# MEDICINA / SP, CPC 2023 — 5-FOLD REPEATED CROSS-VALIDATION.
#
# Response is OPEN on (0,1): tx = Conceito Enade / 5 (strictly interior on
# MEDICINA/SP, no boundary obs, so no squeeze needed). No point mass, so we drop
# the old one-inflated
# apparatus (IPL c=1 / BEOI / inflation-stratified folds) and use:
#
#   * PowerLogit(): the 3-parameter PL family (mu, sigma, lambda), each
#     regressed on all predictors -> 3 * p coefficients. With p=19 this is
#     3p = 57 > n = 55, i.e. the high-dimensional regime the method targets
#     (and 3p = 57 > n_train ~ 44 inside every fold).
#   * Beta (BE) as the classical bounded-response competitor (mu, sigma).
#
#   Outer  : R_REPEATS x 5-fold partitions, stratified on tx quantiles.
#   Inner  : cvrisk() selects mstop on each fold's training data (nested CV).
#   Zeta   : best zeta per PL sub-family chosen by Upsilon-on-train.
#
# n(MEDICINA/SP) = 55 -> 5-fold gives ~11 test obs/fold and ~44 train obs/fold.
# ==============================================================================

library(data.table)
library(mboost)
library(gamboostLSS)
library(gamlss.dist)
library(PLreg)
library(EnvStats)
library(caret)
# --- resolucao de caminho relativa a este arquivo (nao depende de outras pastas) --
.script_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(dirname(normalizePath(o)))
  }
  a  <- commandArgs(trailingOnly = FALSE)
  sp <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(sp)) return(dirname(normalizePath(sp)))
  getwd()
}

app_path     <- Sys.getenv("PLBOOST_OUT", unset = ".")   # set PLBOOST_OUT (must contain essentials/ and data/)
results_path <- file.path(app_path,"cross_val", "results")

source(file.path(app_path, "essentials", "families_PL.R"))
source(file.path(app_path, "essentials", "cv_summary.R"))
load(file.path(app_path, "data", "df_medicina_SP.RData"))

cat("Loaded df_inep — rows:", nrow(df_inep), " cols:", ncol(df_inep), "\n")

# ── BOOSTING SETTINGS ─────────────────────────────────────────────────────────
mstop_init <- 500
m_stop_cv  <- 800
nu         <- 0.001

# ── PL FAMILY LIST + ZETA GRIDS ──────────────────────────────────────────────
pl_families <- list(
  NO    = 2,
  TF    = c(3, 4, 5, 6, 8, 10, 15, 20),
  LO    = 2,
  PE    = seq(0.5, 5, by = 0.5),
  SN    = seq(0.5, 5, by = 0.5),
  SLASH = seq(0.5, 5, by = 0.5),
  Hyp   = seq(0.5, 5, by = 0.5)
)

grid_pl <- make.grid(
  max           = c(mu = m_stop_cv, sigma = m_stop_cv, lambda = m_stop_cv),
  min           = mstop_init,
  length.out    = 2,
  dense_mu_grid = TRUE
)

grid_beta <- make.grid(
  max           = c(mu = m_stop_cv, sigma = m_stop_cv),
  min           = mstop_init,
  length.out    = 2,
  dense_mu_grid = TRUE
)

# ── UPSILON FUNCTIONS (whole sample: tx already open on (0,1)) ────────────────
compute_upsilon_pl <- function(y, mu, sigma, lambda, zeta, family) {
  n   <- length(y)
  fda <- sort(PLreg:::pPL(y,
                          mu     = mu,
                          sigma  = sigma,
                          lambda = lambda[1],   # lambda ~ 1 -> constant scalar
                          zeta   = zeta,
                          family = family))
  mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
}

compute_upsilon_beta <- function(y, mu, sigma) {
  n   <- length(y)
  fda <- sort(gamlss.dist::pBE(y, mu = mu, sigma = sigma))
  mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
}

# ==============================================================================
# CROSS-VALIDATION SETTINGS
# ==============================================================================
K_FOLDS   <- 5
R_REPEATS <- 10

# ── RUN IN THREE BATCHES — quit & reopen R between them to reset memory creep ─
#   Batch 1: START=1, END=3,  TAG="1"  -> saves *_1.RData
#   Batch 2: START=4, END=6,  TAG="2"  -> saves *_2.RData
#   Batch 3: START=7, END=10, TAG="3"  -> saves *_3.RData
#
# >>> TO SWITCH BATCH: uncomment ONE of the three blocks below (lines 100-110)
#     and comment the others. Batch 1 is active by default.

# REPEAT_START <- 1
# REPEAT_END   <- 3
# BATCH_TAG    <- "1"

# REPEAT_START <- 4
# REPEAT_END   <- 6
# BATCH_TAG    <- "2"

REPEAT_START <- 7
REPEAT_END   <- 10
BATCH_TAG    <- "3"

options(mc.cores = 8)

TOTAL_ENTRIES <- K_FOLDS * R_REPEATS
BATCH_ENTRIES <- (REPEAT_END - REPEAT_START + 1L) * K_FOLDS

cat("\n################################################################\n")
cat("# INEP MEDICINA/SP —", K_FOLDS, "-fold CV (repeated)\n")
cat("#   folds (k)   :", K_FOLDS, " -> test ~", round(nrow(df_inep)/K_FOLDS),
    "obs/fold, train ~", round(nrow(df_inep) - nrow(df_inep)/K_FOLDS), "obs/fold\n")
cat("#   repeats (R) :", R_REPEATS, "  (each = one full", K_FOLDS, "-fold partition)\n")
cat("#   metrics averaged over", TOTAL_ENTRIES, "fold(s)  <- CV result = mean of these\n")
cat("#   THIS BATCH  : repeats", REPEAT_START, "-", REPEAT_END,
    "= entries", (REPEAT_START - 1L) * K_FOLDS + 1L, "-",
    REPEAT_END * K_FOLDS, "(", BATCH_ENTRIES, "fits ), tag =", BATCH_TAG, "\n")
cat("################################################################\n")

# Stratify folds on tx quantiles so each fold keeps a comparable spread of the
# (open) response (caret bins a numeric outcome into quantile groups).
cat("Response tx — range [", round(min(df_inep$tx), 4), ",",
    round(max(df_inep$tx), 4), "]  mean", round(mean(df_inep$tx), 4), "\n")

fold_seeds <- 1206 + seq_len(R_REPEATS) # credit base seed 123; repeat r -> 123 + r

make_folds <- function(seed) {
  set.seed(seed)
  caret::createFolds(df_inep$tx, k = K_FOLDS, list = TRUE, returnTrain = FALSE)
}

# ── RESULTS CONTAINERS — one list per family, K*R slots ──────────────────────
n_slots <- K_FOLDS * R_REPEATS
slot_id <- function(r, f) (r - 1L) * K_FOLDS + f

results_containers <- list(
  NO    = vector("list", n_slots),
  TF    = vector("list", n_slots),
  LO    = vector("list", n_slots),
  PE    = vector("list", n_slots),
  SN    = vector("list", n_slots),
  SLASH = vector("list", n_slots),
  Hyp   = vector("list", n_slots),
  Beta  = vector("list", n_slots)
)

# ==============================================================================
# OUTER LOOP — repeats x folds
# ==============================================================================
for (r in REPEAT_START:REPEAT_END) {
  
  folds <- make_folds(fold_seeds[r])

  for (f in seq_len(K_FOLDS)) {

    test_idx  <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(df_inep)), test_idx)

    df_train <- df_inep[train_idx]
    df_test  <- df_inep[test_idx]

    sid <- slot_id(r, f)
    cat("\n========== entry", sid, "/", TOTAL_ENTRIES,
        "| repeat", r, "/", R_REPEATS,
        "| fold", f, "/", K_FOLDS,
        "| n_train", nrow(df_train),
        "n_test", nrow(df_test), "==========\n")

    # ── PL FAMILIES ─────────────────────────────────────────────────────────────
    for (fam_name in names(pl_families)) {

      zeta_grid            <- pl_families[[fam_name]]
      zeta_results         <- vector("list", length(zeta_grid))
      upsilon_zeta_profile <- numeric(length(zeta_grid))

      for (z in seq_along(zeta_grid)) {

        zeta_z <- zeta_grid[z]

        tryCatch({

          fam_z <- PowerLogit(family = fam_name, stabilization = "MAD", zeta = zeta_z)

          model_z <- glmboostLSS(
            formula = list(mu     = tx ~ .,
                           sigma  = tx ~ .,
                           lambda = tx ~ 1),
            data     = df_train,
            families = fam_z,
            control  = boost_control(
              mstop = c(mstop_init, mstop_init, mstop_init),
              nu = nu, trace = T)
          )

          # Inner CV to find best mstop (nested — this fold's TRAIN only).
          # Seed the inner folds so mstop selection is reproducible across runs.
          # cv() runs in the main process before mclapply forks, so this alone
          # pins it (no need for L'Ecuyer). One seed per fold slot -> every
          # family/zeta in the same fold is compared on the same inner folds.
          set.seed(50000L + slot_id(r, f))
          cvr_z <- cvrisk(model_z, grid = grid_pl,
                          folds = mboost::cv(model.weights(model_z)),
                          papply = mclapply)
          opt_z <- mstop(cvr_z)

          fam_z <- PowerLogit(family = fam_name, stabilization = "MAD", zeta = zeta_z)

          model_z <- glmboostLSS(
            formula = list(mu     = tx ~ .,
                           sigma  = tx ~ .,
                           lambda = tx ~ 1),
            data     = df_train,
            families = fam_z,
            control  = boost_control(
              mstop = c(opt_z["mu"], opt_z["sigma"], opt_z["lambda"]),
              nu = nu, trace = FALSE)
          )

          # Predict on TRAIN — Upsilon (whole sample; response open on (0,1))
          mu_train_z     <- as.numeric(predict(model_z$mu,     newdata = df_train, type = "response"))
          sigma_train_z  <- as.numeric(predict(model_z$sigma,  newdata = df_train, type = "response"))
          lambda_train_z <- as.numeric(predict(model_z$lambda, newdata = df_train, type = "response"))

          upsilon_z <- compute_upsilon_pl(
            y      = df_train$tx,
            mu     = mu_train_z,
            sigma  = sigma_train_z,
            lambda = lambda_train_z,
            zeta   = zeta_z,
            family = fam_name
          )
          upsilon_zeta_profile[z] <- upsilon_z

          # Predict on TEST fold — MSE, MAE
          mu_hat_z     <- as.numeric(predict(model_z$mu,     newdata = df_test, type = "response"))
          sigma_hat_z  <- as.numeric(predict(model_z$sigma,  newdata = df_test, type = "response"))
          lambda_hat_z <- as.numeric(predict(model_z$lambda, newdata = df_test, type = "response"))

          mse_z <- mean((df_test$tx - mu_hat_z)^2)
          mae_z <- mean(abs(df_test$tx - mu_hat_z))

          zeta_results[[z]] <- list(
            zeta         = zeta_z,
            upsilon      = upsilon_z,
            mse          = mse_z,
            mae          = mae_z,
            mstop_mu     = opt_z["mu"],
            mstop_sigma  = opt_z["sigma"],
            mstop_lambda = opt_z["lambda"],
            mu_hat       = mu_hat_z,
            sigma_hat    = sigma_hat_z,
            lambda_hat   = lambda_hat_z,
            coef_mu      = coef(model_z$mu,     which = ""),
            coef_sigma   = coef(model_z$sigma,  which = ""),
            coef_lambda  = coef(model_z$lambda, which = "")
          )

          cat("  r", r, "f", f, "| fam =", fam_name, "| zeta =", zeta_z,
              "| Upsilon =", round(upsilon_z, 5),
              "| MSE x100 =", round(mse_z * 100, 3),
              "| MAE =", round(mae_z, 4), "\n")

          rm(model_z, cvr_z, fam_z)
          gc()

        }, error = function(e) {
          cat("  r", r, "f", f, "| fam =", fam_name, "| zeta =", zeta_z,
              "| FAILED:", conditionMessage(e), "\n")
        })
      }

      valid <- which(!sapply(zeta_results, is.null))

      if (length(valid) == 0) {
        cat("  r", r, "f", f, "| fam =", fam_name,
            "| ALL ZETAS FAILED — recorded as NULL\n")
        next
      }

      best_z <- valid[which.min(upsilon_zeta_profile[valid])]
      best   <- zeta_results[[best_z]]

      cat(">>> r", r, "f", f, "| fam =", fam_name, "| zeta_best =", best$zeta,
          "| Upsilon_best =", round(best$upsilon, 5), "\n")

      results_containers[[fam_name]][[sid]] <- list(
        repeat_id            = r,
        fold_id              = f,
        slot                 = sid,
        zeta_best            = best$zeta,
        upsilon              = best$upsilon,
        upsilon_zeta_profile = setNames(upsilon_zeta_profile, zeta_grid),
        error_pl             = best$mse,
        mae_pl               = best$mae,
        mstop_mu             = best$mstop_mu,
        mstop_sigma          = best$mstop_sigma,
        mstop_lambda         = best$mstop_lambda,
        y_test               = df_test$tx,
        mu_hat               = best$mu_hat,
        sigma_hat            = best$sigma_hat,
        lambda_hat           = best$lambda_hat,
        coef_mu              = best$coef_mu,
        coef_sigma           = best$coef_sigma,
        coef_lambda          = best$coef_lambda,
        test_idx             = test_idx
      )

      rm(zeta_results)
      gc()
    }

    # ── BETA (classical bounded-response competitor) ────────────────────────────
    tryCatch({

      fam_beta <- as.families(gamlss.dist::BE(), stabilization = "MAD")

      model_beta <- glmboostLSS(
        formula = list(mu    = tx ~ .,
                       sigma = tx ~ .),
        data     = df_train,
        families = fam_beta,
        control  = boost_control(
          mstop = c(mstop_init, mstop_init),
          nu = nu, trace = T)
      )

      set.seed(60000L + slot_id(r, f))   # reproducible inner folds for the beta fit
      cvr_beta <- cvrisk(model_beta, grid = grid_beta,
                         folds = mboost::cv(model.weights(model_beta)),
                         papply = mclapply)
      opt_beta <- mstop(cvr_beta)

      fam_beta <- as.families(gamlss.dist::BE(), stabilization = "MAD")

      model_beta <- glmboostLSS(
        formula = list(mu    = tx ~ .,
                       sigma = tx ~ .),
        data     = df_train,
        families = fam_beta,
        control  = boost_control(
          mstop = c(opt_beta["mu"], opt_beta["sigma"]),
          nu = nu, trace = FALSE)
      )

      mu_hat_beta    <- as.numeric(predict(model_beta$mu,    newdata = df_test, type = "response"))
      sigma_hat_beta <- as.numeric(predict(model_beta$sigma, newdata = df_test, type = "response"))

      mse_beta <- mean((df_test$tx - mu_hat_beta)^2)
      mae_beta <- mean(abs(df_test$tx - mu_hat_beta))

      mu_train_beta    <- as.numeric(predict(model_beta$mu,    newdata = df_train, type = "response"))
      sigma_train_beta <- as.numeric(predict(model_beta$sigma, newdata = df_train, type = "response"))

      upsilon_beta <- compute_upsilon_beta(
        y     = df_train$tx,
        mu    = mu_train_beta,
        sigma = sigma_train_beta
      )

      cat("  r", r, "f", f, "| fam = Beta",
          "| Upsilon(train) =", round(upsilon_beta, 5),
          "| MSE x100 =", round(mse_beta * 100, 3),
          "| MAE =", round(mae_beta, 4), "\n")

      results_containers[["Beta"]][[sid]] <- list(
        repeat_id   = r,
        fold_id     = f,
        slot        = sid,
        upsilon     = upsilon_beta,
        error_beta  = mse_beta,
        mae_beta    = mae_beta,
        mstop_mu    = opt_beta["mu"],
        mstop_sigma = opt_beta["sigma"],
        y_test      = df_test$tx,
        mu_hat      = mu_hat_beta,
        sigma_hat   = sigma_hat_beta,
        coef_mu     = coef(model_beta$mu,    which = ""),
        coef_sigma  = coef(model_beta$sigma, which = ""),
        test_idx    = test_idx
      )

      rm(model_beta, cvr_beta, fam_beta)
      gc()

    }, error = function(e) {
      cat("  r", r, "f", f, "| fam = Beta | FAILED:", conditionMessage(e),
          "— recorded as NULL\n")
    })

    cat(">>> DONE entry", sid, "/", TOTAL_ENTRIES,
        "(repeat", r, "fold", f, ") — this batch:",
        sid - (REPEAT_START - 1L) * K_FOLDS, "/", BATCH_ENTRIES, "\n")
    gc()
  }
}

# ==============================================================================
# RESULTS SUMMARY (this batch's filled slots only)
# ==============================================================================
summarise_cv(
  results_containers = results_containers,
  family_order = c(names(pl_families), "Beta"),
  meta = list(app = "INEP MEDICINA/SP", k_folds = K_FOLDS, r_repeats = R_REPEATS,
              repeat_start = REPEAT_START, repeat_end = REPEAT_END,
              batch_entries = BATCH_ENTRIES, total_entries = TOTAL_ENTRIES,
              batch_tag = BATCH_TAG),
  out_file = file.path(results_path,
                       paste0("summary_medicina_SP_cv_", BATCH_TAG, ".md"))
)

# ── SAVE — one .RData per family, suffixed by batch tag ───────────────────────
for (fam_name in names(pl_families)) {
  sim_results <- results_containers[[fam_name]]
  save(sim_results, var_dict,
       file = file.path(results_path,
                        paste0("sim_results_medicina_SP_cv_", fam_name, "_", BATCH_TAG, ".RData")))
  cat("Saved:", paste0("sim_results_medicina_SP_cv_", fam_name, "_", BATCH_TAG, ".RData"), "\n")
}

sim_results <- results_containers[["Beta"]]
save(sim_results, var_dict,
     file = file.path(results_path,
                      paste0("sim_results_medicina_SP_cv_Beta_", BATCH_TAG, ".RData")))
cat("Saved:", paste0("sim_results_medicina_SP_cv_Beta_", BATCH_TAG, ".RData"), "\n")
