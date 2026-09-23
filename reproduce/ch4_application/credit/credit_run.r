# ==============================================================================
# credit_run_cv.R
# IPL boosted regression on LGD (Vermont) — 5-FOLD CROSS-VALIDATION
# Design: one random partition into 5 groups; each group is the test
# set once; the CV result is the MEAN of the 5 fold metrics (RMSE, MAE, Upsilon).
#
#   Outer  : R_REPEATS x 5-fold partitions, STRATIFIED on the inflation indicator
#            (lgd == 1 vs lgd in (0,1)) so every fold keeps the ~60/40 mix.
#            Folds disjoint & exhaustive. DEFAULT R_REPEATS = 1 (single pass,
#            as specified for the dissertation). The repeat structure is left settable:
#            raise R_REPEATS only if repeated CV is agreed with him.
#   Inner  : cvrisk() still selects mstop on each fold's 4/5 training data
#            (nested CV — no leakage into the held-out fold).
#   Zeta   : best zeta per family chosen by Upsilon-on-train (continuous part).
#
# n(VT) = 320  -> 5-fold gives ~64 test obs/fold and ~256 train obs/fold.
#
# Memory discipline (from the IPL runs): options(mc.cores = 8) + rm()/gc().
# Batching by REPEAT range only matters if R_REPEATS > 1.
# ==============================================================================

library(data.table)
library(mboost)
library(gamboostLSS)
library(gamlss.dist)
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


app_path     <- Sys.getenv("PLBOOST_OUT", unset = ".")   # set PLBOOST_OUT (must contain data/)
results_path <- file.path(app_path,"cross_val", "results")

library(PLboost); library(PLreg); library(gamboostLSS)   # PL / IPL families and distributions
source(file.path(.script_dir(), "..", "essentials", "cv_summary.R"))
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

# ── NA CHECK ──────────────────────────────────────────────────────────────────
na_counts <- sort(colSums(is.na(df_lgd)), decreasing = TRUE)
na_counts <- na_counts[na_counts > 0]
print(na_counts)
cat("Rows:", nrow(df_lgd), "  Cols:", ncol(df_lgd), "\n")

# ── RENAME PREDICTORS TO V1, V2, ... + BUILD DICTIONARY ──────────────────────
pred_cols <- setdiff(names(df_lgd), "lgd")
new_names <- paste0("V", seq_along(pred_cols))

var_dict <- data.table(
  new_name  = new_names,
  orig_name = pred_cols
)

setnames(df_lgd, old = pred_cols, new = new_names)

# ── BOOSTING SETTINGS ─────────────────────────────────────────────────────────
mstop_init <- 500
m_stop_cv  <- 800
nu         <- 0.001

# ── IPL FAMILY LIST + ZETA GRIDS ──────────────────────────────────────────────
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

# ── UPSILON FUNCTIONS (continuous part only: lgd > 0 & lgd < 1) ───────────────
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
# CROSS-VALIDATION SETTINGS
# ==============================================================================
K_FOLDS      <- 5           # 5-fold (test ~64, train ~256 / fold)
R_REPEATS    <- 10          # REPEATED CV — 10 passes of the 5-fold partition.

# ── RUN IN TWO BATCHES — quit & reopen R between them to reset memory creep ──
#   Batch 1 (set below): START=1,  END=3,  TAG="1"  -> saves *_1.RData
#   Batch 2 (then swap): START=4,  END=6,  TAG="2"  -> saves *_2.RData
#   Batch 3 (set below): START=7,  END=9,  TAG="3"  -> saves *_3.RData
#   Batch 4 (then swap): START=10, END=10, TAG="4"  -> saves *_4.RData

REPEAT_START <- 1           # batch 1: 1   | batch 2: 1
REPEAT_END   <- 3           # batch 1: 5   | batch 2: 3
BATCH_TAG    <- "1"         # batch 1: "1" | batch 2: "2"

# REPEAT_START <- 4           # batch 1: 1   | batch 2: 4
# REPEAT_END   <- 6           # batch 1: 5   | batch 2: 6
# BATCH_TAG    <- "2"         # batch 1: "1" | batch 2: "2"
# 
# REPEAT_START <- 7           # batch 1: 1   | batch 3: 7
# REPEAT_END   <- 10           # batch 1: 5   | batch 3: 9
# BATCH_TAG    <- "3"         # batch 1: "1" | batch 3: 32"
# 
# REPEAT_START <- 10          # batch 1: 1   | batch 4: 10
# REPEAT_END   <- 10          # batch 1: 5   | batch 4: 10
# BATCH_TAG    <- "4"         # batch 1: "1" | batch 4: "4"

# IPL fits 4 parameters (mu, sigma, lambda, alpha) and credit-VT is the larger
# data (n = 320 vs itagrade's 171), so each cvrisk() worker is much heavier than
# in the PL run. 8 parallel mclapply forks OOM'd a 16 GB machine here. 4 halves
# the peak with NO effect on results (identical folds/seeds — only fewer workers
# run at once). Drop to 2 if it still struggles; raise back toward 8 on a
# bigger-RAM machine.  NB: 10 repeats x 5 folds is 10x the work of 1rep.
options(mc.cores = 8)

TOTAL_ENTRIES <- K_FOLDS * R_REPEATS                       # error bars computed over this
BATCH_ENTRIES <- (REPEAT_END - REPEAT_START + 1L) * K_FOLDS # entries produced by THIS batch

cat("\n################################################################\n")
cat("# credit VT —", K_FOLDS, "-fold CV", if (R_REPEATS > 1) "(repeated)" else "(single pass)", "\n")
cat("#   folds (k)   :", K_FOLDS, " -> test ~", round(nrow(df_lgd)/K_FOLDS),
    "obs/fold, train ~", round(nrow(df_lgd) - nrow(df_lgd)/K_FOLDS), "obs/fold\n")
cat("#   repeats (R) :", R_REPEATS, "  (each = one full", K_FOLDS, "-fold partition)\n")
cat("#   metrics averaged over", TOTAL_ENTRIES, "fold(s)  <- CV result = mean of these\n")
cat("#   THIS BATCH  : repeats", REPEAT_START, "-", REPEAT_END,
    "= entries", (REPEAT_START - 1L) * K_FOLDS + 1L, "-",
    REPEAT_END * K_FOLDS, "(", BATCH_ENTRIES, "fits ), tag =", BATCH_TAG, "\n")
cat("################################################################\n")

# Stratify folds on the inflation indicator so each fold keeps the ~60/40
# continuous/inflated mix (lgd == 1 is the inflated point mass).
strat_factor <- factor(ifelse(df_lgd$lgd >= 1, "infl", "cont"))
cat("Stratification — continuous:", sum(strat_factor == "cont"),
    " inflated(lgd=1):", sum(strat_factor == "infl"), "\n")

fold_seeds <-  1206 + seq_len(R_REPEATS)        # credit base seed 123; repeat r -> 123 + r

make_folds <- function(seed) {
  set.seed(seed)
  caret::createFolds(strat_factor, k = K_FOLDS, list = TRUE, returnTrain = FALSE)
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
  BEOI  = vector("list", n_slots)
)

# ==============================================================================
# OUTER LOOP — repeats x folds
# ==============================================================================
for (r in REPEAT_START:REPEAT_END) {

  folds <- make_folds(fold_seeds[r])

  for (f in seq_len(K_FOLDS)) {

    test_idx  <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(df_lgd)), test_idx)

    df_train <- df_lgd[train_idx]
    df_test  <- df_lgd[test_idx]

    sid <- slot_id(r, f)
    cat("\n========== entry", sid, "/", TOTAL_ENTRIES,
        "| repeat", r, "/", R_REPEATS,
        "| fold", f, "/", K_FOLDS,
        "| n_train", nrow(df_train),
        "n_test", nrow(df_test), "==========\n")

    # ── IPL FAMILIES ───────────────────────────────────────────────────────────
    for (fam_name in names(ipl_families)) {

      zeta_grid            <- ipl_families[[fam_name]]
      zeta_results         <- vector("list", length(zeta_grid))
      upsilon_zeta_profile <- numeric(length(zeta_grid))

      for (z in seq_along(zeta_grid)) {

        zeta_z <- zeta_grid[z]

        tryCatch({

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

          # Inner CV to find best mstop (nested — this fold's TRAIN only).
          # Seed the inner folds so mstop selection is reproducible across runs.
          # cv() runs in the main process before mclapply forks, so this alone
          # pins it (no need for L'Ecuyer). One seed per fold slot.
          set.seed(50000L + sid)
          cvr_z <- cvrisk(model_z, grid = grid_ipl,
                          folds = mboost::cv(model.weights(model_z)),
                          papply = mclapply)
          opt_z <- mstop(cvr_z)

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

          # Predict on TEST fold — MSE, MAE
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

    # ── BEOI ───────────────────────────────────────────────────────────────────
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

      set.seed(60000L + sid)   # reproducible inner folds for the BEOI fit
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

      mu_hat_beoi    <- as.numeric(predict(model_beoi$mu,    newdata = df_test, type = "response"))
      sigma_hat_beoi <- as.numeric(predict(model_beoi$sigma, newdata = df_test, type = "response"))
      nu_hat_beoi    <- as.numeric(predict(model_beoi$nu,    newdata = df_test, type = "response"))

      mse_beoi <- mean((df_test$lgd - mu_hat_beoi)^2)
      mae_beoi <- mean(abs(df_test$lgd - mu_hat_beoi))

      mu_train_beoi    <- as.numeric(predict(model_beoi$mu,    newdata = df_train, type = "response"))
      sigma_train_beoi <- as.numeric(predict(model_beoi$sigma, newdata = df_train, type = "response"))

      upsilon_beoi <- compute_upsilon_beoi(
        y     = df_train$lgd,
        mu    = mu_train_beoi,
        sigma = sigma_train_beoi
      )

      cat("  r", r, "f", f, "| fam = BEOI",
          "| Upsilon(train) =", round(upsilon_beoi, 5),
          "| MSE x100 =", round(mse_beoi * 100, 3),
          "| MAE =", round(mae_beoi, 4), "\n")

      results_containers[["BEOI"]][[sid]] <- list(
        repeat_id   = r,
        fold_id     = f,
        slot        = sid,
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
      cat("  r", r, "f", f, "| fam = BEOI | FAILED:", conditionMessage(e),
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
  family_order = c(names(ipl_families), "BEOI"),
  meta = list(app = "credit VT", k_folds = K_FOLDS, r_repeats = R_REPEATS,
              repeat_start = REPEAT_START, repeat_end = REPEAT_END,
              batch_entries = BATCH_ENTRIES, total_entries = TOTAL_ENTRIES,
              batch_tag = BATCH_TAG),
  out_file = file.path(results_path,
                       paste0("summary_credit_VT_cv_", BATCH_TAG, ".md"))
)

# ── SAVE — one .RData per family, suffixed by batch tag ───────────────────────
for (fam_name in names(ipl_families)) {
  sim_results <- results_containers[[fam_name]]
  save(sim_results, var_dict,
       file = file.path(results_path,
                        paste0("sim_results_credit_VT_cv_", fam_name, "_", BATCH_TAG, ".RData")))
  cat("Saved:", paste0("sim_results_credit_VT_cv_", fam_name, "_", BATCH_TAG, ".RData"), "\n")
}

sim_results <- results_containers[["BEOI"]]
save(sim_results, var_dict,
     file = file.path(results_path,
                      paste0("sim_results_credit_VT_cv_BEOI_", BATCH_TAG, ".RData")))
cat("Saved:", paste0("sim_results_credit_VT_cv_BEOI_", BATCH_TAG, ".RData"), "\n")
