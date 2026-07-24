# ==============================================================================
# credit_combine.R
# Combines results from 2.1 (reps 1-50, seed 123) and
#                         2.2 (reps 51-100, seed 1234)
# into single .RData files, then prints full performance summary
# ==============================================================================
app_path <- Sys.getenv("PLBOOST_OUT", unset = ".")  # set PLBOOST_OUT or run from the repo root
results_path <- file.path(app_path, "01_replications", "results")

all_families <- c("NO", "TF", "LO", "PE", "SN", "SLASH", "Hyp", "BEOI")

for (fam_name in all_families) {

  f1 <- file.path(results_path, paste0("sim_results_credit_VT_", fam_name, "_1.RData"))
  f2 <- file.path(results_path, paste0("sim_results_credit_VT_", fam_name, "_2.RData"))

  load(f1); res1 <- sim_results; dict1 <- var_dict
  load(f2); res2 <- sim_results

  sim_results <- c(res1, res2)

  save(sim_results, var_dict,
       file = file.path(results_path,
                        paste0("sim_results_credit_VT_", fam_name, ".RData")))

  cat("Combined and saved:", paste0("sim_results_credit_VT_", fam_name, ".RData"),
      "— total entries:", length(sim_results),
      "| successes:", sum(!sapply(sim_results, is.null)), "\n")
}

# ==============================================================================
# PERFORMANCE SUMMARY
# ==============================================================================

cat("\n================================================================\n")
cat("  CREDIT LGD APPLICATION (Vermont) — FULL PERFORMANCE SUMMARY\n")
cat("================================================================\n\n")

for (fam_name in all_families) {

  fpath <- file.path(results_path,
                     paste0("sim_results_credit_VT_", fam_name, ".RData"))

  if (!file.exists(fpath)) {
    cat("---", fam_name, "--- FILE NOT FOUND\n\n")
    next
  }

  load(fpath)

  success <- which(!sapply(sim_results, is.null))
  n_ok    <- length(success)
  n_total <- length(sim_results)
  pct     <- round(n_ok / n_total * 100, 1)

  cat("---", fam_name, "---\n")
  cat("  Successes    :", n_ok, "/", n_total,
      paste0("(", pct, "%)"),
      "  |  Failures:", n_total - n_ok, "\n")

  if (n_ok > 0) {
    err_field <- if (fam_name == "BEOI") "error_beoi" else "error_ipl"
    mae_field <- if (fam_name == "BEOI") "mae_beoi"   else "mae_ipl"

    errors <- sapply(sim_results[success], function(x) x[[err_field]]) * 100
    mae    <- sapply(sim_results[success], function(x) x[[mae_field]])
    ups    <- sapply(sim_results[success], function(x) x$upsilon)

    cat(sprintf("  MSE x100     :  Mean = %.3f   SD = %.3f   Min = %.3f   Max = %.3f\n",
                mean(errors), sd(errors), min(errors), max(errors)))
    cat(sprintf("  MAE          :  Mean = %.4f  SD = %.4f  Min = %.4f  Max = %.4f\n",
                mean(mae), sd(mae), min(mae), max(mae)))
    cat(sprintf("  Upsilon      :  Mean = %.4f  SD = %.4f  Min = %.4f  Max = %.4f\n",
                mean(ups), sd(ups), min(ups), max(ups)))
  }
  cat("\n")
}

cat("================================================================\n")
