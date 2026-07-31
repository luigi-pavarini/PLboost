# ==============================================================================
# cv_summary.R — summarise a finished CV run.
#
# Replaces the per-family summary block that used to be copy-pasted at the end of
# every *_run_cv.R. It does two things:
#   1. prints the same per-family Mean/SD lines to the console (interactive runs
#      look exactly as before);
#   2. writes a Markdown table next to the .RData results, so the summary no
#      longer has to be saved by hand (this is what the old "1rep results.rtf"
#      was — now produced automatically and in a paste-ready format).
#
# One helper serves both runs (PL = itagrade, IPL = credit). Metric fields are
# coalesced by name, so it does not matter whether an entry stores its squared
# error under error_pl / error_ipl / error_beta / error_beoi (same for MAE).
#
# Arguments
#   results_containers : named list of per-family result lists (NULL = failed fit)
#   family_order       : family names in display order, INCLUDING the extra
#                        family (e.g. c(names(pl_families), "Beta"))
#   meta               : list(app, k_folds, r_repeats, repeat_start, repeat_end,
#                             batch_entries, total_entries, batch_tag)
#   out_file           : full path of the .md summary to write
#
# Returns (invisibly) the summary as a data.frame.
# ==============================================================================

# first non-NULL of the named fields on one result entry (NA if none present)
.cv_first <- function(entry, fields) {
  for (nm in fields) if (!is.null(entry[[nm]])) return(entry[[nm]])
  NA_real_
}

# (mse x100, mae, upsilon) for one (repeat, fold) result entry
cv_metrics <- function(entry) {
  c(mse100  = .cv_first(entry, c("error_pl", "error_ipl",
                                 "error_beta", "error_beoi")) * 100,
    mae     = .cv_first(entry, c("mae_pl", "mae_ipl",
                                 "mae_beta", "mae_beoi")),
    upsilon = entry$upsilon)
}

summarise_cv <- function(results_containers, family_order, meta, out_file) {

  fmt <- function(x, d)
    if (is.na(x)) "—" else formatC(round(x, d), format = "f", digits = d)

  # ── aggregate per family ────────────────────────────────────────────────────
  rows <- lapply(family_order, function(fam) {
    res     <- results_containers[[fam]]
    success <- which(!vapply(res, is.null, logical(1)))
    n_ok    <- length(success)
    base    <- data.frame(Family = fam,
                          Success = paste0(n_ok, "/", meta$batch_entries),
                          MSE100_mean = NA_real_, MSE100_sd = NA_real_,
                          MAE_mean = NA_real_, MAE_sd = NA_real_,
                          Upsilon_mean = NA_real_, Upsilon_sd = NA_real_,
                          stringsAsFactors = FALSE)
    if (n_ok == 0) return(base)
    M <- vapply(res[success], cv_metrics, numeric(3))   # 3 x n_ok
    base$MSE100_mean  <- mean(M["mse100", ]);  base$MSE100_sd  <- sd(M["mse100", ])
    base$MAE_mean     <- mean(M["mae", ]);     base$MAE_sd     <- sd(M["mae", ])
    base$Upsilon_mean <- mean(M["upsilon", ]); base$Upsilon_sd <- sd(M["upsilon", ])
    base
  })
  tab <- do.call(rbind, rows)

  # ── console (same information as the old inline block) ───────────────────────
  cat("\n================================================================\n")
  cat("  SUMMARY —", meta$app, meta$k_folds, "-fold CV\n")
  cat("  THIS BATCH: repeats", meta$repeat_start, "to", meta$repeat_end,
      "=", meta$batch_entries, "entries (of", meta$total_entries, "total)\n")
  cat("================================================================\n\n")
  for (i in seq_len(nrow(tab))) {
    cat("---", tab$Family[i], "--- successes:", tab$Success[i], "(this batch)\n")
    if (!is.na(tab$MSE100_mean[i])) {
      cat("  MSE x100  Mean:", fmt(tab$MSE100_mean[i], 3),
          "  SD:", fmt(tab$MSE100_sd[i], 3), "\n")
      cat("  MAE       Mean:", fmt(tab$MAE_mean[i], 4),
          "  SD:", fmt(tab$MAE_sd[i], 4), "\n")
      cat("  Upsilon   Mean:", fmt(tab$Upsilon_mean[i], 4),
          "  SD:", fmt(tab$Upsilon_sd[i], 4), "\n")
    }
  }

  # ── Markdown file ───────────────────────────────────────────────────────────
  run_desc <- if (meta$r_repeats > 1)
    paste0("repeated CV, R = ", meta$r_repeats) else "single pass (R = 1)"
  md <- c(
    paste0("# ", meta$app, " — ", meta$k_folds, "-fold CV summary"),
    "",
    paste0("- Run: ", run_desc),
    paste0("- Batch tag: ", meta$batch_tag,
           "  (repeats ", meta$repeat_start, "–", meta$repeat_end, ")"),
    paste0("- Entries: ", meta$batch_entries, " this batch / ",
           meta$total_entries, " total"),
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "| Family | Success | MSE×100 mean | MSE×100 SD | MAE mean | MAE SD | Upsilon mean | Upsilon SD |",
    "|--------|---------|--------------|-----------|----------|--------|--------------|-----------|"
  )
  for (i in seq_len(nrow(tab))) {
    md <- c(md, paste0(
      "| ", tab$Family[i],
      " | ", tab$Success[i],
      " | ", fmt(tab$MSE100_mean[i], 3),
      " | ", fmt(tab$MSE100_sd[i],   3),
      " | ", fmt(tab$MAE_mean[i],    4),
      " | ", fmt(tab$MAE_sd[i],      4),
      " | ", fmt(tab$Upsilon_mean[i], 4),
      " | ", fmt(tab$Upsilon_sd[i],   4),
      " |"))
  }
  md <- c(md, "",
          "_MSE and MAE are on the held-out TEST fold; Upsilon is on TRAIN (in-sample). Means/SDs are over the successful folds of this batch._",
          "")
  writeLines(md, out_file)
  cat("\nSaved summary table ->", out_file, "\n")

  invisible(tab)
}
