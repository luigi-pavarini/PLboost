# ==============================================================================
# cv_load.R — load & concatenate the CV result files for one family.
#
# Adapts automatically to the number of repeats:
#   * R_REPEATS = 1  -> a single batch file, length k   (5 entries)
#   * R_REPEATS > 1  -> one file per batch; each holds a length-(k*R) list with
#                       only its own repeat-range slots filled (rest NULL).
# Files are named <prefix>_cv_<FAM>_<BATCH>.RData and hold `sim_results` + `var_dict`.
#
# Returns: list(sim_results, var_dict, n_entries, n_files)
#   sim_results — combined list (NULLs = failed/unfilled), exactly what the
#                 plotting code expects; downstream filters NULLs as before.
# ==============================================================================
load_cv_results <- function(results_path, prefix, fam) {
  pat   <- paste0("^", prefix, "_cv_", fam, "_.*\\.RData$")
  files <- sort(list.files(results_path, pattern = pat, full.names = TRUE))
  if (length(files) == 0)
    stop("No CV result files for '", fam, "' in ", results_path,
         "  (pattern: ", pat, ")")

  combined <- NULL
  var_dict <- NULL
  for (f in files) {
    e <- new.env(); load(f, envir = e)
    sr <- e$sim_results
    if (is.null(var_dict)) var_dict <- e$var_dict
    if (is.null(combined)) {
      combined <- sr
    } else {
      if (length(sr) > length(combined)) length(combined) <- length(sr)
      fill <- which(!vapply(sr, is.null, logical(1)))
      combined[fill] <- sr[fill]
    }
  }

  n_entries <- sum(!vapply(combined, is.null, logical(1)))
  cat("Loaded", length(files), "file(s) for", fam,
      "->", n_entries, "non-null CV entries (folds x repeats)\n")
  list(sim_results = combined, var_dict = var_dict,
       n_entries = n_entries, n_files = length(files))
}
