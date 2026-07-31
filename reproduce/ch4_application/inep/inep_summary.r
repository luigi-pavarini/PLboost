# ==============================================================================
# inep_summary_combined.R
# Combine ALL 10rep batch files into ONE summary table over the full 50 fold-
# entries (5 folds x 10 repeats). INEP ES is split into THREE batches:
#   _1 = repeats 1-3, _2 = repeats 4-6, _3 = repeats 7-10  (slots 1-15/16-30/31-50)
#
# load_cv_results() globs EVERY batch tag (_1, _2, _3, ...) and merges by slot,
# so this script is batch-count-agnostic. Pure post-processing — reads results
# only, runs no models. Mirrors credit_summary_combined.r.
#
# Output: summary_inep_ES_cv_ALL.md  (next to the .RData files)
# ==============================================================================

# ── PER-APP PARAMETERS (only this block differs between apps) ─────────────────
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
app_name     <- "inep_appl"
prefix       <- "sim_results_medicina_SP"            # load_cv_results prefix (matches inep_run_cv_10rep.r output)
family_order <- c("NO", "TF", "LO", "PE", "SN", "SLASH", "Hyp", "Beta")
K_FOLDS      <- 5
R_REPEATS    <- 10
# ─────────────────────────────────────────────────────────────────────────────

results_path <- file.path(app_path, "cross_val", "results")
source(file.path(app_path, "essentials", "cv_load.R"))
source(file.path(app_path, "essentials", "cv_summary.R"))

TOTAL_ENTRIES <- K_FOLDS * R_REPEATS                         # 50 = error-bar population

# Load each family's MERGED results (all batch tags _1/_2/_3). If a family has no files
# yet, keep an all-NULL slot vector so the table still renders ( "—" / 0/50 ).
results_containers <- setNames(lapply(family_order, function(fam) {
  tryCatch(
    load_cv_results(results_path, prefix, fam)$sim_results,
    error = function(e) {
      message("  [skip] ", fam, ": ", conditionMessage(e))
      vector("list", TOTAL_ENTRIES)
    }
  )
}), family_order)

# Single summary over ALL entries: batch_entries = total_entries so the Success
# column and footer read out of 50, and SD is over the pooled folds.
summarise_cv(
  results_containers = results_containers,
  family_order       = family_order,
  meta = list(app = app_name, k_folds = K_FOLDS, r_repeats = R_REPEATS,
              repeat_start = 1, repeat_end = R_REPEATS,
              batch_entries = TOTAL_ENTRIES, total_entries = TOTAL_ENTRIES,
              batch_tag = "ALL"),
  out_file = file.path(results_path,
                       paste0("summary_", sub("^sim_results_", "", prefix),
                              "_cv_ALL.md"))
)
