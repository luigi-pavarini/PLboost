# ==============================================================================
# selection_table.R
# Selection table of an application: for every covariate and submodel, the
# proportion of the 50 CV fits that selected it (Freq.), the average
# coefficient over the 50 fits (zeros included) and the variability band
# avg +/- 1.96 sd across fits.
#
# The dissertation tables (Chapter 4) show a hand-picked subset of these rows,
# ordered by frequency; every printed value is reproduced by this script.
#
# Usage (from the repository root, PLBOOST_OUT as for the figure scripts):
#   PLBOOST_OUT=... Rscript reproduce/ch4_application/essentials/selection_table.R \
#       reproduce/ch4_application/credit/credit_img_ipl.r
#   PLBOOST_OUT=... Rscript reproduce/ch4_application/essentials/selection_table.R \
#       reproduce/ch4_application/inep/inep_img.r
#
# The figure script is sourced first because it builds the coefficient matrices
# (coef_mu_mat, coef_sigma_mat and, for the inflated model, coef_alpha_mat);
# its figures are regenerated as a side effect.
# Output: $PLBOOST_OUT/cross_val/selection_table.csv and selection_table.tex
# ==============================================================================

img_script <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(img_script)) stop("pass the figure script (credit_img_ipl.r or inep_img.r)")
source(img_script)

mats <- list(mu = coef_mu_mat, sigma = coef_sigma_mat)
if (exists("coef_alpha_mat")) mats$alpha <- coef_alpha_mat

tab <- do.call(rbind, lapply(names(mats), function(k) {
  m <- mats[[k]]
  m <- m[!grepl("Intercept", rownames(m)), , drop = FALSE]
  avg <- rowMeans(m); s <- apply(m, 1, sd)
  data.frame(submodel = k, variable = map_names(rownames(m), var_dict),
             freq = rowMeans(m != 0), avg = avg,
             lower = avg - 1.96 * s, upper = avg + 1.96 * s)
}))
tab <- tab[tab$freq > 0, ]
tab <- tab[order(match(tab$submodel, names(mats)), -tab$freq, -abs(tab$avg)), ]

out_dir <- file.path(Sys.getenv("PLBOOST_OUT", unset = "."), "cross_val")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(tab, file.path(out_dir, "selection_table.csv"), row.names = FALSE)

num <- function(x) ifelse(x < 0, sprintf("$-$%.3f", abs(x)), sprintf("%.3f", x))
tex <- sprintf("  %s & \\texttt{%s} & %.2f & %s & %s & %s \\\\",
               ifelse(duplicated(tab$submodel), "",
                      sprintf("$\\eta_\\%s$", tab$submodel)),
               gsub("_", "\\\\_", tab$variable), tab$freq,
               num(tab$avg), num(tab$lower), num(tab$upper))
writeLines(tex, file.path(out_dir, "selection_table.tex"))
cat("Wrote", nrow(tab), "rows to", out_dir, "\n")
