# ==============================================================================
# inep_table_pl_cv.r
# Renders the pooled CV comparison table (all batches) as a PNG image, for
# pasting into the email / dissertation. Same numbers as
# summary_medicina_SP_cv_ALL.md, but printed as a figure. Winner row (PE)
# highlighted. Saves to imgs/PE/plot_table_comparison.png.
# ==============================================================================

suppressMessages({library(gridExtra); library(grid)})
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
results_path <- file.path(app_path, "cross_val", "results")
out_png      <- file.path(app_path, "cross_val", "imgs", "NO", "plot_table_comparison.png")

source(file.path(app_path, "essentials", "cv_load.R"))

prefix   <- "sim_results_medicina_SP"
fams     <- c("NO", "TF", "LO", "PE", "SN", "SLASH", "Hyp", "Beta")
WINNER   <- "NO"

err_of <- function(x) (if (!is.null(x$error_pl)) x$error_pl else x$error_beta) * 100
mae_of <- function(x)  if (!is.null(x$mae_pl))   x$mae_pl   else x$mae_beta

rows <- lapply(fams, function(f) {
  sr <- tryCatch(load_cv_results(results_path, prefix, f)$sim_results, error = function(e) list())
  ok <- sr[!vapply(sr, is.null, logical(1))]
  if (length(ok) == 0) return(data.frame(Familia = f, Sucesso = "0/50",
                                          MSEx100 = "—", MAE = "—", Upsilon = "—"))
  mse <- vapply(ok, err_of, numeric(1)); mae <- vapply(ok, mae_of, numeric(1))
  ups <- vapply(ok, function(x) x$upsilon, numeric(1))
  ms  <- function(v, d) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(v), sd(v))
  data.frame(Familia = if (f == "Beta") "Beta" else paste0("PL (", f, ")"),
             Sucesso = paste0(length(ok), "/50"),
             MSEx100 = ms(mse, 3), MAE = ms(mae, 4), Upsilon = ms(ups, 4),
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
# Build accented / special glyphs from code points so the source stays pure ASCII
# and survives a non-UTF-8 Rscript locale.  0xED=i-acute  0xD7=times  0x3A5=Upsilon  0xE9=e-acute
u <- function(cp) intToUtf8(cp)
colnames(tab) <- c(paste0("Fam", u(0xED), "lia"), "Sucesso",
                   paste0("MSE", u(0xD7), "100"), "MAE", "Upsilon")
win_row <- which(fams == WINNER)

# ── Theme: clean, winner row highlighted ─────────────────────────────────────
th <- ttheme_minimal(
  core = list(
    fg_params = list(cex = 0.9, fontface = ifelse(seq_len(nrow(tab)) == win_row, 2, 1)),
    bg_params = list(fill = ifelse(seq_len(nrow(tab)) == win_row, "#F4D9D9", "white"),
                     col = "gray85", lwd = 1)),
  colhead = list(fg_params = list(cex = 0.95, fontface = 2, col = "white"),
                 bg_params = list(fill = "#4A4A4A")))

g   <- tableGrob(tab, rows = NULL, theme = th)
cap_txt <- paste0("M", u(0xE9), "dia (DP) sobre os 50 folds; MSE e MAE no teste, ",
                  "Upsilon no treino. Vencedor: PL (NO).")
cap <- textGrob(cap_txt, gp = gpar(cex = 0.62, col = "gray40"))
full <- arrangeGrob(g, cap, ncol = 1,
                    heights = unit.c(sum(g$heights), unit(2.2, "lines")))

# Size the canvas to the grob so nothing is clipped (measure in a temp device).
tmp <- tempfile(fileext = ".png"); png(tmp)
gw <- convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
gh <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
dev.off(); unlink(tmp)

png(out_png, width = (gw + 0.4) * 300, height = (gh + 0.9) * 300, res = 300,
    type = "quartz")   # quartz renders UTF-8 (accents, ×, Greek Υ) on macOS
grid.newpage(); grid.draw(full)
dev.off()
cat("Saved:", out_png, "  (", round(gw + 0.4, 1), "x", round(gh + 0.9, 1), "in )\n")
