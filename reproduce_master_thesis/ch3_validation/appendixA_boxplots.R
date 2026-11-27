# ==============================================================================
# appendixA_boxplots.R
#
# Figures A.1 to A.5 of the master thesis: boosting and maximum likelihood
# estimates across the 100 replications of the validation in Sections 3.1 and
# 3.2, one figure per scenario and one panel per submodel.
#
# Reads, from $PLBOOST_OUT (default: the working directory), the two files
# written by validation_multirep.R:
#   validation_replications.csv   one row per replication and coefficient
#   validation_results.csv        true values
# Writes $PLBOOST_OUT/appA_validation/box_<scenario>_<submodel>.png
#
# Run:  Rscript appendixA_boxplots.R
# ==============================================================================

base <- Sys.getenv("PLBOOST_OUT", unset = ".")
out  <- file.path(base, "appA_validation")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
reps <- read.csv(file.path(base, "validation_replications.csv"), stringsAsFactors = FALSE)
res  <- read.csv(file.path(base, "validation_results.csv"),      stringsAsFactors = FALSE)

# one entry per submodel: coefficients, axis labels and y-axis symbol
B_mu <- list(id = "mu", sym = quote(eta[mu]),
  ks  = c("mu.(Intercept)", "mu.V1", "mu.V2", "mu.V3", "mu.V4"),
  rot = list(expression("Intercept"), expression(x[1]), expression(x[2]),
             expression(x[3]), expression(x[4])))
B_sg <- list(id = "sigma", sym = quote(eta[sigma]),
  ks  = c("sig.(Intercept)", "sig.V3", "sig.V4", "sig.V5"),
  rot = list(expression("Intercept"), expression(x[3]), expression(x[4]), expression(x[5])))
B_lb <- list(id = "lambda", sym = quote(log ~ lambda), ks = "lambda", rot = list(expression(lambda)))
B_al <- list(id = "alpha", sym = quote(eta[alpha]),
  ks  = c("al.(Intercept)", "al.U1", "al.U2"),
  rot = list(expression("Intercept"), expression(t[1]), expression(t[2])))

scenarios <- list(
  list(name = "PL Scenario 1",            tag = "pl1",      bls = list(B_mu, B_sg, B_lb)),
  list(name = "PL Scenario 2",            tag = "pl2",      bls = list(B_mu, B_sg, B_lb)),
  list(name = "APL High boundary mass",   tag = "apl_high", bls = list(B_mu, B_sg, B_lb, B_al)),
  list(name = "APL Medium boundary mass", tag = "apl_med",  bls = list(B_mu, B_sg, B_lb, B_al)),
  list(name = "APL Low boundary mass",    tag = "apl_low",  bls = list(B_mu, B_sg, B_lb, B_al)))

# each coefficient gets a cell with the boosting and MLE boxes side by side and
# its own dashed line at the true value; the y range is common to the panel.
# The log(lambda) panel has a single cell, drawn square and centred.
panel <- function(scen, bl, file) {
  d  <- reps[reps$scenario == scen, ]
  tv <- res[res$scenario == scen, c("coefficient", "true")]
  n  <- length(bl$ks); ncell <- max(n, 2)
  yl <- extendrange(d$estimate[d$coefficient %in% bl$ks], f = if (bl$id == "lambda") 0.25 else 0.03)
  boxes <- list(); at <- c(); cols <- c()
  for (j in seq_len(n)) {
    dk <- d[d$coefficient == bl$ks[j], ]
    boxes[[length(boxes) + 1]] <- dk$estimate[dk$estimator == "boost"]
    boxes[[length(boxes) + 1]] <- dk$estimate[dk$estimator == "mle"]
    at <- c(at, 2.4 * j - 1.5, 2.4 * j - 0.65); cols <- c(cols, "#cfdced", "#eadcc8")
  }
  png(file, width = if (bl$id == "lambda") 1200 else round(300 * ncell + 300), height = 1200, res = 170)
  op <- par(mar = c(2.9, 4.2, 1.4, 1), mgp = c(2.7, 0.7, 0))
  boxplot(boxes, at = at, xaxt = "n", las = 1, col = cols, border = "grey25", outcex = 0.4,
          boxwex = 0.80, xlim = if (bl$id == "lambda") 1.325 + c(-2.4, 2.4) else c(0.35, 2.4 * ncell + 0.35),
          ylim = yl, ylab = bl$sym, cex.axis = 0.85, cex.lab = 1.05)
  for (j in seq_len(n)) {
    t0 <- tv$true[tv$coefficient == bl$ks[j]]
    if (bl$id == "lambda") segments(1.325 - 1.2, t0, 1.325 + 1.2, t0, lty = 2, col = "red", lwd = 1.3)
    else segments(2.4 * j - 2.05, t0, 2.4 * j + 0.35, t0, lty = 2, col = "red", lwd = 1.3)
    mtext(bl$rot[[j]], side = 1, line = 0.5, at = 2.4 * j - 1.075, cex = 0.85)
  }
  if (n > 1) for (j in 1:(n - 1)) abline(v = 2.4 * j + 0.35, col = "grey82", lwd = 0.8)
  legend("topright", c("Boosting", "MLE", "true"), pch = c(22, 22, NA),
         pt.bg = c("#cfdced", "#eadcc8", NA), pt.cex = 1.4, lty = c(NA, NA, 2),
         col = c("grey25", "grey25", "red"), lwd = c(NA, NA, 1.4), bty = "n", cex = 0.75,
         horiz = TRUE, x.intersp = 0.7, seg.len = 1.5)
  par(op); invisible(dev.off())
}

for (sc in scenarios) for (bl in sc$bls)
  panel(sc$name, bl, file.path(out, sprintf("box_%s_%s.png", sc$tag, bl$id)))
cat("wrote", length(list.files(out, pattern = "^box_")), "figures to", out, "\n")
