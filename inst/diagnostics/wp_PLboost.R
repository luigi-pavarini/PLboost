# ==============================================================================
# Worm plot for quantile residuals of PL and APL fits
# Experimental, not part of the master thesis and not exported by the package.
# Source: gamlss::wp (van Buuren and Fredriks, 2001), single panel, no xvar.
# Requirements: none. The residuals come from residuals_boosting.R, next to
# this file, or from any other quantile residual.
# Load with:
#   source(system.file("diagnostics", "wp_PLboost.R", package = "PLboost"))
# ==============================================================================

# ==============================================================================
# WORM PLOT
# Detrended normal QQ plot: x is the N(0,1) quantile, y is the ordered residual
# minus that quantile. A good fit stays around zero.
# The band is the one in gamlss::wp. If r ~ N(0,1), the sample quantile of order
# p has standard error sqrt(p (1 - p) / n) / dnorm(z), with z = qnorm(p). It
# depends only on n, so it is the same yardstick for every model.
#   res   quantile residuals
#   ylim  NULL to fit the band in [-3, 3], or a common range for several panels
# Returns, invisibly, the share of points outside the band.
# ==============================================================================

plot_worm <- function(res, main = "", conf = 0.95, ylim = NULL) {
  n  <- length(res)
  qq <- qqnorm(res, plot = FALSE)
  d  <- qq$y - qq$x
  k  <- qnorm(1 - (1 - conf) / 2)
  z  <- seq(-4, 4, by = 0.05)
  se <- sqrt(pnorm(z) * (1 - pnorm(z)) / n) / dnorm(z)
  if (is.null(ylim)) ylim <- range(d, k * se[abs(z) <= 3], -k * se[abs(z) <= 3])

  plot(qq$x, d, pch = "+", las = 1, cex = 0.8, ylim = ylim, xlim = c(-3, 3),
       xlab = "Quantile N(0,1)", ylab = "Deviation", main = main)
  lines(z,  k * se, lty = 2)
  lines(z, -k * se, lty = 2)
  abline(h = 0, col = "grey60")
  legend("bottom", bty = "n",
         legend = c(sprintf("mean %.3f", mean(res)), sprintf("sd %.3f", sd(res))))

  p_x <- pnorm(qq$x)
  invisible(mean(abs(d) > k * sqrt(p_x * (1 - p_x) / n) / dnorm(qq$x)))
}
