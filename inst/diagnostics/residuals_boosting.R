# ==============================================================================
# Quantile residuals and simulated envelope for PL and APL fits
# Experimental, not part of the master thesis and not exported by the package.
# Source: PLreg_0.4.1/R/PLreg.R  (residuals.PLreg and envelope)
# Requirements: PLreg and PLboost (pAPL, rAPL). A maximum likelihood APL fit is
# the IPLreg() object defined in example_residuals.R, next to this file.
# Load with:
#   source(system.file("diagnostics", "residuals_boosting.R", package = "PLboost"))
# ==============================================================================

# ==============================================================================
# QUANTILE RESIDUAL, PL
# PLreg.R, residuals.PLreg, branch type = "quantile"
# ==============================================================================

res_quantile_PL <- function(y, mu, sigma, lambda, zeta = 2, family) {
  qnorm(pPL(y, mu, sigma, lambda, zeta = zeta, family = family))
}

# ==============================================================================
# RANDOMIZED QUANTILE RESIDUAL, APL
# Queiroz & Ferrari (2024, bimj), eq. for r_i^q. Not in PLreg.R.
# ==============================================================================

res_quantile_APL <- function(y, alpha, mu, sigma, lambda, zeta = 2, family, c) {
  at <- (y == c)
  u  <- numeric(length(y))
  u[!at] <- pAPL(y[!at], alpha[!at], mu[!at], sigma[!at], lambda,
                 zeta = zeta, family = family, c = c)
  u[at]  <- if (c == 1) runif(sum(at), 1 - alpha[at], 1)
            else        runif(sum(at), 0, alpha[at])
  qnorm(u)
}

# ==============================================================================
# PARAMETERS FROM A FITTED MODEL
# PLreg.R, envelope: the muhat/sigmahat/lambdahat block
# ==============================================================================

par_from_fit <- function(object, family, zeta = 2, c = NULL, data = NULL) {
  if (inherits(object, "PLreg")) {
    S <- model.matrix(object, model = "dispersion")
    return(list(mu     = object$fitted.values,
                sigma  = object$link$dispersion$linkinv(
                           as.vector(S %*% object$coefficients$dispersion)),
                lambda = as.numeric(object$coefficients$skewness),
                zeta   = object$zeta, family = object$family))
  }
  if (inherits(object, "IPLreg")) {
    pl <- object$fitPL
    return(list(mu     = pl$link$median$linkinv(as.vector(
                           model.matrix(object$f_mu, data) %*% pl$coefficients$median)),
                sigma  = pl$link$dispersion$linkinv(as.vector(
                           model.matrix(object$f_sigma, data) %*% pl$coefficients$dispersion)),
                lambda = as.numeric(pl$coefficients$skewness),
                alpha  = as.vector(predict(object$fitglm, newdata = data,
                                           type = "response")),
                zeta   = pl$zeta, family = pl$family))
  }
  out <- list(
    mu     = as.vector(fitted(object, parameter = "mu",     type = "response")),
    sigma  = as.vector(fitted(object, parameter = "sigma",  type = "response")),
    lambda = as.vector(fitted(object, parameter = "lambda", type = "response"))[1],
    zeta   = zeta,
    family = family)
  if (!is.null(c))
    out$alpha <- as.vector(fitted(object, parameter = "alpha", type = "response"))
  out
}

# ==============================================================================
# REFIT ON A NEW RESPONSE
# PLreg.R, envelope: the PLreg.fit(X, y_env, S, ...) line
# ==============================================================================

refit_from_fit <- function(object, y_env) {
  if (inherits(object, "IPLreg")) {
    dat <- object$data; dat[[1]] <- y_env
    return(IPLreg(object$formula, data = dat, family = object$fitPL$family,
                  zeta = object$fitPL$zeta, c = object$c))
  }
  cl <- if (inherits(object, "PLreg")) object$call else attr(object, "call")
  dat <- if (inherits(object, "PLreg")) eval(cl$data) else attr(object, "data")
  dat[[1]] <- y_env
  cl$data <- quote(dat)
  if (inherits(object, "PLreg")) cl$type <- object$type
  eval(cl)
}

# ==============================================================================
# SIMULATED ENVELOPE
# PLreg.R, envelope. Only type = "quantile" is kept.
#   object  a PLreg fit, a glmboostLSS fit, or an IPLreg list (APL by ML)
#   data    the data frame the model was fitted on
#   family  density generator, and zeta, for boosting fits only
#   c       NULL for PL, 0 or 1 for APL
# ==============================================================================

envelope_PLboost <- function(object, data, family = NULL, zeta = 2, c = NULL,
                             rep = 40, conf = 0.95, main = "", envcol = "black",
                             plot = TRUE) {
  rep <- max(30, floor(rep))
  if (conf <= 0 | conf >= 1)
    stop("The conf argument must be within the interval (0, 1).", call. = FALSE)

  y   <- data[[1]]
  n   <- length(y)
  par <- par_from_fit(object, family, zeta, c, data)

  resRid <- if (is.null(c))
    res_quantile_PL(y, par$mu, par$sigma, par$lambda, par$zeta, par$family)
  else
    res_quantile_APL(y, par$alpha, par$mu, par$sigma, par$lambda, par$zeta,
                     par$family, c)

  resid_env <- matrix(0, n, rep)
  i   <- 1
  bar <- txtProgressBar(min = 0, max = rep, initial = 0, width = 50,
                        char = "+", style = 3)
  while (i <= rep) {
    tryCatch({
      y_env <- if (is.null(c))
                 rPL(n, par$mu, par$sigma, par$lambda, par$zeta, family = par$family)
               else
                 rAPL(n, par$alpha, par$mu, par$sigma, par$lambda, par$zeta,
                      family = par$family, c = c)
      dat_env <- data; dat_env[[1]] <- y_env
      fit_env <- refit_from_fit(object, y_env)
      p <- par_from_fit(fit_env, family, zeta, c, dat_env)
      res_env <- if (is.null(c))
        res_quantile_PL(y_env, p$mu, p$sigma, p$lambda, p$zeta, p$family)
      else
        res_quantile_APL(y_env, p$alpha, p$mu, p$sigma, p$lambda, p$zeta,
                         p$family, c)
      resid_env[, i] <- sort(res_env)
      setTxtProgressBar(bar, i)
      i <- i + 1
    }, error = function(e) cat("ERROR :", conditionMessage(e), "\n"))
  }
  close(bar); cat("\n")

  liml <- apply(resid_env, 1, quantile, prob = (1 - conf) / 2)
  limu <- apply(resid_env, 1, quantile, prob = (1 - (1 - conf) / 2))
  mm   <- apply(resid_env, 1, median)

  out <- list(res = resRid, liml = liml, limu = limu, mm = mm, main = main,
              prop_outside = mean(sort(resRid) < liml | sort(resRid) > limu))
  if (plot) plot_envelope(out, envcol)
  invisible(out)
}

# ==============================================================================
# ENVELOPE PLOT
# PLreg.R, envelope: the qqnorm block
# ==============================================================================

plot_envelope <- function(e, envcol = "black") {
  ylim <- range(e$res, e$liml, e$limu)
  qqnorm(e$res, pch = "+", las = 1, cex = 0.8, lwd = 3, ylim = ylim,
         xlab = "Quantile N(0,1)", ylab = "Quantile residuals", main = e$main)
  par(new = TRUE)
  qqnorm(e$liml, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim,
         lty = 1, main = "", col = envcol)
  par(new = TRUE)
  qqnorm(e$limu, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim,
         lty = 1, main = "", col = envcol)
  par(new = TRUE)
  qqnorm(e$mm, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim,
         lty = 2, main = e$main)
}
