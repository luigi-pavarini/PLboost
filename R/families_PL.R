# ==============================================================================
# Power-Logit families (via PLreg) for use with gamboostLSS/mboost
# Links:
#   - mu     : logit (μ ∈ (0,1))
#   - sigma  : log   (σ > 0)
#   - lambda : log   (λ > 0)  ➜ kept FIXED (no regression on λ)
#
# Notes:
#   * We include λ as a model parameter (log link), but DO NOT fit a
#     regression for it; treat λ as a hyperparameter set once (e.g., λ = 1.2).
#   * Density via PLreg:::dPL (advisor's implementation).
# Requirements: PLreg, mboost, gamboostLSS
# ==============================================================================

#' Power logit families of each parameter
#'
#' The \code{mboost} \code{Family} objects behind \code{\link{PowerLogit}}
#' and \code{\link{GJS}}, one per parameter. Most users call the
#' constructors. These are exported for custom combinations.
#'
#' @param mu,sigma,lambda value of the other parameters, or \code{NULL} to
#'   take them from the current fit.
#' @param family density generator of the power logit class: \code{"NO"}
#'   (normal), \code{"LO"} (logistic), \code{"TF"} (Student t), \code{"PE"}
#'   (power exponential), \code{"Hyp"} (hyperbolic), \code{"SN"}
#'   (sinh-normal) or \code{"SLASH"} (slash).
#' @param zeta extra shape parameter of the generator, used by \code{"TF"},
#'   \code{"PE"}, \code{"Hyp"}, \code{"SN"} and \code{"SLASH"}.
#' @param stabilization scaling of the negative gradient: \code{"none"},
#'   \code{"MAD"} or \code{"L2"}, as in \code{gamboostLSS}.
#' @return A \code{Family} object.
#' @name PL_submodels
NULL

# ==============================================================================
# Custom Family for mu (location)
# eta_mu = f; mu = plogis(f)
# Chain rule: dℓ/df = (dℓ/dμ) * (dμ/df), with dμ/df = logistic'(f)
# ==============================================================================
#' @rdname PL_submodels
#' @export
PLMu <- function(sigma = NULL,
                 lambda = NULL,
                 family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                 zeta = 2,
                 stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  # Loss for mu
  loss <- function(y, f, sigma, lambda) {
    mu <- plogis(f)

    -PLreg:::dPL(y,
      mu     = mu,
      sigma  = sigma,
      lambda = lambda[1],
      family = family,
      zeta   = zeta,
      log    = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    # print(loss(y = y, f = f, sigma = sigma, lambda = lambda))
    sum(w * loss(y = y, f = f, sigma = sigma, lambda = lambda), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_mu
  ngradient <- function(y, f, w = 1) {
    if (length(unique(lambda)) > 1) {
      stop("lambda must be constant across observations — use lambda ~ 1 in the formula.")
    }
    mu <- plogis(f)
    v_qts <- v_function(y, mu, sigma, lambda, family, zeta)
    dl_dmu <- (lambda / (sigma * mu * (1 - mu^lambda))) * v_qts$z * v_qts$v
    dmu_df <- exp(-f) / (1 + exp(-f))^2
    grad <- dl_dmu * dmu_df
    stabilized <- gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
    # print(stabilized)
    stabilized
  }

  response <- function(f) plogis(f)

  offset <- function(y, w = 1) {
    RET <<- qlogis(median(y))
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("Power-Logit: mu (family=%s, logit link)", family)
  )
}

# ==============================================================================
# Custom Family for sigma (dispersion)
# eta_sigma = f; sigma = exp(f)
# Chain rule: dℓ/df = (dℓ/dσ) * (dσ/df), with dσ/df = exp(f) = σ
# ==============================================================================
#' @rdname PL_submodels
#' @export
PLSigma <- function(mu = NULL,
                    lambda = NULL,
                    family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                    zeta = 2,
                    stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  # Loss for sigma
  loss <- function(y, f, mu, lambda) {
    sigma <- exp(f)

    -PLreg:::dPL(y,
      mu     = mu,
      sigma  = sigma,
      lambda = lambda[1],
      family = family,
      zeta   = zeta,
      log    = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    # print(loss(y = y, f = f, mu = mu, lambda = lambda))
    sum(w * loss(y = y, f = f, mu = mu, lambda = lambda), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_sigma
  ngradient <- function(y, f, w = 1) {
    if (length(unique(lambda)) > 1) {
      stop("lambda must be constant across observations — use lambda ~ 1 in the formula.")
    }
    sigma <- exp(f)
    v_qts <- v_function(y, mu, sigma, lambda, family, zeta)
    dl_dsigma <- (1 / sigma) * (v_qts$z^2 * v_qts$v - 1)
    dsigma_df <- sigma
    grad <- dl_dsigma * dsigma_df
    stabilized <- gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
    stabilized
  }

  response <- function(f) exp(f)

  offset <- function(y, w = 1) {
    if (is.null(mu)) {
      mu <<- median(y)
    } # start for mu
    if (is.null(lambda)) {
      lambda <<- 1
    } # start for lambda
    RET <<- optimize(risk,
      interval = log(c(0.1, 10)), # sigma in [0.1, 10]
      y = y, w = w
    )$minimum
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("Power-Logit: sigma (family=%s, log link)", family)
  )
}

# ==============================================================================
# Custom Family for lambda (shape / asymmetry)
# eta_lambda = g; lambda = exp(g)
# Chain rule: dℓ/df = (dℓ/dλ) * (dλ/df), with dλ/df = exp(g) = λ
# ==============================================================================
#' @rdname PL_submodels
#' @export
PLLambda <- function(mu = NULL,
                     sigma = NULL,
                     family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                     zeta = 2,
                     stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  # Loss for lambda
  loss <- function(y, f, mu, sigma) {
    lam <- exp(f)

    -PLreg:::dPL(y,
      mu     = mu,
      sigma  = sigma,
      lambda = lam[1],
      family = family,
      zeta   = zeta,
      log    = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, sigma = sigma), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_lambda
  ngradient <- function(y, f, w = 1) {
    lam <- exp(f)
    if (length(unique(lam)) > 1) {
      stop("lambda must be constant across observations — use lambda ~ 1 in the formula.")
    }
    v_qts <- v_function(y, mu, sigma, lam, family, zeta)
    dl_dlam <- 1 / lam +
      (y^lam * log(y)) / (1 - y^lam) -
      (1 / sigma) * v_qts$z * v_qts$v *
        (log(y) / (1 - y^lam) - log(mu) / (1 - mu^lam))
    dlam_df <- exp(f)
    grad <- dl_dlam * dlam_df
    gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
  }

  response <- function(f) exp(f)

  offset <- function(y, w = 1) {
    if (is.null(mu)) {
      mu <<- median(y)
    } # start for mu
    if (is.null(sigma)) {
      sigma <<- 10
    } # start for sigma
    RET <<- optimize(risk,
      interval = log(c(0.1, 2)), # lambda in [0.1, 2]
      y = y, w = w
    )$minimum
    # RET <<- log(0.1)   # removed: use the optimized offset (like APL), do not override
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("Power-Logit: lambda (family=%s, log link)", family)
  )
}

#' Power logit families for gamboostLSS
#'
#' Build the \code{Families} object of the power logit (PL) distribution for
#' \code{\link[gamboostLSS]{glmboostLSS}} and
#' \code{\link[gamboostLSS]{gamboostLSS}}. \code{PowerLogit()} has three
#' submodels: the median \code{mu} (logit link), the dispersion \code{sigma}
#' (log link) and the skewness \code{lambda} (log link). \code{GJS()} is the
#' special case \code{lambda = 1} (generalized Johnson SB), with submodels
#' \code{mu} and \code{sigma}.
#'
#' @param family density generator of the power logit class: \code{"NO"}
#'   (normal), \code{"LO"} (logistic), \code{"TF"} (Student t), \code{"PE"}
#'   (power exponential), \code{"Hyp"} (hyperbolic), \code{"SN"}
#'   (sinh-normal) or \code{"SLASH"} (slash).
#' @param zeta extra shape parameter of the generator, used by \code{"TF"},
#'   \code{"PE"}, \code{"Hyp"}, \code{"SN"} and \code{"SLASH"}.
#' @param stabilization scaling of the negative gradient: \code{"none"},
#'   \code{"MAD"} or \code{"L2"}, as in \code{gamboostLSS}.
#'
#' @details \code{lambda} is a scalar shared by all observations, so it must
#'   enter the formula as \code{lambda = y ~ 1}. The families stop otherwise.
#'
#' @return A \code{Families} object.
#' @seealso \code{\link{APL}} for a response with a point mass at 0 or 1,
#'   \code{\link{PL_submodels}} for the families of each parameter.
#' @examples
#' library(gamboostLSS)
#' set.seed(1206)
#' n <- 200
#' x <- rnorm(n)
#' y <- PLreg::rPL(n, mu = plogis(0.5 + x), sigma = 1, lambda = 1.2,
#'                 family = "NO")
#' fit <- glmboostLSS(
#'   list(mu = y ~ x, sigma = y ~ x, lambda = y ~ 1),
#'   data = data.frame(y, x),
#'   families = PowerLogit(family = "NO", stabilization = "MAD"),
#'   control = boost_control(mstop = 300, nu = 0.01)
#' )
#' coef(fit)
#' @export
PowerLogit <- function(stabilization = c("none", "MAD", "L2"),
                       family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                       zeta = 2) {
  stabilization <- match.arg(stabilization)
  family <- match.arg(family)
  Families(
    mu = PLMu(
      sigma = NULL,
      lambda = NULL,
      family = family,
      zeta = zeta,
      stabilization = stabilization
    ),
    sigma = PLSigma(
      mu = NULL,
      lambda = NULL,
      family = family,
      zeta = zeta,
      stabilization = stabilization
    ),
    lambda = PLLambda(
      mu = NULL,
      sigma = NULL,
      family = family,
      zeta = zeta,
      stabilization = stabilization
    )
  )
}

#' @rdname PowerLogit
#' @export
GJS <- function(stabilization = c("none", "MAD", "L2"),
                family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                zeta = 2) {
  stabilization <- match.arg(stabilization)
  family <- match.arg(family)

  Families(
    mu = PLMu(
      sigma = NULL,
      lambda = 1,
      family = family,
      zeta = zeta,
      stabilization = stabilization
    ),
    sigma = PLSigma(
      mu = NULL,
      lambda = 1,
      family = family,
      zeta = zeta,
      stabilization = stabilization
    )
  )
}
