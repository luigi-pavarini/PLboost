# ==============================================================================
# Adjusted Power Logit (APL) Families for gamboostLSS/mboost
# ==============================================================================
#
# Author: Luigi Pavarini de Lima (MSc Probability & Statistics, IME-USP)
# Framework: Adjusted Power Logit regression (Queiroz & Ferrari, 2024a)
#
# Implementation of APL distribution family for bounded data in [0,1]
# with potential probability mass at boundaries (0 or 1).
#
# Available family:
# - APL(): Full 4-parameter model (mu, sigma, lambda, alpha)
#
# ==============================================================================

# ==============================================================================
# APL Family for mu (location parameter)
# ==============================================================================
APLMu <- function(sigma = NULL,
                  lambda = NULL,
                  alpha = NULL,
                  family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                  c = 0, # which boundary gets mass
                  zeta = 2,
                  stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  # Loss for mu
  loss <- function(y, f, sigma, lambda, alpha) {
    mu <- plogis(f)
    -dAPL(y,
      alpha = alpha,
      mu = mu,
      sigma = sigma,
      lambda = lambda[1],
      c = c,
      family = family,
      zeta = zeta,
      log = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, sigma = sigma, lambda = lambda, alpha = alpha), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_mu
  ngradient <- function(y, f, w = 1) {
    if (length(unique(lambda)) > 1) {
      stop("lambda must be constant across observations — use lambda ~ 1 in the formula.")
    }
    mu <- plogis(f)

    # Only compute gradients for continuous observations
    continuous <- (y > 0) & (y < 1)
    grad <- rep(0, length(y))

    y_cont <- y[continuous]
    mu_cont <- mu[continuous]
    sigma_cont <- sigma[continuous]
    lambda_cont <- lambda[1]

    v_qts <- v_function(y_cont, mu_cont, sigma_cont, lambda_cont, family, zeta)

    dl_dmu <- (lambda_cont / (sigma_cont * mu_cont * (1 - mu_cont^lambda_cont))) *
      v_qts$z * v_qts$v
    dmu_df <- exp(-f[continuous]) / (1 + exp(-f[continuous]))^2

    grad[continuous] <- dl_dmu * dmu_df

    stabilized <- gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
    stabilized
  }

  response <- function(f) plogis(f)

  offset <- function(y, w = 1) {
    continuous_y <- y[(y > 0) & (y < 1)]
    RET <- qlogis(median(continuous_y))
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("APL: mu (family=%s, c=%d, logit link)", family, c)
  )
}

# ==============================================================================
# APL Family for sigma (dispersion parameter)
# ==============================================================================
APLSigma <- function(mu = NULL,
                     lambda = NULL,
                     alpha = NULL,
                     family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                     c = 0,
                     zeta = 2,
                     stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  # Loss for sigma
  loss <- function(y, f, mu, lambda, alpha) {
    sigma <- exp(f)

    -dAPL(y,
      alpha = alpha,
      mu = mu,
      sigma = sigma,
      lambda = lambda[1],
      c = c,
      family = family,
      zeta = zeta,
      log = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, lambda = lambda, alpha = alpha), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_sigma
  ngradient <- function(y, f, w = 1) {
    if (length(unique(lambda)) > 1) {
      stop("lambda must be constant across observations — use lambda ~ 1 in the formula.")
    }
    sigma <- exp(f)

    # Only compute gradients for continuous observations
    continuous <- (y > 0) & (y < 1)
    grad <- rep(0, length(y))

    y_cont <- y[continuous]
    sigma_cont <- sigma[continuous]
    mu_cont <- mu[continuous]
    lambda_cont <- lambda[1]

    v_qts <- v_function(y_cont, mu_cont, sigma_cont, lambda_cont, family, zeta)

    dl_dsigma <- (1 / sigma_cont) * (v_qts$z^2 * v_qts$v - 1)
    dsigma_df <- sigma_cont

    grad[continuous] <- dl_dsigma * dsigma_df

    stabilized <- gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
    stabilized
  }

  response <- function(f) exp(f)

  offset <- function(y, w = 1) {
    if (is.null(mu)) {
      mu <<- median(median(y[(y > 0) & (y < 1)]))
    } # start for mu
    if (is.null(lambda)) {
      lambda <<- 1
    } # start for lambda
    if (is.null(alpha)) {
      alpha <<- mean(y == c)
    }
    continuous_y <- y[(y > 0) & (y < 1)]
    RET <- optimize(risk,
      interval = log(c(0.1, 10)),
      y = continuous_y, w = w
    )$minimum
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("APL: sigma (family=%s, c=%d, log link)", family, c)
  )
}

# ==============================================================================
# APL Family for lambda (shape parameter)
# ==============================================================================
APLLambda <- function(mu = NULL,
                      sigma = NULL,
                      alpha = NULL,
                      family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                      c = 0,
                      zeta = 2,
                      stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  # Loss for lambda
  loss <- function(y, f, mu, sigma, alpha) {
    lam <- exp(f)

    -dAPL(y,
      alpha = alpha,
      mu = mu,
      sigma = sigma,
      lambda = lam[1],
      c = c,
      family = family,
      zeta = zeta,
      log = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, sigma = sigma, alpha = alpha), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_lambda
  ngradient <- function(y, f, w = 1) {
    lam <- exp(f)
    if (length(unique(lam)) > 1) {
      stop("lambda must be constant across observations — use lambda ~ 1 in the formula.")
    }

    # Only compute gradients for continuous observations
    continuous <- (y > 0) & (y < 1)
    grad <- rep(0, length(y))

    y_cont <- y[continuous]
    lam_cont <- lam[1]
    mu_cont <- mu[continuous]
    sigma_cont <- sigma[continuous]

    v_qts <- v_function(y_cont, mu_cont, sigma_cont, lam_cont, family, zeta)

    dl_dlam <- 1 / lam_cont +
      (y_cont^lam_cont * log(y_cont)) / (1 - y_cont^lam_cont) -
      (1 / sigma_cont) * v_qts$z * v_qts$v *
        (log(y_cont) / (1 - y_cont^lam_cont) - log(mu_cont) / (1 - mu_cont^lam_cont))

    dlam_df <- lam_cont

    grad[continuous] <- dl_dlam * dlam_df

    gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
  }

  response <- function(f) exp(f)

  offset <- function(y, w = 1) {
    if (is.null(mu)) {
      mu <<- median(y[(y > 0) & (y < 1)])
    } # start for mu
    if (is.null(sigma)) {
      sigma <<- 10
    } # start for sigma
    if (is.null(alpha)) {
      alpha <<- mean(y == c)
    }
    continuous_y <- y[(y > 0) & (y < 1)]
    RET <- optimize(risk,
      interval = log(c(0.1, 2)),
      y = continuous_y, w = w
    )$minimum
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("APL: lambda (family=%s, c=%d, log link)", family, c)
  )
}

# ==============================================================================
# APL Family for alpha (boundary probability)
# ==============================================================================
APLAlpha <- function(mu = NULL,
                     sigma = NULL,
                     lambda = NULL,
                     family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                     c = 0,
                     zeta = 2,
                     stabilization = c("none", "MAD", "L2")) {
  family <- match.arg(family)
  stabilization <- match.arg(stabilization)

  if (!(c %in% c(0, 1))) {
    stop("c must be 0 (mass at 0) or 1 (mass at 1)")
  }

  # Loss for alpha
  loss <- function(y, f, mu, sigma, lambda) {
    alpha <- plogis(f)

    -dAPL(y,
      alpha = alpha,
      mu = mu,
      sigma = sigma,
      lambda = lambda[1],
      c = c,
      family = family,
      zeta = zeta,
      log = TRUE
    )
  }

  # Empirical risk
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, mu = mu, sigma = sigma, lambda = lambda), na.rm = TRUE)
  }

  # Negative gradient wrt f = eta_alpha
  ngradient <- function(y, f, w = 1) {
    alpha <- plogis(f)

    # Identify boundary and continuous observations
    if (c == 0) {
      at_boundary <- (y == 0)
      ind_c <- ifelse(y == 0, 1, 0)
    } else {
      at_boundary <- (y == 1)
      ind_c <- ifelse(y == 1, 1, 0)
    }

    dl_dalpha <- (ind_c - alpha) / (alpha * (1 - alpha))
    dalpha_df <- exp(-f) / (1 + exp(-f))^2

    grad <- dl_dalpha * dalpha_df

    stabilized <- gamboostLSS:::stabilize_ngradient(grad, w = w, stabilization = stabilization)
    stabilized
  }

  response <- function(f) plogis(f)

  offset <- function(y, w = 1) {
    # Calculate proportion at boundary
    if (c == 0) {
      prop_boundary <- mean(y == 0, na.rm = TRUE)
    } else {
      prop_boundary <- mean(y == 1, na.rm = TRUE)
    }

    # Start with observed proportion, bounded away from 0 and 1
    prop_boundary <- pmax(pmin(prop_boundary, 0.99), 0.01)

    RET <- qlogis(prop_boundary)
    return(RET)
  }

  mboost::Family(
    ngradient = ngradient,
    risk = risk,
    loss = loss,
    response = response,
    offset = offset,
    name = sprintf("APL: alpha (family=%s, c=%d, logit link)", family, c)
  )
}

# ==============================================================================
# APL Family Constructor (Full 4-parameter model)
# ==============================================================================
# 4-parameter (mu, sigma, lambda, alpha) Adjusted Power Logit
# Example:
#   fam <- APL(family = "NO", c = 0)  # Normal errors, mass at 0
#   fit <- glmboostLSS(y ~ X | X | X | X, families = fam, data = df)
# ==============================================================================
APL <- function(stabilization = c("none", "MAD", "L2"),
                family = c("NO", "LO", "TF", "PE", "Hyp", "SN", "SLASH"),
                c = 0, # which boundary gets mass: 0 or 1
                zeta = 2) {
  stabilization <- match.arg(stabilization)
  family <- match.arg(family)

  if (!(c %in% c(0, 1))) {
    stop("c must be 0 (mass at 0) or 1 (mass at 1)")
  }

  Families(
    mu = APLMu(
      sigma = NULL,
      lambda = NULL,
      alpha = NULL,
      family = family,
      c = c,
      zeta = zeta,
      stabilization = stabilization
    ),
    sigma = APLSigma(
      mu = NULL,
      lambda = NULL,
      alpha = NULL,
      family = family,
      c = c,
      zeta = zeta,
      stabilization = stabilization
    ),
    lambda = APLLambda(
      mu = NULL,
      sigma = NULL,
      alpha = NULL,
      family = family,
      c = c,
      zeta = zeta,
      stabilization = stabilization
    ),
    alpha = APLAlpha(
      mu = NULL,
      sigma = NULL,
      lambda = NULL,
      family = family,
      c = c,
      zeta = zeta,
      stabilization = stabilization
    )
  )
}

# ==============================================================================
# USAGE GUIDE: Adjusted Power Logit (APL) Models
# ==============================================================================
#
# Available Model:
# - APL():  Full 4-parameter model (mu, sigma, lambda, alpha)
#
# Parameter c controls which boundary gets probability mass:
# - c = 0: Mass at 0 (point mass at zero)
# - c = 1: Mass at 1 (point mass at one)
#
# Available error distributions:
# "NO" = Normal, "LO" = Logistic, "TF" = Student-t, "PE" = Power Exponential,
# "SN" = Sinh-Normal, "Hyp" = Hyperbolic, "SLASH" = Slash
# ==============================================================================

# ==============================================================================
# BASIC EXAMPLES
# ==============================================================================

# # Data with exact 0s
# fam_apl_0 <- APL(family = "NO", c = 0)  # 4-parameter, mass at 0
# fit_apl_0 <- glmboostLSS(
#   formula = list(mu = y ~ ., sigma = y ~ ., lambda = y ~ ., alpha = y ~ .),
#   data = data,
#   families = fam_apl_0,
#   control = boost_control(mstop = 150, nu = 0.01)
# )
#
# # Data with exact 1s
# fam_apl_1 <- APL(family = "TF", c = 1, zeta = 3)  # 4-parameter, mass at 1, robust
# fit_apl_1 <- glmboostLSS(
#   formula = list(mu = y ~ ., sigma = y ~ ., lambda = y ~ ., alpha = y ~ .),
#   data = data,
#   families = fam_apl_1,
#   control = boost_control(mstop = 150, nu = 0.01)
# )

# ==============================================================================
# PARAMETER INTERPRETATION
# ==============================================================================

# μ (mu):     Location parameter for continuous part
# σ (sigma):  Dispersion parameter for continuous part
# λ (lambda): Shape/asymmetry parameter for continuous part
# α (alpha):  Probability of boundary values (0 or 1)

# Transform predictions:
# plogis(pred$mu)     # mu ∈ (0,1)
# exp(pred$sigma)     # sigma > 0
# exp(pred$lambda)    # lambda > 0
# plogis(pred$alpha)  # alpha ∈ (0,1)
