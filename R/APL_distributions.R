# ==============================================================================
# Adjusted Power Logit (APL) Distribution Functions
# ==============================================================================

#' Adjusted Power Logit Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Adjusted Power Logit (APL) distribution.
#'
#' @param x,q vector of quantiles in [0,1].
#' @param p vector of probabilities in [0,1].
#' @param n number of observations.
#' @param alpha vector of boundary probabilities (point mass at c), in (0,1).
#' @param mu vector of medians of the PL component, in (0,1).
#' @param sigma vector of dispersion parameters of the PL component (> 0).
#' @param lambda scalar skewness parameter of the PL component (>= 0).
#' @param zeta scalar extra parameter (> 0). Default is 2.
#' @param family string specifying the PL family: "NO", "TF", "LO", "PE", "SN", "Hyp", "SLASH".
#' @param c boundary point: either 0 or 1.
#' @param log,log.p logical; if TRUE, probabilities are given as log(p). Default is FALSE.
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x).
#'
#' @details
#' The APL distribution is a mixture:
#' \itemize{
#'   \item With probability \code{alpha}, Y = c (point mass).
#'   \item With probability \code{1 - alpha}, Y follows a PL distribution.
#' }
#'
#' The CDF is:
#' \itemize{
#'   \item If c = 0: F(y) = alpha + (1 - alpha) * pPL(y, ...) for y in (0,1], and F(0) = alpha.
#'   \item If c = 1: F(y) = (1 - alpha) * pPL(y, ...) for y in [0,1), and F(1) = 1.
#' }
#'
#' For \code{dAPL}, observations at y = c return \code{alpha} (the point mass probability).
#'
#' @name APL_distribution
#' @examples
#' # PL with a point mass at zero, normal family
#' dAPL(c(0, 0.3, 0.5, 1),
#'   alpha = 0.3, mu = 0.5, sigma = 1, lambda = 1,
#'   family = "NO", c = 0
#' )
#'
#' # CDF
#' pAPL(0.5, alpha = 0.3, mu = 0.5, sigma = 1, lambda = 1, family = "NO", c = 0)
#'
#' # Quantile
#' qAPL(0.7, alpha = 0.3, mu = 0.5, sigma = 1, lambda = 1, family = "NO", c = 0)
#'
#' # Random generation
#' set.seed(42)
#' y <- rAPL(1000, alpha = 0.3, mu = 0.5, sigma = 1, lambda = 1, family = "NO", c = 0)
#' hist(y[y > 0], prob = TRUE, main = "Continuous component", breaks = 20)
#' mean(y == 0) # should be close to 0.3
# ==============================================================================
# dAPL: density (with point mass at c)
# ==============================================================================
#' @rdname APL_distribution
#' @export
dAPL <- function(x, alpha, mu, sigma, lambda, zeta = 2, family, c,
                 log = FALSE) {
  if (!(c %in% c(0, 1))) {
    stop("c must be either 0 or 1.")
  }
  if (any(alpha <= 0) || any(alpha >= 1)) {
    stop("alpha must be in (0, 1).")
  }

  # Recycle parameters to length of x (lambda and zeta are always scalars)
  n <- max(length(x), length(alpha), length(mu), length(sigma))
  x <- rep_len(x, n)
  alpha <- rep_len(alpha, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)

  dens <- numeric(n)

  at_c <- (x == c)
  in_cont <- (x > 0 & x < 1 & !at_c)
  outside <- (x < 0 | x > 1)

  dens[at_c] <- alpha[at_c]

  dens[in_cont] <- (1 - alpha[in_cont]) *
    dPL(x[in_cont],
      mu     = mu[in_cont],
      sigma  = sigma[in_cont],
      lambda = lambda,
      zeta   = zeta,
      family = family,
      log    = FALSE
    )

  dens[outside] <- 0

  if (log) log(dens) else dens
}

# ==============================================================================
# pAPL: CDF
# ==============================================================================
#' @rdname APL_distribution
#' @export
pAPL <- function(q, alpha, mu, sigma, lambda, zeta = 2, family, c,
                 lower.tail = TRUE, log.p = FALSE) {
  if (!(c %in% c(0, 1))) {
    stop("c must be either 0 or 1.")
  }
  if (any(alpha <= 0) || any(alpha >= 1)) {
    stop("alpha must be in (0, 1).")
  }

  # lambda and zeta are always scalars
  n <- max(length(q), length(alpha), length(mu), length(sigma))
  q <- rep_len(q, n)
  alpha <- rep_len(alpha, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)

  cdf <- numeric(n)

  if (c == 0) {
    # F(y) = alpha                         for y = 0
    # F(y) = alpha + (1-alpha)*pPL(y, ...) for y in (0, 1]
    # F(y) = 0                             for y < 0
    # F(y) = 1                             for y > 1

    below <- q < 0
    at_c <- q == 0
    middle <- q > 0 & q < 1
    above <- q >= 1

    cdf[below] <- 0
    cdf[at_c] <- alpha[at_c]
    cdf[above] <- 1

    cdf[middle] <- alpha[middle] +
      (1 - alpha[middle]) *
        pPL(q[middle],
          mu = mu[middle],
          sigma = sigma[middle],
          lambda = lambda,
          zeta = zeta,
          family = family,
          lower.tail = TRUE, log.p = FALSE
        )
  } else {
    # c == 1
    # F(y) = 0                         for y < 0
    # F(y) = (1-alpha)*pPL(y, ...)     for y in [0, 1)
    # F(y) = 1                         for y >= 1

    below <- q < 0
    middle <- q >= 0 & q < 1
    above <- q >= 1

    cdf[below] <- 0
    cdf[above] <- 1

    cdf[middle] <- (1 - alpha[middle]) *
      pPL(q[middle],
        mu = mu[middle],
        sigma = sigma[middle],
        lambda = lambda,
        zeta = zeta,
        family = family,
        lower.tail = TRUE, log.p = FALSE
      )
  }

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)

  cdf
}

# ==============================================================================
# qAPL: quantile function
# ==============================================================================
#' @rdname APL_distribution
#' @export
qAPL <- function(p, alpha, mu, sigma, lambda, zeta = 2, family, c,
                 lower.tail = TRUE, log.p = FALSE) {
  if (!(c %in% c(0, 1))) {
    stop("c must be either 0 or 1.")
  }
  if (any(alpha <= 0) || any(alpha >= 1)) {
    stop("alpha must be in (0, 1).")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (any(p < 0) || any(p > 1)) {
    stop("p must be between 0 and 1.")
  }

  # lambda and zeta are always scalars
  n <- max(length(p), length(alpha), length(mu), length(sigma))
  p <- rep_len(p, n)
  alpha <- rep_len(alpha, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)

  q <- numeric(n)

  if (c == 0) {
    at_c <- p <= alpha
    middle <- p > alpha & p < 1
    at_1 <- p == 1

    q[at_c] <- 0
    q[at_1] <- 1

    p_cont <- (p[middle] - alpha[middle]) / (1 - alpha[middle])
    q[middle] <- qPL(p_cont,
      mu = mu[middle],
      sigma = sigma[middle],
      lambda = lambda,
      zeta = zeta,
      family = family,
      lower.tail = TRUE
    )
  } else {
    at_1 <- p >= (1 - alpha)
    middle <- p < (1 - alpha) & p > 0
    at_0 <- p == 0

    q[at_1] <- 1
    q[at_0] <- 0

    p_cont <- p[middle] / (1 - alpha[middle])
    q[middle] <- qPL(p_cont,
      mu = mu[middle],
      sigma = sigma[middle],
      lambda = lambda,
      zeta = zeta,
      family = family,
      lower.tail = TRUE
    )
  }

  q
}

# ==============================================================================
# rAPL: random generation
# ==============================================================================
#' @rdname APL_distribution
#' @export
rAPL <- function(n, alpha, mu, sigma, lambda, zeta = 2, family, c) {
  if (!(c %in% c(0, 1))) {
    stop("c must be either 0 or 1.")
  }
  if (any(alpha <= 0) || any(alpha >= 1)) {
    stop("alpha must be in (0, 1).")
  }
  if (n <= 0) {
    stop("n must be a positive integer.")
  }

  # lambda and zeta are always scalars
  alpha <- rep_len(alpha, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)

  at_boundary <- rbinom(n, size = 1, prob = alpha) == 1

  r <- numeric(n)

  r[at_boundary] <- c

  n_cont <- sum(!at_boundary)
  r[!at_boundary] <- rPL(n_cont,
    mu     = mu[!at_boundary],
    sigma  = sigma[!at_boundary],
    lambda = lambda,
    zeta   = zeta,
    family = family
  )

  r
}
