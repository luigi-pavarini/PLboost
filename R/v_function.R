#' Score helpers of the power logit density
#'
#' Standardized latent residual \code{z}, weight function \code{v(z)} of the
#' density generator and its derivative, used by the gradients of the PL and
#' APL families.
#'
#' @param y response values in (0, 1).
#' @param mu median, in (0, 1).
#' @param sigma dispersion, a scalar or one value per observation.
#' @param lambda skewness, a scalar. With \code{lambda = 0} the log-log limit
#'   of the power logit is used.
#' @param family density generator: \code{"NO"}, \code{"TF"}, \code{"LO"},
#'   \code{"SN"}, \code{"PE"}, \code{"Hyp"} or \code{"SLASH"}.
#' @param zeta extra shape parameter, used by \code{"TF"}, \code{"SN"},
#'   \code{"PE"}, \code{"Hyp"} and \code{"SLASH"}.
#'
#' @return A list with \code{z}, \code{v} and \code{dv}, each of the length
#'   of \code{y}.
#' @examples
#' v_function(c(0.2, 0.5, 0.8), mu = 0.5, sigma = 1, lambda = 1.2,
#'            family = "NO")
#' @export
v_function <- function(y, mu, sigma, lambda, family, zeta = 2) {
  # === Input validation =======================================================
  valid_families <- c("NO", "TF", "LO", "SN", "PE", "Hyp", "SLASH")
  if (!(family %in% valid_families)) {
    stop("Error: 'family' must be one of: NO, TF, LO, SN, PE, Hyp, SLASH.")
  }

  # === Compute z ==============================================================
  if (lambda[1] == 0) {
    z <- (-1 / sigma) * (log(log(y) / log(mu)))
  } else {
    z <- (1 / sigma) * (VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
  }

  # === Compute v(z) and v'(z) by family =======================================
  if (family == "NO") {
    vz <- rep(1, length(y))
    dvz <- rep(0, length(y))
  }

  if (family == "TF") {
    vz <- (zeta + 1) / (zeta + z^2)
    dvz <- -2 * (zeta + 1) * z / ((zeta + z^2)^2)
  }

  if (family == "LO") {
    vz <- (exp(abs(z)) - 1) / (abs(z) * (1 + exp(abs(z))))
    dvz <- sign(z) * (2 * abs(z) * exp(abs(z)) - exp(2 * abs(z)) + 1) /
      (z^2 * (1 + exp(abs(z)))^2)
  }

  if (family == "SN") {
    vz <- 4 * sinh(z) * cosh(z) / (zeta^2 * z) - tanh(z) / z
    dvz <- ((cosh(z) * z - sinh(z)) / z^2) * (4 * cosh(z) / zeta^2 - 1 / cosh(z)) +
      (sinh(z) / z) * (4 * sinh(z) / zeta^2 + sinh(z) / (cosh(z)^2))
  }

  if (family == "PE") {
    pzeta <- sqrt(2^(-2 / zeta) * gamma(1 / zeta) * (gamma(3 / zeta)^(-1)))
    vz <- (zeta * (z^2)^(zeta / 2 - 1)) / (2 * pzeta^zeta)
    dvz <- ((zeta^2 / 2 - zeta) * (pzeta^(-zeta)) * (z^2)^(zeta / 2)) / (z^3)
  }

  if (family == "Hyp") {
    vz <- zeta / sqrt(1 + z^2)
    dvz <- -(z * zeta) / (1 + z^2)^(3 / 2)
  }

  if (family == "SLASH") {
    s_aux <- z^2 / 2
    beta_aux <- zeta + (1 / 2)
    gama_aux <- zipfR::Igamma(beta_aux, s_aux)
    vz <- (2 / z^2) * zipfR::Igamma(zeta + (3 / 2), s_aux) / gama_aux
    dvz <- (-2 / z) * vz +
      (2 / z) * (1 / gama_aux^2) * exp(-s_aux) * (s_aux^beta_aux) *
        (gama_aux * (1 - beta_aux / s_aux) + exp(-s_aux) * (s_aux^(beta_aux - 1)))
  }

  list(z = z, v = vz, dv = dvz)
}
