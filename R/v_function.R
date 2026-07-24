# ─────────────────────────────────────────────────────────────
# v_function: Self-sufficient score helpers (z, v, v')
#
# Arguments:
#   y      : response vector, values in (0,1)
#   mu     : fitted mean vector, values in (0,1)
#   sigma  : dispersion — scalar or vector (one per observation)
#   lambda : power parameter — scalar
#   family : character flag — one of "NO", "TF", "LO", "SN",
#            "PE", "Hyp", "SLASH"
#   zeta   : extra shape parameter — scalar (default = 2)
#            used by TF, SN, PE, Hyp, SLASH
#
# Returns: list with named elements
#   z  — standardised latent residual (vector, length = length(y))
#   v  — weight function v(z)        (vector)
#   dv — derivative v'(z)            (vector)
#
# Notes:
#   * When lambda = 0 the Power-Logit collapses to the Log-Log
#     model; z is computed via the limit form (advisor's version).
#   * When lambda != 0, z uses the logit-link of y^lambda and
#     mu^lambda (VGAM::logitlink).
#
# Requirements: VGAM, zipfR
# ─────────────────────────────────────────────────────────────

v_function <- function(y, mu, sigma, lambda, family, zeta = 2) {
  
  # ── Input validation ────────────────────────────────────────
  valid_families <- c("NO", "TF", "LO", "SN", "PE", "Hyp", "SLASH")
  if (!(family %in% valid_families)) {
    stop("Error: 'family' must be one of: NO, TF, LO, SN, PE, Hyp, SLASH.")
  }
  
  # ── Compute z ───────────────────────────────────────────────
  if (lambda[1] == 0) {
    z <- (-1 / sigma) * (log(log(y) / log(mu)))
  } else {
    z <- (1 / sigma) * (VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
  }
  
  # ── Compute v(z) and v'(z) by family ────────────────────────
  if (family == "NO") {
    vz  <- rep(1, length(y))
    dvz <- rep(0, length(y))
  }
  
  if (family == "TF") {
    vz  <- (zeta + 1) / (zeta + z^2)
    dvz <- -2 * (zeta + 1) * z / ((zeta + z^2)^2)
  }
  
  if (family == "LO") {
    vz  <- (exp(abs(z)) - 1) / (abs(z) * (1 + exp(abs(z))))
    dvz <- sign(z) * (2 * abs(z) * exp(abs(z)) - exp(2 * abs(z)) + 1) /
      (z^2 * (1 + exp(abs(z)))^2)
  }
  
  if (family == "SN") {
    vz  <- 4 * sinh(z) * cosh(z) / (zeta^2 * z) - tanh(z) / z
    dvz <- ((cosh(z) * z - sinh(z)) / z^2) * (4 * cosh(z) / zeta^2 - 1 / cosh(z)) +
      (sinh(z) / z) * (4 * sinh(z) / zeta^2 + sinh(z) / (cosh(z)^2))
  }
  
  if (family == "PE") {
    pzeta <- sqrt(2^(-2 / zeta) * gamma(1 / zeta) * (gamma(3 / zeta)^(-1)))
    vz    <- (zeta * (z^2)^(zeta / 2 - 1)) / (2 * pzeta^zeta)
    dvz   <- ((zeta^2 / 2 - zeta) * (pzeta^(-zeta)) * (z^2)^(zeta / 2)) / (z^3)
  }
  
  if (family == "Hyp") {
    vz  <- zeta / sqrt(1 + z^2)
    dvz <- -(z * zeta) / (1 + z^2)^(3 / 2)
  }
  
  if (family == "SLASH") {
    s_aux    <- z^2 / 2
    beta_aux <- zeta + (1 / 2)
    gama_aux <- zipfR::Igamma(beta_aux, s_aux)
    vz  <- (2 / z^2) * zipfR::Igamma(zeta + (3 / 2), s_aux) / gama_aux
    dvz <- (-2 / z) * vz +
      (2 / z) * (1 / gama_aux^2) * exp(-s_aux) * (s_aux^beta_aux) *
      (gama_aux * (1 - beta_aux / s_aux) + exp(-s_aux) * (s_aux^(beta_aux - 1)))
  }
  
  list(z = z, v = vz, dv = dvz)
}