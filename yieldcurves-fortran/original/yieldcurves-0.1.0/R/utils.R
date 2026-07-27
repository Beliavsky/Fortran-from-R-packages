# Internal validation helpers and Nelson-Siegel factor loadings

#' Validate maturities vector
#' @noRd
validate_maturities <- function(m, arg = "maturities") {
  if (!is.numeric(m)) {
    cli_abort("{.arg {arg}} must be a numeric vector.")
  }
  if (length(m) < 1) {
    cli_abort("{.arg {arg}} must have at least one element.")
  }
  if (any(is.na(m))) {
    cli_abort("{.arg {arg}} must not contain NA values.")
  }
  if (any(m <= 0)) {
    cli_abort("{.arg {arg}} must contain only positive values.")
  }
  invisible(m)
}

#' Validate rates vector
#' @noRd
validate_rates <- function(r, n, arg = "rates") {
  if (!is.numeric(r)) {
    cli_abort("{.arg {arg}} must be a numeric vector.")
  }
  if (length(r) != n) {
    cli_abort(
      "{.arg {arg}} must have the same length as maturities ({n}), not {length(r)}."
    )
  }
  if (any(is.na(r))) {
    cli_abort("{.arg {arg}} must not contain NA values.")
  }
  invisible(r)
}

#' Validate a yc_curve object
#' @noRd
validate_yc_curve <- function(x, arg = "curve") {
  if (!inherits(x, "yc_curve")) {
    cli_abort("{.arg {arg}} must be a {.cls yc_curve} object.")
  }
  invisible(x)
}

#' Validate positive scalar
#' @noRd
validate_positive_scalar <- function(x, arg = "x") {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    cli_abort("{.arg {arg}} must be a single positive number.")
  }
  invisible(x)
}

#' Nelson-Siegel factor loadings
#'
#' Compute the three NS factor loadings for given maturities and tau.
#' @param m Numeric vector of maturities.
#' @param tau Decay parameter.
#' @return Matrix with columns: level, slope, curvature.
#' @noRd
ns_loadings <- function(m, tau) {
  x <- m / tau
  level <- rep(1, length(m))
  slope <- (1 - exp(-x)) / x
  curvature <- slope - exp(-x)
  cbind(level = level, slope = slope, curvature = curvature)
}

#' Nelson-Siegel rate computation
#' @noRd
ns_rate <- function(m, beta0, beta1, beta2, tau) {
  L <- ns_loadings(m, tau)
  as.numeric(beta0 * L[, 1] + beta1 * L[, 2] + beta2 * L[, 3])
}

#' Svensson factor loadings (adds second hump)
#' @noRd
sv_loadings <- function(m, tau1, tau2) {
  L <- ns_loadings(m, tau1)
  x2 <- m / tau2
  curvature2 <- (1 - exp(-x2)) / x2 - exp(-x2)
  cbind(L, curvature2 = curvature2)
}

#' Svensson rate computation
#' @noRd
sv_rate <- function(m, beta0, beta1, beta2, beta3, tau1, tau2) {
  L <- sv_loadings(m, tau1, tau2)
  as.numeric(beta0 * L[, 1] + beta1 * L[, 2] + beta2 * L[, 3] + beta3 * L[, 4])
}
