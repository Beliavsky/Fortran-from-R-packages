#' Extract Forward Rates
#'
#' Compute forward rates from a yield curve. Can compute either
#' instantaneous forward rates or forward-forward rates between two tenors.
#'
#' The instantaneous forward rate is derived as:
#' \deqn{f(m) = r(m) + m \cdot r'(m)}
#'
#' @param curve A `yc_curve` object.
#' @param maturities Optional numeric vector of maturities at which to
#'   compute forward rates. If NULL, uses the curve's own maturities.
#' @param horizon Optional numeric. If provided, computes the forward rate
#'   from each maturity to maturity + horizon (forward-forward rate).
#'
#' @return A data frame with columns `maturity` and `forward_rate`.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' yc_forward(fit)
#' yc_forward(fit, maturities = c(1, 5, 10), horizon = 1)
yc_forward <- function(curve, maturities = NULL, horizon = NULL) {
  validate_yc_curve(curve)

  if (is.null(maturities)) {
    maturities <- curve$maturities
  } else {
    validate_maturities(maturities)
  }

  if (!is.null(horizon)) {
    validate_positive_scalar(horizon, "horizon")
    # Forward-forward rate: f(t, t+h) = [r(t+h)*(t+h) - r(t)*t] / h
    r_t <- yc_predict(curve, maturities)$rate
    r_th <- yc_predict(curve, maturities + horizon)$rate
    fwd <- (r_th * (maturities + horizon) - r_t * maturities) / horizon
  } else {
    # Instantaneous forward: f(t) = r(t) + t * r'(t)
    fwd <- switch(curve$method,
      nelson_siegel = ns_instantaneous_forward(maturities, curve$params),
      svensson = sv_instantaneous_forward(maturities, curve$params),
      {
        # Numerical derivative for spline/observed curves
        r <- yc_predict(curve, maturities)$rate
        dm <- pmax(maturities * 1e-6, 1e-8)
        r_plus <- yc_predict(curve, maturities + dm)$rate
        dr <- (r_plus - r) / dm
        r + maturities * dr
      }
    )
  }

  data.frame(
    maturity = maturities,
    forward_rate = fwd,
    stringsAsFactors = FALSE
  )
}

#' Nelson-Siegel instantaneous forward rate (analytical)
#' @noRd
ns_instantaneous_forward <- function(m, params) {
  x <- m / params$tau
  params$beta0 + params$beta1 * exp(-x) +
    params$beta2 * x * exp(-x)
}

#' Svensson instantaneous forward rate (analytical)
#' @noRd
sv_instantaneous_forward <- function(m, params) {
  x1 <- m / params$tau1
  x2 <- m / params$tau2
  params$beta0 + params$beta1 * exp(-x1) +
    params$beta2 * x1 * exp(-x1) +
    params$beta3 * x2 * exp(-x2)
}
