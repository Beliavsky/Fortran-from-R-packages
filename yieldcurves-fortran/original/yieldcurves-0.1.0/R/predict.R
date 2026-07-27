#' Predict Rates from a Fitted Yield Curve
#'
#' Evaluate a fitted yield curve at new maturities.
#'
#' @param curve A `yc_curve` object from a fitting function.
#' @param maturities Numeric vector of maturities at which to predict rates.
#'
#' @return A data frame with columns `maturity` and `rate`.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' yc_predict(fit, c(3, 7, 15, 20))
yc_predict <- function(curve, maturities) {
  validate_yc_curve(curve)
  validate_maturities(maturities)

  predicted <- switch(curve$method,
    nelson_siegel = ns_rate(maturities, curve$params$beta0, curve$params$beta1,
                            curve$params$beta2, curve$params$tau),
    svensson = sv_rate(maturities, curve$params$beta0, curve$params$beta1,
                        curve$params$beta2, curve$params$beta3,
                        curve$params$tau1, curve$params$tau2),
    cubic_spline = curve$params$splinefun(maturities),
    observed = {
      # Linear interpolation for observed curves
      approx(curve$maturities, curve$rates, xout = maturities, rule = 2)$y
    }
  )

  data.frame(
    maturity = maturities,
    rate = predicted,
    stringsAsFactors = FALSE
  )
}
