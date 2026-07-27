#' Interpolate Yield Curve
#'
#' Interpolate rates at arbitrary maturities from an observed or fitted
#' yield curve.
#'
#' @param curve A `yc_curve` object.
#' @param maturities Numeric vector of maturities at which to interpolate.
#' @param method Character. Interpolation method: `"linear"` (default),
#'   `"log_linear"`, or `"cubic"`.
#'
#' @return A data frame with columns `maturity` and `rate`.
#'
#' @export
#' @examples
#' maturities <- c(1, 2, 5, 10, 30)
#' rates <- c(0.045, 0.043, 0.042, 0.040, 0.043)
#' curve <- yc_curve(maturities, rates)
#' yc_interpolate(curve, c(3, 7, 15, 20))
yc_interpolate <- function(curve, maturities,
                            method = c("linear", "log_linear", "cubic")) {
  validate_yc_curve(curve)
  validate_maturities(maturities)
  method <- match.arg(method)

  # For fitted curves, use the model directly
  if (curve$method %in% c("nelson_siegel", "svensson", "cubic_spline")) {
    return(yc_predict(curve, maturities))
  }

  # For observed curves, interpolate
  obs_m <- curve$maturities
  obs_r <- curve$rates

  interp_rates <- switch(method,
    linear = {
      approx(obs_m, obs_r, xout = maturities, rule = 2)$y
    },
    log_linear = {
      # Interpolate in log-discount-factor space
      log_df <- -obs_r * obs_m
      log_df_interp <- approx(obs_m, log_df, xout = maturities, rule = 2)$y
      -log_df_interp / maturities
    },
    cubic = {
      sfun <- splinefun(obs_m, obs_r, method = "natural")
      sfun(maturities)
    }
  )

  data.frame(
    maturity = maturities,
    rate = interp_rates,
    stringsAsFactors = FALSE
  )
}
