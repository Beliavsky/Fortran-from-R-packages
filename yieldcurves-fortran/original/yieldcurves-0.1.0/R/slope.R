#' Yield Curve Slope Measures
#'
#' Compute common slope and curvature measures from a yield curve.
#'
#' @param curve A `yc_curve` object.
#'
#' @return A named list with slope measures:
#'   \describe{
#'     \item{spread_2s10s}{10-year minus 2-year rate (the most common slope
#'       measure).}
#'     \item{spread_2s30s}{30-year minus 2-year rate.}
#'     \item{spread_5s30s}{30-year minus 5-year rate.}
#'     \item{spread_3m10y}{10-year minus 3-month rate (term premium proxy).}
#'     \item{butterfly_2s5s10s}{2 * 5-year minus 2-year minus 10-year
#'       (curvature measure).}
#'   }
#'   Returns NA for any measure whose required tenors fall outside the
#'   curve range.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' yc_slope(fit)
yc_slope <- function(curve) {
  validate_yc_curve(curve)

  min_m <- min(curve$maturities)
  max_m <- max(curve$maturities)

  # Helper to safely get a rate, returning NA if out of range
  safe_rate <- function(m) {
    if (m < min_m || m > max_m) return(NA_real_)
    yc_predict(curve, m)$rate
  }

  r_3m  <- safe_rate(0.25)
  r_2y  <- safe_rate(2)
  r_5y  <- safe_rate(5)
  r_10y <- safe_rate(10)
  r_30y <- safe_rate(30)

  list(
    spread_2s10s     = r_10y - r_2y,
    spread_2s30s     = r_30y - r_2y,
    spread_5s30s     = r_30y - r_5y,
    spread_3m10y     = r_10y - r_3m,
    butterfly_2s5s10s = 2 * r_5y - r_2y - r_10y
  )
}

#' Extract Level, Slope, and Curvature Factors
#'
#' For Nelson-Siegel or Svensson curves, extracts the estimated factors
#' directly from the model parameters. For other curves, computes empirical
#' measures.
#'
#' @param curve A `yc_curve` object.
#'
#' @return A named list with:
#'   \describe{
#'     \item{level}{Long-run level (beta0 for NS/Svensson, or mean rate).}
#'     \item{slope}{Slope factor (beta1 for NS/Svensson, or short - long
#'       rate).}
#'     \item{curvature}{Curvature factor (beta2 for NS/Svensson, or
#'       2*mid - short - long rate).}
#'   }
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' yc_level_slope_curvature(fit)
yc_level_slope_curvature <- function(curve) {
  validate_yc_curve(curve)

  if (curve$method %in% c("nelson_siegel", "svensson")) {
    list(
      level     = curve$params$beta0,
      slope     = curve$params$beta1,
      curvature = curve$params$beta2
    )
  } else {
    # Empirical measures
    r <- if (!is.null(curve$fitted)) curve$fitted else curve$rates
    list(
      level     = mean(r),
      slope     = r[1] - r[length(r)],
      curvature = 2 * r[ceiling(length(r) / 2)] - r[1] - r[length(r)]
    )
  }
}
