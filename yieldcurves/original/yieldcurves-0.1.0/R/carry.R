#' Carry and Roll-Down Analysis
#'
#' Decompose expected return from holding a bond into carry (yield income
#' minus financing cost) and roll-down (capital gain from sliding down the
#' curve).
#'
#' @param curve A `yc_curve` object.
#' @param maturities Numeric vector of bond maturities to analyse. If NULL,
#'   uses the curve's own maturities (excluding the shortest).
#' @param horizon Numeric. Holding period in years. Default is `1/12`
#'   (one month).
#' @param funding_rate Optional numeric. Overnight funding rate as a
#'   decimal. If NULL, uses the shortest rate on the curve.
#'
#' @return A data frame with columns `maturity`, `carry`, `rolldown`, and
#'   `total`.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 1, 2, 5, 10, 30)
#' rates <- c(0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' yc_carry(fit)
yc_carry <- function(curve, maturities = NULL, horizon = 1 / 12,
                      funding_rate = NULL) {
  validate_yc_curve(curve)
  validate_positive_scalar(horizon, "horizon")

  if (is.null(maturities)) {
    maturities <- curve$maturities[curve$maturities > horizon]
  } else {
    validate_maturities(maturities)
    if (any(maturities <= horizon)) {
      cli_abort("All {.arg maturities} must be greater than {.arg horizon}.")
    }
  }

  if (length(maturities) == 0) {
    cli_abort("No maturities are greater than the horizon ({horizon} years).")
  }

  if (is.null(funding_rate)) {
    funding_rate <- yc_predict(curve, horizon)$rate
  }

  # Current yields at each maturity
  r_current <- yc_predict(curve, maturities)$rate

  # Yield at maturity - horizon (where the bond rolls to)
  r_rolled <- yc_predict(curve, maturities - horizon)$rate

  # Carry = (yield - funding rate) * horizon
  carry <- (r_current - funding_rate) * horizon

  # Roll-down = price gain from sliding down the curve
  # For a zero-coupon bond: rolldown approx = (r_current - r_rolled) * remaining_duration
  rolldown <- (r_current - r_rolled) * (maturities - horizon)

  data.frame(
    maturity = maturities,
    carry = carry,
    rolldown = rolldown,
    total = carry + rolldown,
    stringsAsFactors = FALSE
  )
}
