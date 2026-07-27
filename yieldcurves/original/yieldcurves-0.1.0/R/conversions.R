#' Convert Par Rates to Zero Rates
#'
#' Bootstrap zero (spot) rates from par (coupon) rates using iterative
#' stripping.
#'
#' @param maturities Numeric vector of maturities in years (must be
#'   positive integers or half-years).
#' @param par_rates Numeric vector of par rates as decimals.
#' @param frequency Integer. Coupon frequency per year: 1 for annual
#'   (default) or 2 for semi-annual.
#'
#' @return A data frame with columns `maturity` and `zero_rate`.
#'
#' @export
#' @examples
#' maturities <- c(1, 2, 3, 5, 10)
#' par_rates <- c(0.040, 0.042, 0.043, 0.044, 0.045)
#' yc_par_to_zero(maturities, par_rates)
#'
#' # Semi-annual coupons
#' yc_par_to_zero(c(0.5, 1, 2), c(0.04, 0.042, 0.043), frequency = 2)
yc_par_to_zero <- function(maturities, par_rates, frequency = 1) {
  validate_maturities(maturities)
  validate_rates(par_rates, length(maturities), "par_rates")

  if (!frequency %in% c(1, 2)) {
    cli_abort("{.arg frequency} must be 1 (annual) or 2 (semi-annual).")
  }

  # Sort by maturity
  ord <- order(maturities)
  maturities <- maturities[ord]
  par_rates <- par_rates[ord]

  step <- 1 / frequency

  if (frequency == 1 && any(maturities != floor(maturities) & maturities > 1)) {
    cli_warn("Non-integer maturities > 1 may produce inaccurate results. The bootstrap assumes annual coupon bonds with integer coupon dates.")
  }

  n <- length(maturities)
  zero_rates <- numeric(n)

  # First maturity: zero = par
  zero_rates[1] <- par_rates[1]

  for (i in seq_len(n)[-1]) {
    m <- maturities[i]
    c_rate <- par_rates[i]
    coupon <- c_rate / frequency

    # Generate coupon dates at intervals of step up to m
    coupon_dates <- seq(step, m, by = step)
    # Handle floating point: snap final date to m
    if (length(coupon_dates) > 0 && abs(coupon_dates[length(coupon_dates)] - m) < 1e-10) {
      coupon_dates[length(coupon_dates)] <- m
    }

    if (length(coupon_dates) <= 1) {
      # Only one payment: zero = par
      zero_rates[i] <- par_rates[i]
      next
    }

    # Interpolate zero rates at intermediate coupon dates
    intermediate <- coupon_dates[-length(coupon_dates)]
    if (i == 2) {
      z_interp <- approx(maturities[1:i], zero_rates[1:i],
                          xout = intermediate, rule = 2)$y
    } else {
      z_interp <- approx(maturities[1:(i - 1)], zero_rates[1:(i - 1)],
                          xout = intermediate, rule = 2)$y
    }

    # PV of intermediate coupons
    pv_coupons <- sum(coupon / (1 + z_interp / frequency)^(intermediate * frequency))

    # Solve for zero rate at final maturity
    # 1 = pv_coupons + (coupon + 1) / (1 + z_m/freq)^(m*freq)
    pv_final <- 1 - pv_coupons
    periods <- m * frequency
    zero_rates[i] <- frequency * (((1 + coupon) / pv_final)^(1 / periods) - 1)
  }

  data.frame(
    maturity = maturities,
    zero_rate = zero_rates,
    stringsAsFactors = FALSE
  )
}

#' Convert Zero Rates to Par Rates
#'
#' Compute par (coupon) rates from zero (spot) rates. The par rate for
#' maturity T is the coupon rate that makes a bond price equal to par.
#'
#' @param maturities Numeric vector of maturities in years.
#' @param zero_rates Numeric vector of zero rates as decimals.
#' @param frequency Integer. Coupon frequency per year: 1 for annual
#'   (default) or 2 for semi-annual.
#'
#' @return A data frame with columns `maturity` and `par_rate`.
#'
#' @export
#' @examples
#' maturities <- c(1, 2, 3, 5, 10)
#' zero_rates <- c(0.040, 0.042, 0.043, 0.044, 0.045)
#' yc_zero_to_par(maturities, zero_rates)
#'
#' # Semi-annual coupons
#' yc_zero_to_par(c(0.5, 1, 2), c(0.04, 0.042, 0.043), frequency = 2)
yc_zero_to_par <- function(maturities, zero_rates, frequency = 1) {
  validate_maturities(maturities)
  validate_rates(zero_rates, length(maturities), "zero_rates")

  if (!frequency %in% c(1, 2)) {
    cli_abort("{.arg frequency} must be 1 (annual) or 2 (semi-annual).")
  }

  ord <- order(maturities)
  maturities <- maturities[ord]
  zero_rates <- zero_rates[ord]

  step <- 1 / frequency

  if (frequency == 1 && any(maturities != floor(maturities) & maturities > 1)) {
    cli_warn("Non-integer maturities > 1 may produce inaccurate results. The par rate formula assumes annual coupon bonds with integer coupon dates.")
  }

  n <- length(maturities)
  par_rates <- numeric(n)

  for (i in seq_len(n)) {
    m <- maturities[i]
    coupon_dates <- seq(step, m, by = step)
    # Snap final date
    if (length(coupon_dates) > 0 && abs(coupon_dates[length(coupon_dates)] - m) < 1e-10) {
      coupon_dates[length(coupon_dates)] <- m
    }

    if (length(coupon_dates) <= 1) {
      # Single payment: par = zero
      par_rates[i] <- zero_rates[i]
      next
    }

    # Interpolate zero rates at all coupon dates
    z_at_coupons <- approx(maturities[1:i], zero_rates[1:i],
                            xout = coupon_dates, rule = 2)$y

    # Discount factors at coupon dates
    df <- (1 + z_at_coupons / frequency)^(-coupon_dates * frequency)

    # Par rate = frequency * (1 - df_T) / sum(df)
    par_rates[i] <- frequency * (1 - df[length(df)]) / sum(df)
  }

  data.frame(
    maturity = maturities,
    par_rate = par_rates,
    stringsAsFactors = FALSE
  )
}
