#' Duration and Convexity
#'
#' Compute Macaulay duration, modified duration, and convexity for
#' zero-coupon bonds at each maturity on the curve.
#'
#' @param curve A `yc_curve` object.
#' @param maturities Optional numeric vector of maturities. If NULL,
#'   uses the curve's own maturities.
#' @param compounding Character. Compounding convention: `"continuous"`
#'   (default), `"annual"`, or `"semi_annual"`.
#'
#' @return A data frame with columns `maturity`, `macaulay_duration`,
#'   `modified_duration`, and `convexity`.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 1, 2, 5, 10, 30)
#' rates <- c(0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_nelson_siegel(maturities, rates)
#' yc_duration(fit)
yc_duration <- function(curve, maturities = NULL,
                         compounding = c("continuous", "annual",
                                         "semi_annual")) {
  validate_yc_curve(curve)
  compounding <- match.arg(compounding)

  if (is.null(maturities)) {
    maturities <- curve$maturities
  } else {
    validate_maturities(maturities)
  }

  rates <- yc_predict(curve, maturities)$rate

  if (compounding == "continuous") {
    # For continuous compounding: Macaulay duration = maturity (zero-coupon)
    # Modified duration = maturity (same as Macaulay for continuous)
    # Convexity = maturity^2
    mac_dur <- maturities
    mod_dur <- maturities
    convexity <- maturities^2
  } else if (compounding == "annual") {
    # Macaulay duration for zero-coupon = maturity
    # Modified duration = maturity / (1 + r)
    # Convexity = maturity * (maturity + 1) / (1 + r)^2
    mac_dur <- maturities
    mod_dur <- maturities / (1 + rates)
    convexity <- maturities * (maturities + 1) / (1 + rates)^2
  } else {
    # Semi-annual: Modified duration = maturity / (1 + r/2)
    # Convexity = maturity * (maturity + 0.5) / (1 + r/2)^2
    mac_dur <- maturities
    mod_dur <- maturities / (1 + rates / 2)
    convexity <- maturities * (maturities + 0.5) / (1 + rates / 2)^2
  }

  data.frame(
    maturity = maturities,
    macaulay_duration = mac_dur,
    modified_duration = mod_dur,
    convexity = convexity,
    stringsAsFactors = FALSE
  )
}

#' Coupon Bond Duration and Convexity
#'
#' Compute Macaulay duration, modified duration, and convexity for a
#' coupon-bearing bond.
#'
#' @param face Numeric. Face (par) value of the bond. Default is 100.
#' @param coupon_rate Numeric. Annual coupon rate as a decimal (e.g., 0.05
#'   for 5 percent).
#' @param maturity Numeric. Time to maturity in years.
#' @param yield Numeric. Yield to maturity as a decimal.
#' @param frequency Integer. Coupon frequency per year: 1 for annual or
#'   2 for semi-annual (default).
#' @param compounding Character. Compounding convention: `"semi_annual"`
#'   (default), `"annual"`, or `"continuous"`.
#'
#' @return A list with components `macaulay_duration`, `modified_duration`,
#'   `convexity`, and `price`.
#'
#' @export
#' @examples
#' # 2-year 5% bond at 4% yield, semi-annual coupons
#' yc_bond_duration(face = 100, coupon_rate = 0.05, maturity = 2,
#'                  yield = 0.04, frequency = 2)
yc_bond_duration <- function(face = 100, coupon_rate, maturity, yield,
                              frequency = 2,
                              compounding = c("semi_annual", "annual",
                                              "continuous")) {
  compounding <- match.arg(compounding)

  if (!is.numeric(face) || length(face) != 1 || face <= 0) {
    cli_abort("{.arg face} must be a single positive number.")
  }
  if (!is.numeric(coupon_rate) || length(coupon_rate) != 1) {
    cli_abort("{.arg coupon_rate} must be a single numeric value.")
  }
  if (!is.numeric(maturity) || length(maturity) != 1 || maturity <= 0) {
    cli_abort("{.arg maturity} must be a single positive number.")
  }
  if (!is.numeric(yield) || length(yield) != 1) {
    cli_abort("{.arg yield} must be a single numeric value.")
  }
  if (!frequency %in% c(1, 2)) {
    cli_abort("{.arg frequency} must be 1 (annual) or 2 (semi-annual).")
  }

  if (compounding == "continuous") {
    # Continuous compounding
    n_periods <- maturity * frequency
    periods <- seq_len(n_periods)
    times <- periods / frequency
    coupon_cf <- rep(face * coupon_rate / frequency, n_periods)
    coupon_cf[n_periods] <- coupon_cf[n_periods] + face

    # Price = sum CF_i * exp(-y * t_i)
    disc <- exp(-yield * times)
    price <- sum(coupon_cf * disc)

    # Macaulay = (1/P) * sum(t_i * CF_i * exp(-y*t_i))
    mac_dur <- sum(times * coupon_cf * disc) / price

    # Modified = Macaulay for continuous compounding
    mod_dur <- mac_dur

    # Convexity = (1/P) * sum(t_i^2 * CF_i * exp(-y*t_i))
    convexity <- sum(times^2 * coupon_cf * disc) / price
  } else {
    # Discrete compounding (annual or semi-annual)
    freq <- if (compounding == "annual") 1 else 2
    n_periods <- as.integer(round(maturity * freq))
    periods <- seq_len(n_periods)
    times <- periods / freq

    coupon_cf <- rep(face * coupon_rate / freq, n_periods)
    coupon_cf[n_periods] <- coupon_cf[n_periods] + face

    y_per <- yield / freq
    disc <- (1 + y_per)^(-periods)
    price <- sum(coupon_cf * disc)

    # Macaulay = (1/P) * sum(t_i * CF_i / (1+y/freq)^i)
    mac_dur <- sum(times * coupon_cf * disc) / price

    # Modified = Macaulay / (1 + y/freq)
    mod_dur <- mac_dur / (1 + y_per)

    # Convexity = (1/P) * sum(t_i * (t_i + 1/freq) * CF_i / (1+y/freq)^i) / (1+y/freq)^2
    convexity <- sum(times * (times + 1 / freq) * coupon_cf * disc) / (price * (1 + y_per)^2)
  }

  list(
    macaulay_duration = mac_dur,
    modified_duration = mod_dur,
    convexity = convexity,
    price = price
  )
}
