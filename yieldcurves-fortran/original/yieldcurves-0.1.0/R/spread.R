#' Z-Spread
#'
#' Compute the Z-spread (zero-volatility spread) for a bond. The Z-spread
#' is the constant spread added to each zero rate on the benchmark curve
#' that makes the discounted cash flows equal the market price.
#'
#' @param price Numeric. Market price of the bond.
#' @param coupon_rate Numeric. Annual coupon rate as a decimal.
#' @param maturity Numeric. Time to maturity in years.
#' @param curve Either a `yc_curve` object or a list/data frame with
#'   components `maturities` and `rates`.
#' @param face Numeric. Face value of the bond. Default is 100.
#' @param frequency Integer. Coupon frequency per year: 1 for annual or
#'   2 for semi-annual (default).
#'
#' @return A list with components `zspread` (the Z-spread as a decimal),
#'   `price` (the input price), and `model_price` (the price implied by the
#'   curve with the Z-spread applied).
#'
#' @export
#' @examples
#' # Create a benchmark curve
#' curve <- yc_curve(c(0.5, 1, 2, 5, 10), c(0.03, 0.035, 0.04, 0.042, 0.045))
#'
#' # A bond priced below par (positive Z-spread)
#' yc_zspread(price = 95, coupon_rate = 0.04, maturity = 5,
#'            curve = curve, frequency = 2)
yc_zspread <- function(price, coupon_rate, maturity, curve,
                        face = 100, frequency = 2) {
  if (!is.numeric(price) || length(price) != 1 || price <= 0) {
    cli_abort("{.arg price} must be a single positive number.")
  }
  if (!is.numeric(coupon_rate) || length(coupon_rate) != 1) {
    cli_abort("{.arg coupon_rate} must be a single numeric value.")
  }
  if (!is.numeric(maturity) || length(maturity) != 1 || maturity <= 0) {
    cli_abort("{.arg maturity} must be a single positive number.")
  }
  if (!is.numeric(face) || length(face) != 1 || face <= 0) {
    cli_abort("{.arg face} must be a single positive number.")
  }
  if (!frequency %in% c(1, 2)) {
    cli_abort("{.arg frequency} must be 1 (annual) or 2 (semi-annual).")
  }

  # Extract maturities and rates from curve
  if (inherits(curve, "yc_curve")) {
    curve_m <- curve$maturities
    curve_r <- if (!is.null(curve$fitted)) curve$fitted else curve$rates
  } else if (is.list(curve) || is.data.frame(curve)) {
    curve_m <- curve$maturities
    curve_r <- curve$rates
    if (is.null(curve_m) || is.null(curve_r)) {
      cli_abort("{.arg curve} must have {.field maturities} and {.field rates} components.")
    }
  } else {
    cli_abort("{.arg curve} must be a {.cls yc_curve} object or a list with {.field maturities} and {.field rates}.")
  }

  # Generate cash flow times
  step <- 1 / frequency
  cf_times <- seq(step, maturity, by = step)
  if (length(cf_times) > 0 && abs(cf_times[length(cf_times)] - maturity) < 1e-10) {
    cf_times[length(cf_times)] <- maturity
  }
  n <- length(cf_times)

  # Cash flows
  coupon <- face * coupon_rate / frequency
  cfs <- rep(coupon, n)
  cfs[n] <- cfs[n] + face

  # Interpolate zero rates at cash flow times
  z_rates <- approx(curve_m, curve_r, xout = cf_times, rule = 2)$y

  # Objective: find spread s such that sum(CF_i / (1+z_i+s)^t_i) = price
  obj <- function(s) {
    disc <- (1 + z_rates + s)^(-cf_times)
    sum(cfs * disc) - price
  }

  result <- uniroot(obj, interval = c(-0.05, 0.50), tol = 1e-12,
                     extendInt = "yes")
  zs <- result$root

  # Verify model price
  disc_final <- (1 + z_rates + zs)^(-cf_times)
  model_price <- sum(cfs * disc_final)

  list(
    zspread = zs,
    price = price,
    model_price = model_price
  )
}

#' Key Rate Durations
#'
#' Compute key rate durations by bumping the yield curve at specific tenors.
#' Each bump is triangular: the full shift is applied at the key rate tenor
#' and linearly interpolated to zero at adjacent key rate tenors.
#'
#' @param coupon_rate Numeric. Annual coupon rate as a decimal.
#' @param maturity Numeric. Time to maturity in years.
#' @param curve Either a `yc_curve` object or a list/data frame with
#'   components `maturities` and `rates`.
#' @param key_rates Numeric vector of key rate tenors in years.
#'   Default is `c(1, 2, 5, 10, 30)`.
#' @param shift Numeric. Size of the rate bump in decimal (default 0.0001,
#'   i.e. 1 basis point).
#' @param face Numeric. Face value. Default is 100.
#' @param frequency Integer. Coupon frequency: 1 (annual) or 2
#'   (semi-annual, default).
#'
#' @return A data frame with columns `tenor` and `key_rate_duration`.
#'
#' @export
#' @examples
#' curve <- yc_curve(c(1, 2, 5, 10, 30), c(0.03, 0.035, 0.04, 0.042, 0.045))
#' yc_key_rate_duration(coupon_rate = 0.04, maturity = 10,
#'                      curve = curve, key_rates = c(1, 2, 5, 10, 30))
yc_key_rate_duration <- function(coupon_rate, maturity, curve,
                                  key_rates = c(1, 2, 5, 10, 30),
                                  shift = 0.0001, face = 100,
                                  frequency = 2) {
  if (!is.numeric(coupon_rate) || length(coupon_rate) != 1) {
    cli_abort("{.arg coupon_rate} must be a single numeric value.")
  }
  if (!is.numeric(maturity) || length(maturity) != 1 || maturity <= 0) {
    cli_abort("{.arg maturity} must be a single positive number.")
  }
  if (!is.numeric(key_rates) || length(key_rates) < 1) {
    cli_abort("{.arg key_rates} must be a numeric vector.")
  }
  if (!is.numeric(shift) || length(shift) != 1 || shift <= 0) {
    cli_abort("{.arg shift} must be a single positive number.")
  }
  if (!frequency %in% c(1, 2)) {
    cli_abort("{.arg frequency} must be 1 (annual) or 2 (semi-annual).")
  }

  key_rates <- sort(key_rates)

  # Extract maturities and rates from curve
  if (inherits(curve, "yc_curve")) {
    curve_m <- curve$maturities
    curve_r <- if (!is.null(curve$fitted)) curve$fitted else curve$rates
  } else if (is.list(curve) || is.data.frame(curve)) {
    curve_m <- curve$maturities
    curve_r <- curve$rates
    if (is.null(curve_m) || is.null(curve_r)) {
      cli_abort("{.arg curve} must have {.field maturities} and {.field rates} components.")
    }
  } else {
    cli_abort("{.arg curve} must be a {.cls yc_curve} object or a list with {.field maturities} and {.field rates}.")
  }

  # Generate cash flow times
  step <- 1 / frequency
  cf_times <- seq(step, maturity, by = step)
  if (length(cf_times) > 0 && abs(cf_times[length(cf_times)] - maturity) < 1e-10) {
    cf_times[length(cf_times)] <- maturity
  }
  n <- length(cf_times)

  # Cash flows
  coupon <- face * coupon_rate / frequency
  cfs <- rep(coupon, n)
  cfs[n] <- cfs[n] + face

  # Base zero rates at cash flow times
  z_base <- approx(curve_m, curve_r, xout = cf_times, rule = 2)$y

  # Helper: price bond given zero rates at cf_times
  price_bond <- function(z) {
    disc <- (1 + z)^(-cf_times)
    sum(cfs * disc)
  }

  # Compute base price from the curve
  base_price <- price_bond(z_base)

  krd <- numeric(length(key_rates))

  for (k in seq_along(key_rates)) {
    # Create triangular bump at key_rates[k]
    kr <- key_rates[k]
    # Find adjacent key rates
    kr_left <- if (k > 1) key_rates[k - 1] else 0
    kr_right <- if (k < length(key_rates)) key_rates[k + 1] else kr + (kr - kr_left)

    # Bump weights for each cash flow time
    bump <- numeric(n)
    for (j in seq_len(n)) {
      t_j <- cf_times[j]
      if (t_j >= kr_left && t_j <= kr) {
        # Left side of triangle
        if (kr > kr_left) {
          bump[j] <- (t_j - kr_left) / (kr - kr_left)
        } else {
          bump[j] <- 1
        }
      } else if (t_j > kr && t_j <= kr_right) {
        # Right side of triangle
        if (kr_right > kr) {
          bump[j] <- (kr_right - t_j) / (kr_right - kr)
        } else {
          bump[j] <- 1
        }
      }
    }

    # Compute bumped price
    z_bumped <- z_base + bump * shift
    price_up <- price_bond(z_bumped)

    # Key rate duration = -(1/P) * dP/dr
    krd[k] <- -(price_up - base_price) / (shift * base_price)
  }

  data.frame(
    tenor = key_rates,
    key_rate_duration = krd,
    stringsAsFactors = FALSE
  )
}
