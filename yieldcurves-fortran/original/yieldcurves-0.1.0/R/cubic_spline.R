#' Fit Cubic Spline Yield Curve
#'
#' Fit a yield curve using cubic spline interpolation. Provides an exact
#' fit through all observed data points with smooth interpolation between
#' them.
#'
#' @param maturities Numeric vector of maturities in years.
#' @param rates Numeric vector of observed yields as decimals.
#' @param method Character. Spline method: `"natural"` (default) or
#'   `"fmm"` (Forsythe, Malcolm, and Moler).
#' @param type Character. Rate type: `"zero"` (default), `"par"`, or
#'   `"forward"`.
#' @param date Optional Date for the curve.
#'
#' @return A `yc_curve` object with `method = "cubic_spline"`.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_cubic_spline(maturities, rates)
#' fit
yc_cubic_spline <- function(maturities, rates,
                             method = c("natural", "fmm"),
                             type = c("zero", "par", "forward"),
                             date = NULL) {
  method_spline <- match.arg(method)
  type <- match.arg(type)
  validate_maturities(maturities)
  validate_rates(rates, length(maturities))

  if (length(maturities) < 3) {
    cli_abort("Cubic spline fitting requires at least 3 data points.")
  }

  if (!is.null(date) && !inherits(date, "Date")) {
    cli_abort("{.arg date} must be a {.cls Date} or NULL.")
  }

  # Sort by maturity
  ord <- order(maturities)
  maturities <- maturities[ord]
  rates <- rates[ord]

  # Fit spline
  sfun <- splinefun(maturities, rates, method = method_spline)
  fitted_rates <- sfun(maturities)

  new_yc_curve(
    maturities = maturities,
    rates = rates,
    type = type,
    method = "cubic_spline",
    params = list(spline_method = method_spline, splinefun = sfun),
    fitted = fitted_rates,
    residuals = rates - fitted_rates,
    date = date
  )
}
