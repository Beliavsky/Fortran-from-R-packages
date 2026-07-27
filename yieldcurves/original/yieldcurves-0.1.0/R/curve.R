#' Create a Yield Curve Object
#'
#' Construct a `yc_curve` object from observed maturity-rate pairs. This is
#' the core data structure used throughout the package.
#'
#' @param maturities Numeric vector of maturities in years (e.g., 0.25 for
#'   3 months, 2 for 2 years).
#' @param rates Numeric vector of yields as decimals (e.g., 0.05 for 5\%).
#'   Must be the same length as `maturities`.
#' @param type Character. The type of rate: `"zero"` (default), `"par"`, or
#'   `"forward"`.
#' @param date Optional Date for the curve observation.
#'
#' @return A `yc_curve` object (S3 class) with components:
#'   \describe{
#'     \item{maturities}{Numeric vector of maturities in years.}
#'     \item{rates}{Numeric vector of rates as decimals.}
#'     \item{type}{Character string indicating rate type.}
#'     \item{method}{Character string indicating fitting method.}
#'     \item{params}{List of model parameters (empty for observed curves).}
#'     \item{fitted}{Numeric vector of fitted rates (NULL for observed curves).}
#'     \item{residuals}{Numeric vector of residuals (NULL for observed curves).}
#'     \item{date}{Date of the curve observation.}
#'     \item{n_obs}{Integer count of maturity points.}
#'   }
#'
#' @export
#' @examples
#' # US Treasury yields (2Y, 5Y, 10Y, 30Y)
#' maturities <- c(2, 5, 10, 30)
#' rates <- c(0.045, 0.042, 0.040, 0.043)
#' curve <- yc_curve(maturities, rates)
#' curve
yc_curve <- function(maturities, rates, type = c("zero", "par", "forward"),
                     date = NULL) {
  type <- match.arg(type)
  validate_maturities(maturities)
  validate_rates(rates, length(maturities))

  if (!is.null(date) && !inherits(date, "Date")) {
    cli_abort("{.arg date} must be a {.cls Date} or NULL.")
  }

  # Sort by maturity

  ord <- order(maturities)
  maturities <- maturities[ord]
  rates <- rates[ord]

  new_yc_curve(
    maturities = maturities,
    rates = rates,
    type = type,
    method = "observed",
    params = list(),
    fitted = NULL,
    residuals = NULL,
    date = date
  )
}

#' Internal constructor for yc_curve
#' @noRd
new_yc_curve <- function(maturities, rates, type, method, params = list(),
                          fitted = NULL, residuals = NULL, date = NULL) {
  structure(
    list(
      maturities = maturities,
      rates      = rates,
      type       = type,
      method     = method,
      params     = params,
      fitted     = fitted,
      residuals  = residuals,
      date       = date,
      n_obs      = length(maturities)
    ),
    class = "yc_curve"
  )
}
