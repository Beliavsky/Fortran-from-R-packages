#' Compute Discount Factors
#'
#' Calculate discount factors from a yield curve assuming continuous
#' compounding.
#'
#' @param curve A `yc_curve` object.
#' @param maturities Optional numeric vector of maturities. If NULL, uses
#'   the curve's own maturities.
#' @param compounding Character. Compounding convention: `"continuous"`
#'   (default), `"annual"`, or `"semi_annual"`.
#'
#' @return A data frame with columns `maturity` and `discount_factor`.
#'
#' @export
#' @examples
#' maturities <- c(1, 2, 5, 10)
#' rates <- c(0.045, 0.043, 0.042, 0.040)
#' curve <- yc_curve(maturities, rates)
#' yc_discount(curve)
#' yc_discount(curve, compounding = "annual")
yc_discount <- function(curve, maturities = NULL,
                         compounding = c("continuous", "annual",
                                          "semi_annual")) {
  validate_yc_curve(curve)
  compounding <- match.arg(compounding)

  if (is.null(maturities)) {
    maturities <- curve$maturities
    r <- if (!is.null(curve$fitted)) curve$fitted else curve$rates
  } else {
    validate_maturities(maturities)
    r <- yc_predict(curve, maturities)$rate
  }

  df <- switch(compounding,
    continuous  = exp(-r * maturities),
    annual      = (1 + r)^(-maturities),
    semi_annual = (1 + r / 2)^(-2 * maturities)
  )

  data.frame(
    maturity = maturities,
    discount_factor = df,
    stringsAsFactors = FALSE
  )
}
