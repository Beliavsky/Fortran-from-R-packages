#' Fit a Yield Curve
#'
#' Unified interface for fitting a yield curve using different methods.
#' Dispatches to [yc_nelson_siegel()], [yc_svensson()], or
#' [yc_cubic_spline()].
#'
#' @param maturities Numeric vector of maturities in years.
#' @param rates Numeric vector of observed yields as decimals.
#' @param method Character. Fitting method: `"nelson_siegel"` (default),
#'   `"svensson"`, or `"cubic_spline"`.
#' @param type Character. Rate type: `"zero"` (default), `"par"`, or
#'   `"forward"`.
#' @param date Optional Date for the curve.
#' @param ... Additional arguments passed to the fitting function.
#'
#' @return A `yc_curve` object.
#'
#' @export
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
#' rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)
#' fit <- yc_fit(maturities, rates, method = "nelson_siegel")
#' fit
yc_fit <- function(maturities, rates,
                    method = c("nelson_siegel", "svensson", "cubic_spline"),
                    type = c("zero", "par", "forward"),
                    date = NULL, ...) {
  method <- match.arg(method)
  type <- match.arg(type)

  switch(method,
    nelson_siegel = yc_nelson_siegel(maturities, rates, type = type,
                                      date = date, ...),
    svensson      = yc_svensson(maturities, rates, type = type,
                                 date = date, ...),
    cubic_spline  = yc_cubic_spline(maturities, rates, type = type,
                                     date = date, ...)
  )
}
