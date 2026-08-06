#' Calculate the modified duration of a bond
#'
#' @param yields a series of yields
#' @param maturity constant bond maturity in years
#' @param format_out xts or tibble
#'
#' @details
#' This function is normally used in combination with total_return() to compute bond total returns.
#'
#' The modified duration is the interest rate sensitivity of the price of bond.
#'
#' \deqn{duration = \frac{1}{y_t}{z_t}^{2M}}{duration = 1/(y_t) * z_t^(2M)}
#'
#' with
#'
#' \deqn{z_t = 1+\frac{y_t}{2}}{z_t = 1 + y_t/2}
#'
#' \eqn{M}{M} is the maturity in years (e.g. 10), \eqn{y_t}{y_t} is the yield at time \eqn{t}{t}.
#'
#' @return
#' A series of modified duration
#'
#' @export
mod_duration <- function(yields, maturity, format_out = "xts") {
  if (!(format_out %in% c("xts", "tibble"))) {
    stop("format_out has to be one of xts or tibble")
  }

  out <- (1-(1/(1+0.5*yields)^(2*maturity)))/(yields)

  if (format_out == "tibble") {
    out <- xts_to_tibble(out)
  }

  return(out)
}
