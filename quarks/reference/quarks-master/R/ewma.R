#' Exponentially weighted moving average
#'
#' Estimates volatility of a return series by means of an exponentially
#' weighted moving average.
#'
#' @param x a numeric vector of asset returns
#' @param lambda decay factor for the calculation of weights; default is \code{0.94}
#'
#' @export
#'
#' @return Returns a numerical vector \code{vol} that contains the computed
#' volatility.
#'
#' @examples
#' prices <- DAX$price_close
#' returns <- diff(log(prices))
#' date <- DAX$ref_date[-1]
#' cvar <- ewma(x = returns, lambda = 0.94)
#' csig <- sqrt(cvar)
#' plot(date, csig, type = 'l',
#'      main = 'conditional standard deviations for the DAX30 return series')

ewma <- function(x, lambda = 0.94) {
    if (length(x) <= 1 || !all(!is.na(x)) || !is.numeric(x)) {
        stop("A numeric vector of length > 1 and without NAs must be passed to",
             " 'x'.")
    }
    if (length(lambda) != 1 || is.na(lambda) || !is.numeric(lambda) ||
        lambda < 0 || lambda >= 1) {
      stop("The argument 'lambda' must be a single non-NA double value with ",
           "0 < lambda < 1.")
    }
    n <- length(x)
    vol <- rep(NA, n)
    vol[1] <- stats::var(x)
    for (i in 2:n) {
        vol[i] <- lambda * vol[i - 1] + (1 - lambda) * x[i - 1]^2
    }
    vol
}
