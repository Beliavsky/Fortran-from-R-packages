#' Prior Probability Interval for a Normal Distribution
#'
#' Calculates the mean and variance of a normal prior probability interval.
#' Given a lower and upper bound of a prior probability interval
#' this function recovers the implied normal prior parameters.
#' Useful for specifying informative priors
#' on the steady-state parameters (elements of \eqn{\Psi}).
#'
#' @param l Numeric. The lower bound of the prior probability interval.
#' @param u Numeric. The upper bound of the prior probability interval.
#' @param interval Numeric. The prior probability mass within the interval.
#'   Default \code{0.95}, i.e. 95%.
#' @param annualized_growthrate Logical. If \code{TRUE}, treats \code{l} and \code{u} as
#'   bounds of the annualized steady-state growth rate and calculates the implied mean
#'   and variance on the corresponding quarterly/monthly scale. Useful if you are working with
#'   a variable specified as \code{100*diff(log(x))}. Default \code{FALSE}.
#' @param freq Integer. The data frequency (e.g. \code{4} for quarterly, \code{12} for monthly).
#'   Only used when \code{annualized_growthrate = TRUE}. Default \code{NULL}.
#'
#' @return A list with two elements: \code{mean} and \code{var}, giving the
#'   mean and variance of the implied normal distribution.
#' @export
#' 
#' @details
#' Consider a CPI variable \code{CPI <- data$CPI}, observed at quarterly frequency.
#' In the model, quarter-on-quarter inflation is used \code{x <- 100*diff(log(CPI))}.
#' Suppose the prior belief is that annualized steady-state inflation lies between 1.7 and 2.3
#' with 95% probability (mean 2). On the quarterly scale used by \code{x}, this corresponds to a
#' 95% interval of 0.425 to 0.575 (mean 0.5). Since it is typically more natural to elicit a prior
#' on the annualized scale, the \code{annualized_growthrate} argument performs this conversion.
#' See \code{vignette("Homoscedastic-steady-state-BVAR")} for usage in practice.
#' 
#'
#' @examples
#' ppi(l = 1.7, u = 2.3, interval = 0.95)
#' ppi(l = 1.7, u = 2.3, interval = 0.95, annualized_growthrate = TRUE, freq = 4)
ppi <- function(l, u, interval = 0.95, annualized_growthrate = FALSE, freq = NULL) {
  
  if (l >= u) {
    stop("lower bound must be less than upper bound")
  }
  
  if (interval <= 0 || interval >= 1) {
    stop("interval must be between 0 and 1")
  }
  
  if (annualized_growthrate) {
    if (is.null(freq)) {
      stop("freq must be specified when annualized_growthrate = TRUE")
    }
    if (freq <= 0) {
      stop("freq must be positive")
    }
  }
  
  alpha <- 1 - interval
  z     <- qnorm(1 - alpha / 2)
  
  if (!annualized_growthrate) {
    mu    <- (u + l) / 2
    sigma <- (u - l) / (2 * z)
  } else {
    mu    <- (u + l) / 2 / freq
    sigma <- (u - l) / (2 * z) / freq
  }
  
  return(list(mean = mu, var = sigma^2))
}