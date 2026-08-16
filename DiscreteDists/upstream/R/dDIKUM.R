#' The discrete Inverted Kumaraswamy distribution
#'
#' @author Daniel Felipe Villa Rengifo, \email{dvilla@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the discrete Inverted Kumaraswamy, DIKUM(), distribution
#' with parameters \eqn{\mu} and \eqn{\sigma}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' El-Helbawy, A. A., Hegazy, M. A., Al-Dayian, G. R., & Abd EL-Kader,
#' R. E. (2022). A discrete analog of the inverted Kumaraswamy
#' distribution: properties and estimation with application to
#' COVID-19 data. Pakistan Journal of Statistics and Operation
#' Research, 18(1), 297-328.
#'
#' @seealso \link{DIKUM}.
#'
#' @details
#' The discrete Inverted Kumaraswamy distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = (1-(2+x)^{-\mu})^{\sigma}-(1-(1+x)^{-\mu})^{\sigma}}
#'
#' with \eqn{\mu > 0} and \eqn{\sigma > 0}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\alpha} and \eqn{\beta}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDIKUM} gives the density, \code{pDIKUM} gives the distribution
#' function, \code{qDIKUM} gives the quantile function, \code{rDIKUM}
#' generates random deviates.
#'
#' @example  examples/examples_dDIKUM.R
#'
#' @export
#'
dDIKUM <- function(x, mu=1, sigma=5, log=FALSE){
  if (any(x < 0))       stop("x must be positive")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  p1 <- (1-(2+x)^-mu)^sigma
  p2 <- (1-(1+x)^-mu)^sigma
  res <- p1 - p2
  if(log)
    return(log(res))
  else
    return(res)
}
#' @export
#' @rdname dDIKUM
pDIKUM <- function(q, mu=1, sigma=5,
                  lower.tail=TRUE, log.p=FALSE){
  if(any(q < 0))        stop("x must be positive")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")

  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  cdf <- (1-(2+q)^-mu)^sigma

  if (lower.tail == TRUE) {
    cdf <- cdf
  } else {cdf = 1 - cdf}

  if (log.p == FALSE) {
    cdf <- cdf
  } else {cdf <- log(cdf)}

  return(cdf)
}
#' @importFrom stats runif
#' @export
#' @rdname dDIKUM
rDIKUM <- function(n, mu=1, sigma=5) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))

  u <- runif(n)
  ceiling((1-u^(1/sigma))^(-1/mu) - 2)
}
#' @export
#' @rdname dDIKUM
qDIKUM <- function(p, mu=1, sigma=5,
                  lower.tail=TRUE, log.p=FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  # Calculate the quantile
  quant <- ceiling((1-p^(1/sigma))^(-1/mu) - 2)

  log.quant <- log(quant)

  # Return the quantile
  if(log.p){
    if(lower.tail == TRUE){
      return(log.quant)
    }else{return(1-log.quant)}
  }else{
    if(lower.tail == TRUE){
      return(quant)
    }else{return(1-quant)}
  }
}
