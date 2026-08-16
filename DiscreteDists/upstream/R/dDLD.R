#' The Discrete Lindley distribution
#'
#' @author Yojan Andrés Alcaraz Pérez, \email{yalcaraz@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete Lindley distribution
#' with parameter \eqn{\mu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of positive values of this parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Bakouch, H. S., Jazi, M. A., & Nadarajah, S. (2014). A new
#' discrete distribution. Statistics, 48(1), 200-240.
#'
#' @seealso \link{DLD}.
#'
#' @details
#' The Discrete Lindley distribution with parameters \eqn{\mu} has a support
#' 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu) = \frac{e^{-\mu x}}{1 + \mu} \left[ \mu(1 - 2e^{-\mu}) + (1- e^{-\mu})(1+\mu x)\right]}
#'
#' Note: in this implementation we changed the original parameters \eqn{\theta} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDLD} gives the density, \code{pDLD} gives the distribution
#' function, \code{qDLD} gives the quantile function, \code{rDLD}
#' generates random deviates.
#'
#' @example  examples/examples_dDLD.R
#'
#' @export
#'
dDLD <- function(x,mu, log=FALSE){
  if (any(mu <= 0))  stop("Parameter mu has to be positive!")
  if (any(x < 0))    stop(paste("x must be >=0", "\n", ""))
  res <- -mu*x - log(1+mu) + log(mu*(1-2*exp(-mu)) + (1-exp(-mu))*(1+mu*x))
  if(log)
    result <- res
  else
    result <- exp(res)
  return(result)
}
#' @export
#' @rdname dDLD
pDLD <- function(q, mu, lower.tail = TRUE, log.p=FALSE){
  if (any(mu <= 0))  stop("Parameter mu has to be positive!")
  if (any(q < 0)) stop(paste("q must be >=0", "\n", ""))
  cdf <- 1 - ((1+mu+mu*q)/(1+mu))*exp(-mu*q)

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom pracma lambertWn
#' @export
#' @rdname dDLD
qDLD <- function(p, mu, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0))  stop("Parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1))  stop("Parameter p has to be between 0 and 1")
  x <- floor(-1-(1/mu)-(1/mu)*lambertWn(-(1+mu)*exp(-1-mu)*(1-p)))
  return(x)
}
#' @importFrom pracma lambertWn
#' @export
#' @rdname dDLD
rDLD <- function(n, mu = 0.5) {
  if (any(mu <= 0)) stop("Parameter mu has to be positive!")
  U <- runif(n)
  X <- floor(-1-(1/mu)-(1/mu)*lambertWn(-(1+mu)*exp(-1-mu)*(1-U)))
  return(X)
}
