#' Discrete generalized exponential distribution - a second type
#'
#' @author Valentina Hurtado Sepulveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete generalized exponential distribution
#' a second type with parameters \eqn{\mu} and \eqn{\sigma}.
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
#' Nekoukhou, V., Alamatsaz, M. H., & Bidram, H. (2013). Discrete generalized exponential distribution of a second type. Statistics, 47(4), 876-887.
#'
#' @seealso \link{DGEII}.
#'
#' @details
#' The DGEII distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = (1-\mu^{x+1})^{\sigma}-(1-\mu^x)^{\sigma}}
#'
#' with \eqn{0 < \mu < 1} and \eqn{\sigma > 0}. If \eqn{\sigma=1}, the DGEII distribution
#' reduces to the geometric distribution with success probability \eqn{1-\mu}.
#'
#' Note: in this implementation we changed the original parameters
#' \eqn{p} to \eqn{\mu} and \eqn{\alpha} to \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDGEII} gives the density, \code{pDGEII} gives the distribution
#' function, \code{qDGEII} gives the quantile function, \code{rDGEII}
#' generates random deviates.
#'
#' @example  examples/examples_dDGEII.R
#'
#' @export
#'
dDGEII <- function(x, mu=0.5, sigma=1.5, log=FALSE){
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")

  # Ensure same length vector
  ly    <- max(length(x), length(mu), length(sigma))
  xx    <- rep(x, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Temporal change for invalid x's
  xx[x < 0] <- 0
  xx[is.infinite(x)] <- 1
  xx[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- 2 # No integers

  # pdf in log-scale
  p <- log((1-mu^(xx+1))^sigma - (1-mu^xx)^sigma)

  # Assign -Inf for invalid x's
  p[x < 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  p[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- -Inf

  if (log == FALSE)
    p <- exp(p)
  return(p)
}
#' @export
#' @rdname dDGEII
pDGEII <- function(q, mu=0.5, sigma=1.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))             stop("parameter sigma has to be positive!")

  # Ensure same length vector
  ly    <- max(length(q), length(mu), length(sigma))
  qq    <- rep(q, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Temporal change for invalid x's
  qq[q < 0] <- 0
  qq[q == Inf] <- 0

  # For non-integer x's, the cumulative is the same as the lower integer
  qq <- as.integer(qq)

  # Auxiliary function
  fn <- function(q, mu, sigma) sum(dDGEII(x=0:q, mu=mu, sigma=sigma))
  Vec_fn <- Vectorize(fn)

  # The cumulative
  cdf <- Vec_fn(q=qq, mu=mu, sigma=sigma)

  # Assign values for invalid x's
  cdf[q < 0] <- 0
  cdf[q == Inf] <- 1

  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)

  return(cdf)
}
#' @importFrom stats runif
#' @export
#' @rdname dDGEII
rDGEII <- function(n, mu=0.5, sigma=1.5) {
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")
  if (any(n <= 0))                  stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)
  u <- runif(n=n)
  x <- qDGEII(p=u, mu=mu, sigma=sigma)
  return(x)
}
#' @export
#' @rdname dDGEII
qDGEII <- function(p, mu=0.5, sigma=1.5, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")

  if (log.p == TRUE)
    p <- exp(p)
  if (lower.tail == FALSE)
    p <- 1 - p

  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp <- rep(p, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Temporal change for invalid p's
  pp[p < 0]  <-  0.1
  pp[p > 1]  <-  0.1
  pp[p == 1] <-  0.1
  pp[p == 0] <-  0.1

  # Begin auxiliary function
  one_quantile_DGEII <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      i <- 0
      prob <- dDGEII(x=i, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      while (p >= F) {
        i <- i + 1
        prob <- dDGEII(x=i, mu=mu, sigma=sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_DGEII <- Vectorize(one_quantile_DGEII)
  # End auxiliary function

  # The quantile
  q <- one_quantile_DGEII(p=pp, mu=mu, sigma=sigma)

  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0

  return(q)
}
