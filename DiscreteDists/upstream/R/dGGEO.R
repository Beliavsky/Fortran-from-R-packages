#' The GGEO distribution
#'
#' @author Valentina Hurtado Sepulveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Generalized Geometric distribution
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
#' Gómez-Déniz, E. (2010). Another generalization of the geometric
#' distribution. Test, 19, 399-415.
#'
#' @seealso \link{GGEO}.
#'
#' @details
#' The GGEO distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\sigma \mu^x (1-\mu)}{(1-(1-\sigma) \mu^{x+1})(1-(1-\sigma) \mu^{x})}}
#'
#' with \eqn{0 < \mu < 1} and \eqn{\sigma > 0}. If \eqn{\sigma=1}, the GGEO distribution
#' reduces to the geometric distribution with success probability \eqn{1-\mu}.
#'
#'
#' Note: in this implementation we changed the original parameters
#' \eqn{\theta} for \eqn{\mu} and \eqn{\alpha} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dGGEO} gives the density, \code{pGGEO} gives the distribution
#' function, \code{qGGEO} gives the quantile function, \code{rGGEO}
#' generates random deviates.
#'
#' @example  examples/examples_dGGEO.R
#'
#' @export
#'
dGGEO <- function(x, mu=0.5, sigma=1, log=FALSE){
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")

  res <- ifelse(x < 0,
                -Inf,
                log(sigma) + x*log(mu) + log(1-mu)-log(1-(1-sigma)*mu^(x+1))-log(1-(1-sigma)*mu^x))

  if(log){
    return(res)}
  else{
    return(exp(res))}
}
dGGEO <- Vectorize(dGGEO)
#' @export
#' @rdname dGGEO
pGGEO <- function(q, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")

  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length=ly)
  mu    <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)

  fn <- function(q, mu, sigma) sum(dGGEO(x=0:q, mu=mu, sigma=sigma))
  Vcdf <- Vectorize(fn)
  cdf <- Vcdf(q=q, mu=mu, sigma=sigma)
  cdf <- ifelse(q < 0, 0, cdf)

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom stats runif
#' @export
#' @rdname dGGEO
rGGEO <- function(n, mu=0.5, sigma=1) {
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")
  if (any(n <= 0))                  stop(paste("n must be a positive integer", "\n", ""))

  # Begin auxiliar function
  one_random_GGEO <- function(u, mu, sigma) {
    p <- dGGEO(x=0, mu=mu, sigma=sigma, log=FALSE)
    F <- p
    i <- 0
    while (u >= F) {
      i <- i + 1
      p <- dGGEO(x=i, mu=mu, sigma=sigma, log=FALSE)
      F <- F + p
    }
    return(i)
  }
  one_random_GGEO <- Vectorize(one_random_GGEO)
  # End auxiliar function

  one_random_GGEO(u=runif(n), mu, sigma)
}
#' @export
#' @rdname dGGEO
qGGEO <- function(p, mu=0.5, sigma=1, lower.tail=TRUE,
                  log.p=FALSE) {
  if (any(mu <= 0) | any(mu >= 1))  stop("parameter mu must be in (0, 1)")
  if (any(sigma <= 0))              stop("parameter sigma has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  # Begin auxiliar function
  one_quantile_GGEO <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dGGEO(x=0, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dGGEO(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_GGEO <- Vectorize(one_quantile_GGEO)
  # End auxiliar function

  one_quantile_GGEO(p=p, mu=mu, sigma=sigma)
}
