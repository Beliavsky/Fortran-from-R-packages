#' The Discrete Perks distribution
#'
#' @author Veronica Seguro Varela, \email{vseguro@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete Perks, DPERKS(),
#' distribution
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
#' Tyagi, A., Choudhary, N., & Singh, B. (2020). A new discrete
#' distribution: Theory and applications to discrete failure
#' lifetime and count data. J. Appl. Probab. Statist, 15, 117-143.
#'
#' @seealso \link{DPERKS}.
#'
#' @details
#' The discrete Perks distribution with parameters \eqn{\mu > 0} and \eqn{\sigma > 0}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu(1+\mu)(e^\sigma-1)e^{\sigma x}}{(1+\mu e^{\sigma x})(1+\mu e^{\sigma(x+1)})} }
#'
#' Note: in this implementation we changed the original parameters
#' \eqn{\lambda} for \eqn{\mu} and \eqn{\beta} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDPERKS} gives the density, \code{pDPERKS} gives the distribution
#' function, \code{qDPERKS} gives the quantile function, \code{rDPERKS}
#' generates random deviates.
#'
#' @example examples/examples_dDPERKS.R
#'
#' @export
#'
dDPERKS <- function(x, mu = 0.5, sigma = 0.5, log = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")

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
  part1 <- log(mu) + log(1+mu) + log(exp(sigma)-1)
  part2 <- sigma*xx - log(1+mu*exp(sigma*xx)) - log(1 + mu*exp(sigma*(xx+1)))
  p <- part1 + part2

  # Assign values for invalid x's
  p[x < 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  p[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- -Inf

  if (log == FALSE)
    p <- exp(p)
  return(p)
}
#' @export
#' @rdname dDPERKS
pDPERKS <- function(q, mu = 0.5, sigma = 0.5, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")

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

  # The cumulative
  cdf <- mu*(exp(sigma*(qq+1))- 1)/(1+mu*exp(sigma*(qq+1)))

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
#' @rdname dDPERKS
rDPERKS <- function(n, mu = 0.5, sigma = 0.5) {
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)
  u <- runif(n=n)
  x <- qDPERKS(p=u, mu=mu, sigma=sigma)
  return(x)
}
#' @export
#' @rdname dDPERKS
qDPERKS <- function(p, mu = 0.5, sigma = 0.5, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")

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
  pp[p < 0]  <- 0.1
  pp[p > 1]  <- 0.1
  pp[p == 1] <- 0.1
  pp[p == 0] <- 0.1

  # The quantile
  q <- ceiling(((1/sigma)*log((pp + mu)/(mu*(1-pp)))) - 1)

  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0

  return(q)
}
