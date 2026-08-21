#' Unit Half Logistic-Geometry distribution
#'
#' @author Juan Diego Suarez Hernandez, \email{jsuarezhe@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Unit Half Logistic-Geometry distribution
#' with parameter \eqn{\mu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @seealso \link{UHLG}.
#'
#' @references
#' Ramadan, A. T., Tolba, A. H., & El-Desouky, B. S. (2022).
#' A unit half-logistic geometric distribution and its application
#' in insurance. Axioms, 11(12), 676.
#'
#' @details
#' The Unit Half Logistic-Geometry distribution with parameter \eqn{\mu}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu) = \frac{2 \mu}{(\mu+(2-\mu)x)^2} }
#'
#' for \eqn{0 < x < 1} and \eqn{\mu > 0}.
#'
#' @return
#' \code{dUHLG} gives the density, \code{pUHLG} gives the distribution
#' function, \code{qUHLG} gives the quantile function, \code{rUHLG}
#' generates random deviates.
#'
#' @example examples/examples_dUHLG.R
#'
#' @export
dUHLG <- function(x, mu, log=FALSE){
  if (any(mu <= 0)) stop("mu must be in the interval (0, 1)")
  if (any(x <= 0 | x >= 1)) stop("x must be in the interval (0,1)")

  z <- log(2 * mu)
  w <- 2 * log(mu + (2 - mu) * x)
  res <- z - w
  if (log)
    return(res)
  else
    return(exp(res))
}
#' @export
#' @rdname dUHLG
pUHLG <- function(q, mu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(q <= 0 | q >= 1)) stop("q must be in the interval (0,1)")

  z <- mu * (1 - q)
  w <- mu + (2 - mu) * q
  cdf <- 1 - (z / w)

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  if (log.p) {
    cdf <- log(cdf)
  }
  return(cdf)
}
#' @export
#' @rdname dUHLG
qUHLG <- function(p, mu, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(p <= 0 | p >= 1)) stop("q must be in the interval (0,1)")

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  q <- (p*mu)/(2-2*p+p*mu)
  return(q)
}
#' @export
#' @rdname dUHLG
rUHLG <- function(n, mu) {
  if (any(mu <= 0)) stop("mu must be positive!")

  u <- runif(n)
  x <- qUHLG(u, mu)
  return(x)
}
