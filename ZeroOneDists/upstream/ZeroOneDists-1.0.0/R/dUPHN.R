#' Unit-Power Half-Normal distribution
#'
#' @author Juan Diego Suarez Hernandez, \email{jsuarezhe@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Unit-Power Half-Normal distribution
#' with parameter \eqn{\mu} and \eqn{\sigma}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @seealso \link{dUPHN}.
#'
#' @references
#' Santoro, K. I., Gómez, Y. M., Soto, D., & Barranco-Chamorro, I. (2024).
#' Unit-Power Half-Normal Distribution Including Quantile Regression
#' with Applications to Medical Data. Axioms, 13(9), 599.
#'
#' @details
#' The Unit-Power Half-Normal distribution with parameters
#' \eqn{\mu} and \eqn{\sigma}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu, \sigma) = \frac{2\mu}{\sigma x^2} \phi(\frac{1-x}{\sigma x}) (2 \Phi(\frac{1-x}{\sigma x})-1)^{\mu-1}}
#'
#' for \eqn{0 < x < 1}, \eqn{\mu > 0} and \eqn{\sigma > 0}.
#'
#' @return
#' \code{dUPHN} gives the density, \code{pUPHN} gives the distribution
#' function, \code{qUPHN} gives the quantile function, \code{rUPHN}
#' generates random deviates.
#'
#' @example examples/examples_dUPHN.R
#'
#' @export
#' @importFrom stats dnorm pnorm
dUPHN <- function(x, mu, sigma, log = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  if (any(x <= 0 | x >= 1)) stop("x must be in the interval (0, 1)")

  z <- (1 - x) / (sigma * x)
  part1 <- dnorm(z, log=TRUE)
  part2 <- (mu-1)*log(2 * pnorm(z) - 1)
  log_pdf <- log(2*mu) - log(sigma) - 2 * log(x) + part1 + part2

  if (log) {
    return(log_pdf)
  } else {
    return(exp(log_pdf))
  }
}
#' @importFrom stats pnorm
#' @export
#' @rdname dUPHN
pUPHN <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")

  cdf <- ifelse(q <= 0, 0,
                ifelse(q >= 1, 1,
                       1 - (2 * pnorm((1 - q) / (sigma * q)) - 1)^mu))

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  if (log.p) {
    cdf <- log(cdf)
  }

  return(cdf)
}
#' @importFrom stats qnorm
#' @export
#' @rdname dUPHN
qUPHN <- function(p, mu, sigma,
                  lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p <= 0 | p >= 1)) stop("p must be in (0, 1)")

  inside <- ((1 - p)^(1 / mu)) / 2 + 0.5
  q <- 1 / (sigma * qnorm(inside) + 1)

  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dUPHN
rUPHN <- function(n, mu, sigma) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  u <- runif(n)
  qUPHN(u, mu, sigma)
}
