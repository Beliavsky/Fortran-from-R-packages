#' Bernoulli-geometric distribution
#'
#' @author Hermes Marques, \email{hermes.marques@ufrn.br}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Bernoulli-geometric distribution
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
#' Bourguignon, M., & de Medeiros, R. M. (2022). A simple and useful regression model for fitting count data. Test, 31(3), 790-827.
#'
#' @seealso \link{BerG}.
#'
#' @details
#' The BerG distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{(1-\mu+\sigma)}{(1+\mu+\sigma)}} if \eqn{x=0},
#'
#' \eqn{f(x | \mu, \sigma) = 4 \mu \frac{(\mu+\sigma-1)^{x-1}}{(\mu+\sigma+1)^{x+1}}} if \eqn{x=1, 2, ...},
#'
#' with \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\sigma>|\mu-1|}.
#'
#' @return
#' \code{dBerG} gives the density, \code{pBerG} gives the distribution
#' function, \code{qBerG} gives the quantile function, \code{rBerG}
#' generates random deviates.
#'
#' @example  examples/examples_dBerG.R
#'
#' @export
#'
dBerG <- function(x, mu, sigma, log = FALSE){
  if (any(!is.finite(mu)) || any(mu < 0))
    stop("mu must be finite and non-negative")
  if (any(!is.finite(sigma)) || any(sigma < 0))
    stop("sigma must be finite and non-negative")
  x <- as.integer(x + .5)
  if (any(x < 0))
    stop("'x' must be non-negative")
  dBerG <- ifelse(x == 0,
                  (1 - mu + sigma)/(1 + mu + sigma),
                  4*mu*(mu + sigma - 1)^(x - 1)/
                    (mu + sigma + 1)^(x + 1))

  no_cond <- sigma <= abs(mu-1)
  dBerG[no_cond] <- 0

  if(log){
    dBerG <- ifelse(dBerG < 1e-323, 1e-323, dBerG)
    log(dBerG)
  }else dBerG
}
#dBerG <- Vectorize(dBerG)
#' @export
#' @rdname dBerG
pBerG <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(!is.finite(mu)) || any(mu < 0))
    stop("mu must be finite and non-negative")
  if (any(!is.finite(sigma)) || any(sigma < 0))
    stop("sigma must be finite and non-negative")
  pBerG <- ifelse(q < 0,
                  0,
                  (1 - mu + sigma)/(1 + mu + sigma) +
                    2*mu/(1 + mu + sigma)*(1 - ((mu + sigma - 1)/
                                                  (mu + sigma + 1))^q))
  no_cond <- sigma <= abs(mu-1)
  pBerG[no_cond] <- 0

  if (!lower.tail) pBerG <- 1 - pBerG
  if(log.p)
    log(pBerG)
  else pBerG
}
#' @importFrom stats runif
#' @export
#' @rdname dBerG
rBerG <- function(n, mu, sigma){
  if (any(!is.finite(mu)) || any(mu < 0))
    stop("mu must be finite and non-negative")
  if (any(!is.finite(sigma)) || any(sigma < 0))
    stop("sigma must be finite and non-negative")
  if (!any(sigma > abs(mu - 1)))
    stop("'sigma' must be greater than 'abs(mu - 1)'")
  u <- runif(n)
  qBerG(p = u, mu = mu, sigma = sigma)
}
#' @export
#' @rdname dBerG
qBerG <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(!is.finite(mu)) || any(mu < 0))
    stop("mu must be finite and non-negative")
  if (any(!is.finite(sigma)) || any(sigma < 0))
    stop("sigma must be finite and non-negative")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  numq <- log(1 - p) + log(mu + sigma + 1) - log(2*mu)
  denq <- log(mu + sigma - 1) - log(mu + sigma + 1)
  qBerG <- ifelse(p >= (1 - mu + sigma)/(1 + mu + sigma),
                  ceiling(numq/denq), 0)
  no_cond <- sigma <= abs(mu-1)
  qBerG[no_cond] <- 0
  qBerG
}
