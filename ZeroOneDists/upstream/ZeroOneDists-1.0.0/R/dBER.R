#' Beta Rectangular distribution
#'
#' @author Karina Maria Garay, \email{kgarayo@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Beta Rectangular distribution
#' with parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param nu vector of the nu parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @seealso \link{BER}.
#'
#' @references
#' Bayes, C. L., Bazán, J. L., & García, C. (2012). A new robust
#' regression model for proportions. Bayesian Analysis, 7(4), 841-866.
#'
#' @details
#' The Beta Rectangular distribution with parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu, \sigma, \nu) = \nu + (1 - \nu) b(x| \mu, \sigma)}
#'
#' for \eqn{0 < x < 1}, \eqn{0 < \mu < 1}, \eqn{\sigma > 0} and \eqn{0 < \nu < 1}.
#' The function \eqn{b(.)} corresponds to the traditional beta distribution
#' that can be computed by \code{dbeta(x, shape1=mu*sigma, shape2=(1-mu)*sigma)}.
#'
#' @return
#' \code{dBER} gives the density, \code{pBER} gives the distribution
#' function, \code{qBER} gives the quantile function, \code{rBER}
#' generates random deviates.
#'
#' @example examples/examples_dBER.R
#'
#' @export
#' @importFrom stats dbeta
dBER <- function(x, mu, sigma, nu, log=FALSE) {
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma < 0))       stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  B_value <- dbeta(x, shape1=mu*sigma, shape2=(1-mu)*sigma)
  res <- ifelse(x<0 | x>1, 0, nu + (1 - nu) * B_value)

  if(log)
    return(log(res))
  else
    return(res)
}
#' @export
#' @importFrom stats pbeta
#' @rdname dBER
pBER <- function(q, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma < 0))       stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  B_value <- pbeta(q=q, shape1=mu*sigma, shape2=(1-mu)*sigma)
  cdf <- ifelse(q<0, 0,
                ifelse(q>1, 1, nu * q + (1 - nu) * B_value))

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  return(cdf)
}
#' @export
#' @importFrom stats pbeta uniroot
#' @rdname dBER
qBER <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0 | p > 1))   stop("p must be in the interval (0, 1)")
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma < 0))       stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  aux_fun <- function(x, p, mu, sigma, nu){
    res <- p-nu*x-(1-nu)*pbeta(q=x, shape1=mu*sigma, shape2=(1-mu)*sigma)
    return(res)
  }
  res <- uniroot(aux_fun, interval=c(0, 1),
                 p=p, mu=mu, sigma=sigma, nu=nu)$root

  return(res)
}
qBER <- Vectorize(qBER)
#' @export
#' @importFrom stats runif
#' @rdname dBER
rBER <- function(n, mu, sigma, nu){
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma < 0))       stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")
  u <- runif(n)
  res <- qBER(p=u, mu=mu, sigma=sigma, nu=nu)
  return(res)
}
