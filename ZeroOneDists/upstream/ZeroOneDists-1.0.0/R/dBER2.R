#' Beta Rectangular distribution version 2
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Beta Rectangular distribution
#' with parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}
#' reparameterized to ensure \eqn{E(X)=\mu}.
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
#' @seealso \link{BER2}.
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
#' \code{dBER2} gives the density, \code{pBER2} gives the distribution
#' function, \code{qBER2} gives the quantile function, \code{rBER2}
#' generates random deviates.
#'
#' @example examples/examples_dBER2.R
#'
#' @export
dBER2 <- function(x, mu, sigma, nu, log = FALSE) {
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma <= 0))      stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  par <- gamma_alpha_2_mu_theta(mu, nu)
  res <- dBER(x=x, mu=par[, 1], sigma=sigma, nu=par[, 2])

  if (log)
    res <- log(res)
  return(res)
}
#' @export
#' @rdname dBER2
pBER2 <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma <= 0))      stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  par <- gamma_alpha_2_mu_theta(mu, nu)
  cdf <- pBER(q=q, mu=par[, 1], sigma=sigma, nu=par[, 2])

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  return(cdf)
}
#' @export
#' @rdname dBER2
qBER2 <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0 | p > 1))   stop("x must be in the interval (0, 1)")
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma <= 0))      stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  par <- gamma_alpha_2_mu_theta(mu, nu)
  res <- qBER(p=p, mu=par[, 1], sigma=sigma, nu=par[, 2])

  return(res)
}
#' @export
#' @importFrom stats runif
#' @rdname dBER2
rBER2 <- function(n, mu, sigma, nu) {
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval (0, 1)")
  if (any(sigma < 0))       stop("sigma must be positive")
  if (any(nu < 0 | nu > 1)) stop("nu must be in the interval [0, 1]")

  par <- gamma_alpha_2_mu_theta(mu, nu)
  u <- runif(n)
  res <- qBER(p=u, mu=par[, 1], sigma=sigma, nu=par[, 2])
  return(res)
}
gamma_alpha_2_mu_theta <- function(gamma, alpha) {
  # Compute theta from gamma and alpha
  theta <- alpha * (1 - abs(2 * gamma - 1))

  # Compute mu using the expression involving theta
  numerator <- gamma - 0.5 * theta
  denominator <- 1 - theta
  mu <- numerator / denominator

  # Return a named vector
  return(cbind(mu, theta))
}
