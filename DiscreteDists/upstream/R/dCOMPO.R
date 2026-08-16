#' The COMPO distribution
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the
#' Conway-Maxwell-Poisson distribution
#' with parameters \eqn{\mu} and \eqn{\sigma}.
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
#' @references
#' Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., & Boatwright, P. (2005).
#' A useful distribution for fitting discrete data: revival of the
#' Conway–Maxwell–Poisson distribution. Journal of the Royal Statistical
#' Society Series C: Applied Statistics, 54(1), 127-142.
#'
#' @seealso \link{COMPO}.
#'
#' @details
#' The COMPO distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu^x}{(x!)^{\sigma} Z(\mu, \sigma)} }
#'
#' with \eqn{\mu > 0}, \eqn{\sigma \geq 0} and
#'
#' \eqn{Z(\mu, \sigma)=\sum_{j=0}^{\infty} \frac{\mu^j}{(j!)^\sigma}}.
#'
#' The proposed functions here are based on the functions from
#' the COMPoissonReg package.
#'
#' @return
#' \code{dCOMPO} gives the density, \code{pCOMPO} gives the distribution
#' function, \code{qCOMPO} gives the quantile function, \code{rCOMPO}
#' generates random deviates.
#'
#' @example examples/examples_dCOMPO.R
#'
#' @importFrom COMPoissonReg dcmp
#' @export
dCOMPO <- function(x, mu, sigma, log=FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma < 0))  stop("parameter sigma has to be no negative!")

  temp <- cbind(x, mu, sigma, log)
  dCOMPO_vec(x=temp[, 1], mu=temp[, 2], sigma=temp[, 3], log=temp[,4])
}
#' @importFrom COMPoissonReg pcmp
#' @export
#' @rdname dCOMPO
pCOMPO <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma < 0))  stop("parameter sigma has to be no negative!")
  if (any(q < 0))      stop(paste("q must be >=0", "\n", ""))

  cdf <- pcmp(q, lambda=mu, nu=sigma)

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
pCOMPO <- Vectorize(pCOMPO)
#' @importFrom COMPoissonReg qcmp
#' @export
#' @rdname dCOMPO
qCOMPO <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma < 0))  stop("parameter sigma has to be no negative!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  qcmp(p, lambda=mu, nu=sigma)
}
qCOMPO <- Vectorize(qCOMPO)
#' @importFrom COMPoissonReg rcmp
#' @export
#' @rdname dCOMPO
rCOMPO <- function(n, mu, sigma) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma < 0))  stop("parameter sigma has to be no negative!")

  rcmp(n, lambda=mu, nu=sigma)
}

