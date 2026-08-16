#' The DsPA distribution
#'
#' @author Maria Camila Mena Romana, \email{mamenar@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the
#' discrete power-Ailamujia distribution
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
#' Alghamdi, A. S., Ahsan-ul-Haq, M., Babar, A.,
#' Aljohani, H. M., Afify, A. Z., & Cell, Q. E. (2022).
#' The discrete power-Ailamujia distribution: properties,
#' inference, and applications. AIMS Math, 7(5), 8344-8360.
#'
#' @seealso \link{DsPA}.
#'
#' @details
#' The DsPA distribution with parameters \eqn{\mu} and \eqn{\sigma} has a support
#'  0, 1, 2, ...
#'
#' Note:in this implementation we changed the original parameters
#' \eqn{\beta} and \eqn{\lambda} for \eqn{\mu} and \eqn{\sigma} respectively,
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDsPA} gives the density, \code{pDsPA} gives the distribution
#' function, \code{qDsPA} gives the quantile function.
#'
#' @example examples/examples_dDsPA.R
#'
#' @export
dDsPA <- function(x, mu, sigma, log=FALSE) {
  if (any(x < 0))                          stop(paste("x must be >=0", "\n", ""))
  if (any(mu <= 0))                        stop("parameter mu has to be positive!")
  if (any(sigma <= 0) || any(sigma >= 1))  stop("parameter sigma must be in (0, 1)")

  res <- (sigma^x^mu) * (1-x^mu*log(sigma)) - (sigma^((x+1)^mu)) *
    (1-((x+1)^mu)*log(sigma))

  if (log)
    return(log(res))
  else
    return(res)
}
#' @export
#' @rdname dDsPA
pDsPA <- function(q, mu, sigma,
                  lower.tail=TRUE, log.p=FALSE) {
  if (any(q < 0))                          stop(paste("x must be >=0", "\n", ""))
  if (any(mu <= 0))                        stop("parameter mu has to be positive!")
  if (any(sigma <= 0) || any(sigma >= 1))  stop("parameter sigma must be in (0, 1)")

  cdf <- 1 - sigma^((q + 1)^mu) * (1 - (q + 1)^mu * log(sigma))

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  return(cdf)
}
#' @export
#' @rdname dDsPA
qDsPA <- function(p, mu=1, sigma=1,
                  lower.tail=TRUE, log.p=FALSE) {
  if (any(mu <= 0))                        stop("parameter mu has to be positive!")
  if (any(sigma <= 0) || any(sigma >= 1))  stop("parameter sigma must be in (0, 1)")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin auxiliar function
  one_quantile_DsPA <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dDsPA(x=0, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dDsPA(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_DsPA <- Vectorize(one_quantile_DsPA)
  # End auxiliar function
  one_quantile_DsPA(p=p, mu=mu, sigma=sigma)
}
#' @export
#' @rdname dDsPA
rDsPA <- function(n, mu, sigma){
  if(any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0) || any(sigma >= 1))  stop("parameter sigma must be in (0, 1)")
  u <- runif(n)

  return(qDsPA(p=u, mu=mu, sigma=sigma))
}
