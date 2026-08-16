#' The Discrete Burr Hatke distribution
#'
#' @author Valentina Hurtado Sepulveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete Burr Hatke distribution
#' with parameter \eqn{\mu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param n number of random values to return
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' El-Morshedy, M., Eliwa, M. S., & Altun, E. (2020). Discrete Burr-Hatke distribution with properties, estimation methods and regression model. IEEE access, 8, 74359-74370.
#'
#' @seealso \link{DBH}.
#'
#' @details
#' The Discrete Burr-Hatke distribution with parameters \eqn{\mu} has a support
#' 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu) = (\frac{1}{x+1}-\frac{\mu}{x+2})\mu^{x}}
#'
#'
#' The pmf is log-convex for all values of \eqn{0 < \mu < 1}, where \eqn{\frac{f(x+1;\mu)}{f(x;\mu)}}
#' is an increasing function in \eqn{x} for all values of the parameter \eqn{\mu}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\lambda} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDBH} gives the density, \code{pDBH} gives the distribution
#' function, \code{qDBH} gives the quantile function, \code{rDBH}
#' generates random deviates.
#'
#' @example  examples/examples_dDBH.R
#'
#' @export
#'
dDBH <- function(x, mu, log = FALSE) {
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(x < 0)) stop(paste("x must be >=0", "\n", ""))

  res <- log(1/(x+1)-mu/(x+2)) + x * log(mu)

  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}
#' @export
#' @rdname dDBH
pDBH <- function(q, mu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))
    stop(("parameter mu has to be positive!"))
  if (any(q < 0))
    stop(paste("q must be >=0", "\n", ""))

  cdf <- 1 - (mu^(q+1)/(q+2))

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dDBH
qDBH <- function(p, mu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))  stop("parameter sigma has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))
  ## Auxiliar function
  one_quantile_DBH<- function(p, mu) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dDBH(x=0, mu=mu, log = FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dDBH(x=i, mu, log = FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_DBH <- Vectorize(one_quantile_DBH)
  ## End auxiliar function
  one_quantile_DBH(p=p, mu=mu)
}
#' @export
#' @rdname dDBH
rDBH <- function(n, mu=1) {
  if (any(mu <= 0))
    stop("parameter mu must be positive!")
  if (any(n <= 0))
    stop(paste("x must be a positive integer", "\n", ""))
  ## Begin auxiliar function
  one_random_DBH <- function(u, mu) {
    p <- dDBH(x=0, mu=mu, log = FALSE)
    F <- p
    i <- 0
    while (u >= F) {
      i <- i + 1
      p <- dDBH(x=i, mu=mu, log = FALSE)
      F <- F + p
    }
    return(i)
  }
  one_random_DBH <- Vectorize(one_random_DBH)
  ## End auxiliar function
  one_random_DBH(u=runif(n), mu)
}

