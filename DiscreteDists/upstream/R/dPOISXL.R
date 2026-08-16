#' The Discrete Poisson XLindley
#'
#' @author Mariana Blandon Mejia, \email{mblandonm@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete Poisson XLindley distribution
#' with parameter \eqn{\mu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param n number of random values to return
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \code{P[X > x]}.
#'
#' @references
#' Ahsan-ul-Haq, M., Al-Bossly, A., El-Morshedy, M., & Eliwa,
#' M. S. (2022). Poisson XLindley distribution for count data:
#' statistical and reliability properties with estimation
#' techniques and inference. Computational Intelligence and
#' neuroscience, 2022(1), 6503670.
#'
#' @seealso \link{POISXL}.
#'
#' @details
#' The Discrete Poisson XLindley distribution with parameters \eqn{\mu} has a support
#' 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu) = \frac{\mu^2(x+\mu^2+3(1+\mu))}{(1+\mu)^{4+x}}}; with \eqn{\mu>0}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\alpha} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dPOISXL} gives the density, \code{pPOISXL} gives the distribution
#' function, \code{qPOISXL} gives the quantile function, \code{rPOISXL}
#' generates random deviates.
#'
#' @example  examples/examples_dPOISXL.R
#'
#' @export
dPOISXL <- function(x, mu=0.3, log=FALSE) {
  if (any(x < 0))    stop("x must be positive")
  if (any(mu <= 0))  stop("parameter mu has to be positive!")
  res <- 2*log(mu)+log(x+mu^2+3*(1+mu))-(4+x)*log(1+mu)
  if(log)
    return(res)
  else
    return(exp(res))
}
#' @export
#' @rdname dPOISXL
pPOISXL <- function(q, mu=0.3, lower.tail=TRUE, log.p=FALSE) {
  if (any(q < 0))    stop("x must be positive")
  if (any(mu <= 0))  stop("parameter mu has to be positive!")

  cdf <- 1 - (1 + mu*(q + 4 + mu * (3 + mu))) / (1 + mu)^(4 + q)

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dPOISXL
qPOISXL <- function(p, mu=0.3, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu <= 0))  stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001)) stop("p must be between 0 and 1")

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin auxiliar function
  one_quantile_POISXL <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dPOISXL(x=0, mu=mu, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dPOISXL(x=i, mu=mu, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_POISXL <- Vectorize(one_quantile_POISXL)
  # End auxiliar function
  one_quantile_POISXL(p=p, mu=mu)
}
#' @export
#' @rdname dPOISXL
rPOISXL <- function(n, mu=0.3) {
  if (any(mu <= 0))  stop("parameter mu has to be positive!")
  if (n <= 0)        stop("n must be a positive integer!")

  # Begin auxiliar function
  one_random_POISXL <- function(u, mu) {
    p <- dPOISXL(x=0, mu=mu, log=FALSE)
    F <- p
    i <- 0
    while (u >= F) {
      i <- i + 1
      p <- dPOISXL(x=i, mu=mu, log=FALSE)
      F <- F + p
    }
    return(i)
  }
  one_random_POISXL <- Vectorize(one_random_POISXL)
  # End auxiliar function
  one_random_POISXL(u=runif(n), mu)
}
