#' The hyper-Poisson distribution
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the hyper-Poisson, HYPERPO(), distribution
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
#' Sáez-Castillo, A. J., & Conde-Sánchez, A. (2013).
#' A hyper-Poisson regression model for overdispersed and
#' underdispersed count data. Computational Statistics & Data
#' Analysis, 61, 148-157.
#'
#' @seealso \link{HYPERPO}.
#'
#' @details
#' The hyper-Poisson distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu^x}{_1F_1(1;\mu;\sigma)}\frac{\Gamma(\sigma)}{\Gamma(x+\sigma)}}
#'
#' where the function \eqn{_1F_1(a;c;z)} is defined as
#'
#' \eqn{_1F_1(a;c;z) = \sum_{r=0}^{\infty}\frac{(a)_r}{(c)_r}\frac{z^r}{r!}}
#'
#' and \eqn{(a)_r = \frac{\gamma(a+r)}{\gamma(a)}} for \eqn{a>0} and \eqn{r} positive integer.
#'
#' Note: in this implementation we changed the original parameters \eqn{\lambda} and \eqn{\gamma}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dHYPERPO} gives the density, \code{pHYPERPO} gives the distribution
#' function, \code{qHYPERPO} gives the quantile function, \code{rHYPERPO}
#' generates random deviates.
#'
#' @example examples/examples_dHYPERPO.R
#'
#' @export
#' @useDynLib DiscreteDists
#' @importFrom Rcpp sourceCpp
dHYPERPO <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")

  temp <- cbind(x, mu, sigma, log)
  dHYPERPO_vec(x=temp[, 1], mu=temp[, 2], sigma=temp[, 3], log=temp[,4])
}
#' @export
#' @rdname dHYPERPO
pHYPERPO <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")

  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  # Begin auxiliar function
  aux_func <- function(q, mu, sigma) {
    cdf <- numeric(length(q))
    for (i in 1:length(q)) {
      res <- dHYPERPO(x=-1:q[i], mu=mu[i], sigma=sigma[i], log=FALSE)
      cdf[i] <- sum(res)
    }
    cdf
  }
  # End auxiliar function
  cdf <- aux_func(q=q, mu=mu, sigma=sigma)
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom stats runif
#' @export
#' @rdname dHYPERPO
rHYPERPO <- function(n, mu=1, sigma=1) {
  if (!is.numeric(n) || length(n) != 1 || n < 0)
    stop("invalid arguments")
  if (!(is.double(sigma) || is.integer(sigma)) || !(is.double(mu) || is.integer(mu)))
    stop("Non-numeric argument to mathematical function")
  sigma <- rep(sigma, length.out = n)
  mu <- rep(mu, length.out = n)
  result <- numeric(length = n)
  warn <- FALSE
  for (ind in seq_len(n)) {
    if (sigma[ind] <= 0 || mu[ind] <= 0) {
      result[ind] <- NaN
      warn <- TRUE
    }
    else {
      result[ind] <- simulate_hp(sigma[ind], mu[ind])
    }
  }
  if (warn)
    warning("NaN(s) produced: sigma and mu must be strictly positive")
  result
}
#' @export
#' @rdname dHYPERPO
qHYPERPO <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                     log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin auxiliar function
  one_quantile_hyperpo <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dHYPERPO(x=0, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dHYPERPO(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_hyperpo <- Vectorize(one_quantile_hyperpo)
  # End auxiliar function
  one_quantile_hyperpo(p=p, mu=mu, sigma=sigma)
}
