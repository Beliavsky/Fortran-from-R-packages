#' Unit Maxwell-Boltzmann distribution
#'
#' @author David Villegas Ceballos, \email{david.villegas1@udea.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Unit Maxwell-Boltzmann distribution
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
#' @seealso \link{UMB}.
#'
#' @references
#' Biçer, C., Bakouch, H. S., Biçer, H. D., Alomair, G.,
#' Hussain, T., y Almohisen, A. (2024).
#' Unit Maxwell-Boltzmann Distribution and Its Application to
#' Concentrations Pollutant Data.
#' Axioms, 13(4), 226.
#'
#' @details
#' The Unit Maxwell-Boltzmann distribution with parameter \eqn{\mu}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu) = \frac{\sqrt(2/\pi) \log^2(1/x) \exp(-\frac{\log^2(1/x)}{2\mu^2})}{\mu^3 x} }
#'
#' for \eqn{0 < x < 1} and \eqn{\mu > 0}.
#'
#' @return
#' \code{dUMB} gives the density, \code{pUMB} gives the distribution
#' function, \code{qUMB} gives the quantile function, \code{rUMB}
#' generates random deviates.
#'
#' @example examples/examples_dUMB.R
#'
#' @export
dUMB <- function(x, mu = 1, log = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")

  # log_term <- log(1/x)
  # res <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(mu) - log(x) - (log_term^2) / (2 * mu^2)

  res <- ifelse(x <= 0 | x >= 1, -99999,
                0.5*log(2/pi)+2*log(log(1/x))-3*log(mu)-log(x)-((log(1/x))^2)/(2*mu^2))

  if(log)
    return(res)
  else
    return(exp(res))
}
#' @export
#' @rdname dUMB
pUMB <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) {
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  # Begin aux fun
  auxfun <- function(q) {
    log_term <- log(1 / q)
    erf_arg <- log_term / (sqrt(2) * mu)
    erf_part <- erf(erf_arg)
    term1 <- (sqrt(log_term^2) * erf_part) / log_term
    term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * mu^2))) / mu
    res <- 1 - term1 + term2
    res
  }
  # End aux fun

  cdf <- ifelse(q<=0, 0,
                ifelse(q >=1, 1, auxfun(q)))

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  return(cdf)
}
#' @export
#' @importFrom stats uniroot
#' @rdname dUMB
qUMB <- function(p, mu, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0 | p > 1)) stop("p must be in the interval (0, 1)")
  if (any(mu < 0))        stop("mu must be positive")

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  # Begin auxiliary function
  aux_fun <- function(x, p, mu){
    res <- p - pUMB(x, mu)
    return(res)
  }
  # End auxiliary function

  res <- uniroot(aux_fun, interval=c(0, 1),
                 p=p, mu=mu)$root

  return(res)
}
qUMB <- Vectorize(qUMB)
#' @export
#' @importFrom stats runif
#' @rdname dUMB
rUMB <- function(n = 1, mu = 1) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  u <- runif(n)
  if (length(mu) == 1) mu <- rep(mu, n)
  sapply(1:n, function(i) qUMB(u[i],mu[i]))
}
