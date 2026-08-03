#' The generalized normal distribution
#'
#' @name gnorm
#' @description Density, distribution function and random generation for the generalized normal/exponential power distribution. \cr \cr
#' A generalized normal random variable \eqn{x} with parameters \eqn{\mu}, \eqn{\alpha > 0} and \eqn{\beta > 0} has density:\cr
#' \deqn{p(x) = \beta exp{-(|x - \mu|/\alpha)^\beta}/(2\alpha \Gamma(1/\beta)).} \cr
#' The mean and variance of \eqn{x} are \eqn{\mu} and \eqn{\alpha^2 \Gamma(3/\beta)/\Gamma(1/\beta)}, respectively.
#' @aliases dgnorm pgnorm qgnorm rgnorm
#' @usage
#' dgnorm(x, mu = 0, alpha = 1, beta = 1, log = FALSE)
#' pgnorm(q, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#' qgnorm(p, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#' rgnorm(n, mu = 0, alpha = 1, beta = 1)
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu location parameter
#' @param alpha scale parameter
#' @param beta shape parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise \eqn{P[X> x]}
#' @source \code{dgnorm}, \code{pgnorm}, \code{qgnorm} and\code{rgnorm} are all parametrized as in Version 1 of the \href{https://en.wikipedia.org/wiki/Generalized_normal_distribution}{Generalized Normal Distribution Wikipedia page},
#' which uses the parametrization given by in Nadarajah (2005).
#' The same distribution was described much earlier by Subbotin (1923) and named the exponential power distribution by Box and Tiao (1973). \cr \cr
#' Box, G. E. P. and G. C. Tiao. "Bayesian inference in Statistical Analysis." Addison-Wesley Pub. Co., Reading, Mass (1973). \cr
#' Nadarajah, Saralees. "A generalized normal distribution." Journal of Applied Statistics 32.7 (2005): 685-694. \cr
#' Subbotin, M. T. "On the Law of Frequency of Error." Matematicheskii Sbornik 31.2 (1923):  206-301.
#' @importFrom stats pgamma qgamma rbinom runif
#' @examples
#' # Plot generalized normal/exponential power density
#' # that corresponds to the standard normal density
#' xs <- seq(-1, 1, length.out = 100)
#' plot(xs, dgnorm(xs, mu = 0, alpha = sqrt(2), beta = 2), type = "l",
#'      xlab = "x", ylab = expression(p(x)))
#'
#' # Plot the generalized normal/exponential power CDF
#' # that corresponds to the standard normal CDF
#' s <- seq(-1, 1, length.out = 100)
#' plot(xs, pgnorm(xs, 0, sqrt(2), 2), type = "l", xlab = "q",
#'      ylab = expression(paste("Pr(", x<=q, ")", sep = "")))
#'
#' # Plot the generalized normal/exponential power inverse CDF
#' # that corresponds to the standard normal inverse CDF
#' xs <- seq(0, 1, length.out = 100)
#' plot(xs, qgnorm(xs, 0, sqrt(2), 2), type = "l", xlab = "p",
#'      ylab = expression(paste("q: p = Pr(", x<=q, ")", sep = "")))
#'
#' # Make a histogram of draws from the generalized normal/exponential
#' # power distribution that corresponds to a standard normal distribution
#' xs <- rgnorm(100, 0, sqrt(2), 2)
# 'hist(xs, xlab = "x", freq = FALSE, main = "Histogram of Draws")
#'
#' @export
dgnorm <- function(x, mu = 0, alpha = 1, beta = 1,
                   log = FALSE) {
    # A failsafe for NaN / NAs of alpha / beta
    if(any(is.nan(alpha))){
        alpha[is.nan(alpha)] <- 0
    }
    if(any(is.na(alpha))){
        alpha[is.na(alpha)] <- 0
    }
    if(any(alpha<0)){
        alpha[alpha<0] <- 0
    }
    if(any(is.nan(beta))){
        beta[is.nan(beta)] <- 0
    }
    if(any(is.na(beta))){
        beta[is.na(beta)] <- 0
    }
    gnormValues <- (exp(-(abs(x-mu)/ alpha)^beta)* beta/(2*alpha*gamma(1/beta)))
    if(log){
        gnormValues[] <- log(gnormValues)
    }
    
    return(gnormValues)
}

#' @export
pgnorm <- function(q, mu = 0, alpha = 1, beta = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0
  }
  if(any(alpha<0)){
    alpha[alpha<0] <- 0
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0
  }
  
  p <- (1/2 + sign(q - mu[])* pgamma(abs(q - mu)^beta, shape = 1/beta, rate = (1/alpha)^beta)/2)
  if (lower.tail) {
    if (!log.p) {
      return(p)
    } else {
      return(log(p))
    }
  } else if (!lower.tail) {
    if (!log.p) {
      return(1 - p)
    } else {
      return(log(1 - p))
    }
  }
}

#' @export
qgnorm <- function(p, mu = 0, alpha = 1, beta = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0
  }
  if(any(alpha<0)){
    alpha[alpha<0] <- 0
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0
  }
  
  if (lower.tail & !log.p) {
    p <- p
  } else if (lower.tail & log.p) {
    p <- exp(p)
  } else if (!lower.tail & !log.p) {
    p <- 1 - p
  } else {
    p <- log(1 - p)
  }
  
  lambda <- (1/alpha)^beta
  gnormValues <- (sign(p-0.5)*qgamma(abs(p - 0.5)*2, shape = 1/beta, scale = 1/lambda)^(1/beta) + mu)
  
  return(gnormValues)
}

#' @export
rgnorm <- function(n, mu = 0, alpha = 1, beta = 1) {
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0
  }
  if(any(alpha<0)){
    alpha[alpha<0] <- 0
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0
  }
  
  lambda <- (1/alpha)^beta
  gnormValues <- qgamma(runif(n), shape = 1/beta, scale = alpha^beta)^(1/beta)*((-1)^rbinom(n, 1, 0.5)) + mu
  return(gnormValues)
}

