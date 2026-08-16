#' The hyper-Poisson distribution (with mu as mean)
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the hyper-Poisson in
#' the second parameterization with parameters \eqn{\mu} (as mean) and
#' \eqn{\sigma} as the dispersion parameter.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of positive values of this parameter.
#' @param sigma vector of positive values of this parameter.
#' @param n number of random values to return
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Sáez-Castillo, A. J., & Conde-Sánchez, A. (2013).
#' A hyper-Poisson regression model for overdispersed and
#' underdispersed count data. Computational Statistics & Data
#' Analysis, 61, 148-157.
#'
#' @seealso \link{HYPERPO2}, \link{HYPERPO}.
#'
#' @details
#' The hyper-Poisson distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ...
#'
#' Note: in this implementation the parameter \eqn{\mu} is the mean
#' of the distribution and \eqn{\sigma} corresponds to
#' the dispersion parameter. If you fit a model with this parameterization,
#' the time will increase because an internal procedure to convert \eqn{\mu}
#' to \eqn{\lambda} parameter.
#'
#' @return
#' \code{dHYPERPO2} gives the density, \code{pHYPERPO2} gives the distribution
#' function, \code{qHYPERPO2} gives the quantile function, \code{rHYPERPO2}
#' generates random deviates.
#'
#' @example  examples/examples_dHYPERPO2.R
#'
#' @export
#' @useDynLib DiscreteDists
#' @importFrom Rcpp sourceCpp
dHYPERPO2 <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")

  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_vec_cpp(media=mu, gamma=sigma)

  temp <- cbind(x, mu, sigma, log)
  dHYPERPO_vec(x=temp[, 1], mu=temp[, 2], sigma=temp[, 3], log=temp[,4])
}
#' @export
#' @rdname dHYPERPO2
pHYPERPO2 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")

  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_vec_cpp(media=mu, gamma=sigma)

  pHYPERPO(q=q, mu=mu, sigma=sigma, lower.tail=lower.tail, log.p=log.p)
}
#' @importFrom stats runif
#' @export
#' @rdname dHYPERPO2
rHYPERPO2 <- function(n, mu=1, sigma=1) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))

  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_vec_cpp(media=mu, gamma=sigma)

  if (!is.numeric(n) || length(n) != 1 || n < 0)
    stop("invalid arguments")
  if (!(is.double(sigma) || is.integer(sigma)) || !(is.double(mu) ||
                                                    is.integer(mu)))
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
#' @rdname dHYPERPO2
qHYPERPO2 <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                      log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda_vec_cpp(media=mu, gamma=sigma)

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  qHYPERPO(p=p, mu=mu, sigma=sigma)
}
