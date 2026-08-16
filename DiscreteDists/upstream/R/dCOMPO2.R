#' The COMPO2 distribution (with mu as mean)
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the
#' Comway-Maxwell-Poisson distribution
#' with parameters \eqn{\mu} and \eqn{\sigma}.
#' This parameterization was
#' proposed by Ribeiro et al. (2020) and the main
#' characteristic is that \eqn{E(X)=\mu}.
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
#' Ribeiro Jr, Eduardo E., et al. "Reparametrization of COM–Poisson regression
#' models with applications in the analysis of experimental data."
#' Statistical Modelling 20.5 (2020): 443-466.
#'
#' @seealso \link{COMPO2}.
#'
#' @details
#' The COMPO2 distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \left(\mu + \frac{\exp(\sigma)-1}{2 \exp(\sigma)} \right)^{x \exp(\sigma)} \frac{(x!)^{\exp(\sigma)}}{Z(\mu, \sigma)} }
#'
#' with \eqn{\mu > 0}, \eqn{\sigma \in \Re} and
#'
#' \eqn{Z(\mu, \sigma)=\sum_{j=0}^{\infty} \frac{\mu^j}{(j!)^\sigma}}.
#'
#' The proposed functions here are based on the functions from
#' the COMPoissonReg package.
#'
#' @return
#' \code{dCOMPO2} gives the density, \code{pCOMPO2} gives the distribution
#' function, \code{qCOMPO2} gives the quantile function, \code{rCOMPO2}
#' generates random deviates.
#'
#' @example examples/examples_dCOMPO2.R
#'
#' @importFrom COMPoissonReg dcmp
#' @export
dCOMPO2 <- function(x, mu, sigma, log=FALSE) {
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(x < 0))   stop(paste("x must be >=0", "\n", ""))

  par <- mu_phi_2_lambda_nu_COMPO2(mu, sigma)
  # res <- dcmp(x, lambda=par$lambda, nu=par$nu)
  #
  # if (log)
  #   res <- log(res)
  #
  # return(res)

  temp <- cbind(x, par$lambda, par$nu, log)
  dCOMPO_vec(x=temp[, 1], mu=temp[, 2], sigma=temp[, 3], log=temp[,4])
}
#' @importFrom COMPoissonReg pcmp
#' @export
#' @rdname dCOMPO2
pCOMPO2 <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(q < 0))   stop(paste("q must be >=0", "\n", ""))

  par <- mu_phi_2_lambda_nu_COMPO2(mu, sigma)
  cdf <- pcmp(q, lambda=par$lambda, nu=par$nu)

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
pCOMPO2 <- Vectorize(pCOMPO2)
#' @importFrom COMPoissonReg qcmp
#' @export
#' @rdname dCOMPO2
qCOMPO2 <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  par <- mu_phi_2_lambda_nu_COMPO2(mu, sigma)
  qcmp(p, lambda=par$lambda, nu=par$nu)
}
qCOMPO2 <- Vectorize(qCOMPO2)
#' @importFrom COMPoissonReg rcmp
#' @export
#' @rdname dCOMPO2
rCOMPO2 <- function(n, mu, sigma) {
  par <- mu_phi_2_lambda_nu_COMPO2(mu, sigma)
  rcmp(n, lambda=par$lambda, nu=par$nu)
}

mu_phi_2_lambda_nu_COMPO2 <- function(mu, phi) {
  nu <- exp(phi)
  lambda <- (mu + (nu-1)/(2*nu))^nu
  list("lambda" = lambda, "nu" = nu)
}
