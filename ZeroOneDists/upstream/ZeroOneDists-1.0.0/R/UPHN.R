#' The Unit-Power Half-Normal family
#'
#' @author Juan Diego Suarez Hernandez, \email{jsuarezhe@unal.edu.co}
#'
#' @description
#' The function \code{UPHN()} defines the Unit-Power Half-Normal
#' distribution, a two parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Santoro, K. I., Gomez, Y. M., Soto, D., & Barranco-Chamorro, I. (2024).
#' Unit-Power Half-Normal Distribution Including Quantile Regression
#' with Applications to Medical Data. Axioms, 13(9), 599.
#'
#' @seealso \link{dUPHN}.
#'
#' @details
#' The UPHN distribution with parameters
#' \eqn{\mu} and \eqn{\sigma}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu, \sigma) = \frac{2\mu}{\sigma x^2} \phi(\frac{1-x}{\sigma x}) (2 \Phi(\frac{1-x}{\sigma x})-1)^{\mu-1}}
#'
#' for \eqn{0 < x < 1}, \eqn{\mu > 0} and \eqn{\sigma > 0}.
#'
#' @returns
#' Returns a \code{gamlss.family} object which can be used to fit a
#' UPHN distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_UPHN.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
UPHN <- function(mu.link = "log", sigma.link = "log") {
  mstats <- checklink("mu.link", "UPHN",
                      substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "UPHN",
                      substitute(sigma.link), c("log", "identity"))

  structure(
    list(
      family = c("UPHN", "Unit-Power Half-Normal"),
      parameters = list(mu = TRUE, sigma = TRUE),
      nopar = 2,
      type = "Continuous",
      mu.link = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),
      mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
      mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
      mu.dr = mstats$mu.eta,
      sigma.dr = dstats$mu.eta,

      # First derivates

      dldm = function(y, mu, sigma) {
        dm   <- gamlss::numeric.deriv(dUPHN(y, mu, sigma, log=TRUE),
                                      theta="mu",
                                      delta=0.01)
        dldm <- as.vector(attr(dm, "gradient"))
        dldm
      },

      dldd = function(y, mu, sigma) {
        dd   <- gamlss::numeric.deriv(dUPHN(y, mu, sigma, log=TRUE),
                                      theta="sigma",
                                      delta=0.01)
        dldd <- as.vector(attr(dd, "gradient"))
        dldd
      },

      # Second derivates

      d2ldm2 = function(y, mu, sigma) {
        dm   <- gamlss::numeric.deriv(dUPHN(y, mu, sigma, log=TRUE),
                                      theta="mu",
                                      delta=0.01)
        dldm <- as.vector(attr(dm, "gradient"))
        d2ldm2 <- - dldm * dldm
        d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
        d2ldm2
      },

      d2ldmdd = function(y, mu, sigma) {
        dm   <- gamlss::numeric.deriv(dUPHN(y, mu, sigma, log=TRUE),
                                      theta="mu",
                                      delta=0.01)
        dldm <- as.vector(attr(dm, "gradient"))
        dd   <- gamlss::numeric.deriv(dUPHN(y, mu, sigma, log=TRUE),
                                      theta="sigma",
                                      delta=0.01)
        dldd <- as.vector(attr(dd, "gradient"))
        d2ldmdd <- - dldm * dldd
        d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
        d2ldmdd
      },

      d2ldd2  = function(y, mu, sigma) {
        dd   <- gamlss::numeric.deriv(dUPHN(y, mu, sigma, log=TRUE),
                                      theta="sigma",
                                      delta=0.01)
        dldd <- as.vector(attr(dd, "gradient"))
        d2ldd2 <- - dldd * dldd
        d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
        d2ldd2
      },

      G.dev.incr = function(y, mu, sigma, ...) -2 * dUPHN(y, mu, sigma, log = TRUE),
      rqres = expression(rqres(pfun = "pUPHN", type = "Continuous", y = y, mu = mu, sigma = sigma)),

      mu.initial = expression(mu <- rep(1, length(y))),
      sigma.initial = expression(sigma <- rep(1, length(y))),

      mu.valid = function(mu) all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0),
      y.valid = function(y) all(y > 0 & y < 1)
    ),
    class = c("gamlss.family", "family")
  )
}
