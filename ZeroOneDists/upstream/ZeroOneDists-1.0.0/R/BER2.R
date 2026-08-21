#' Beta Rectangular distribution version 2
#'
#' @description
#' The Beta Rectangular family
#'
#' @param mu.link defines the mu.link, with "logit" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma parameter.
#' @param nu.link defines the nu.link, with "logit" link as the default for the nu parameter.
#'
#' @seealso \link{dBER2}
#'
#' @details
#' The Beta Rectangular distribution with parameters \code{mu},
#' \code{sigma} and \code{nu} has density given by
#'
#' \eqn{f(x| \mu, \sigma, \nu) = \nu + (1 - \nu) b(x| \mu, \sigma)}
#'
#' for \eqn{0 < x < 1}, \eqn{0 < \mu < 1}, \eqn{\sigma > 0} and \eqn{0 < \nu < 1}.
#' The function \eqn{b(.)} corresponds to the traditional beta distribution
#' that can be computed by \code{dbeta(x, shape1=mu*sigma, shape2=(1-mu)*sigma)}.
#'
#' @returns
#' Returns a gamlss.family object which can be used to fit a
#' BER2 distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_BER2.R
#'
#' @references
#' Bayes, C. L., Bazán, J. L., & García, C. (2012). A new robust regression model for proportions. Bayesian Analysis, 7(4), 841-866.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BER2 <- function (mu.link="logit", sigma.link="log", nu.link="logit"){
  mstats <- checklink("mu.link", "Beta Rectangular",
                      substitute(mu.link), c("logit", "own"))
  dstats <- checklink("sigma.link", "Beta Rectangular",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Beta Rectangular",
                      substitute(nu.link), c("logit", "own"))

  structure(list(family=c("BER2", "Beta Rectangular"),
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE),
                 nopar=3,
                 type="Continuous",

                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 nu.link    = as.character(substitute(nu.link)),

                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 nu.linkfun    = vstats$linkfun,

                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 nu.linkinv    = vstats$linkinv,

                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 nu.dr    = vstats$mu.eta,

                 # First derivates
                 dldm = function(y, mu, sigma, nu) {
                   dm   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma, nu) {
                   dd   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 dldv = function(y, mu, sigma, nu) {
                   dv   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.00001)
                   dldv <- as.vector(attr(dv, "gradient"))
                   dldv
                 },

                 # Second derivates

                 d2ldm2 = function(y, mu, sigma, nu) {
                   dm   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldd2  = function(y, mu, sigma, nu) {
                   dd   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 d2ldv2  = function(y, mu, sigma, nu) {
                   dv   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.00001)
                   dldv <- as.vector(attr(dv, "gradient"))
                   d2ldv2 <- - dldv * dldv
                   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2, -1e-15)
                   d2ldv2
                 },

                 d2ldmdd = function(y, mu, sigma, nu) {
                   dm   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldmdv = function(y, mu, sigma, nu) {
                   dm   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dv   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.00001)
                   dldv <- as.vector(attr(dv, "gradient"))

                   d2ldmdv <- - dldm * dldv
                   d2ldmdv <- ifelse(d2ldmdv < -1e-15, d2ldmdv, -1e-15)
                   d2ldmdv
                 },

                 d2ldddv = function(y, mu, sigma, nu) {
                   dd   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dv   <- gamlss::numeric.deriv(dBER2(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.00001)
                   dldv <- as.vector(attr(dv, "gradient"))

                   d2ldmdv <- - dldd * dldv
                   d2ldmdv <- ifelse(d2ldmdv < -1e-15, d2ldmdv, -1e-15)
                   d2ldmdv
                 },

                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dBER2(y, mu, sigma, nu, log=TRUE),
                 rqres      = expression(rqres(pfun="pBER2", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),

                 mu.initial    = expression(mu    <- rep(0.5, length(y))),
                 sigma.initial = expression(sigma <- rep(1.5, length(y))),
                 nu.initial    = expression(nu    <- rep(0.5, length(y))),

                 mu.valid    = function(mu)    all(mu >= 0 & mu <= 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 nu.valid    = function(nu)    all(nu >= 0 & nu <= 1),

                 mean = function(mu) mu,

                 y.valid = function(y) all(y > 0 & y < 1)

  ),
  class=c("gamlss.family", "family"))
}
