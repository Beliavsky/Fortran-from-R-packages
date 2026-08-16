#' The DGEII distribution
#'
#' @author Valentina Hurtado Sepúlveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' The function \code{DGEII()} defines the
#' Discrete generalized exponential distribution of a Second type
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "logit" link as the default for the mu parameter. Other links are "probit" and "cloglog"'(complementary log-log).
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Nekoukhou, V., Alamatsaz, M. H., & Bidram, H. (2013). Discrete generalized exponential distribution of a second type. Statistics, 47(4), 876-887.
#'
#' @seealso \link{dDGEII}.
#'
#' @details
#' The DGEII distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = (1-\mu^{x+1})^{\sigma}-(1-\mu^x)^{\sigma}}
#'
#' with \eqn{0 < \mu < 1} and \eqn{\sigma > 0}. If \eqn{\sigma=1}, the DGEII distribution
#' reduces to the geometric distribution with success probability \eqn{1-\mu}.
#'
#' Note: in this implementation we changed the original parameters
#' \eqn{p} to \eqn{\mu} and \eqn{\alpha} to \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a DGEII distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_DGEII.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DGEII <- function (mu.link="logit", sigma.link="log") {

  mstats <- checklink("mu.link", "DGEII",
                      substitute(mu.link), c("logit", "probit", "cloglog", "cauchit"))
  dstats <- checklink("sigma.link", "DGEII",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("DGEII", "Discrete generalized exponential distribution of a second type II"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Discrete",

                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),

                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,

                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,

                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,

                 # Primeras derivadas por ahora computacionales

                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },


                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Segundas derivadas por ahora computacionales

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))

                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))

                   dd   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dDGEII(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pDGEII", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_DGEII(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_DGEII(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(0 < mu & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}
