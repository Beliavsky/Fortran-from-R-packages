#' The COMPO2 family (with mu as mean)
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' The function \code{COMPO2()} defines the
#' Conway-Maxwell-Poisson distribution
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#' This parameterization was
#' proposed by Ribeiro et al. (2020) and the main
#' characteristic is that \eqn{E(X)=\mu}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "identity" link as the default for the sigma.
#'
#' @references
#' Ribeiro Jr, Eduardo E., et al. "Reparametrization of COM–Poisson regression
#' models with applications in the analysis of experimental data."
#' Statistical Modelling 20.5 (2020): 443-466.
#'
#' @seealso \link{dCOMPO2}.
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
#' Returns a \code{gamlss.family} object which can be used
#' to fit a COMPO2 distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_COMPO2.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
COMPO2 <- function (mu.link="log", sigma.link="identity") {

  mstats <- checklink("mu.link", "COMPO2",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "COMPO2",
                      substitute(sigma.link), c("identity"))

  structure(list(family=c("COMPO2", "Conway Maxwell Poisson second parameterization"),
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

                 # First derivatives

                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dCOMPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dCOMPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Second derivatives

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dCOMPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dCOMPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dCOMPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dCOMPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dCOMPO2(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pCOMPO2", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(mean(y), length(y)) ),
                 sigma.initial = expression(sigma <- rep(0, length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) TRUE,

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu, sigma) {mu}

  ),
  class=c("gamlss.family", "family"))
}
