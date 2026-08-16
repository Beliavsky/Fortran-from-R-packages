#' The hyper Poisson family (with mu as mean)
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' The function \code{HYPERPO2()} defines the
#' hyper Poisson distribution (with mu as mean)
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Sáez-Castillo, A. J., & Conde-Sánchez, A. (2013). A hyper-Poisson regression
#' model for overdispersed and underdispersed count data.
#' Computational Statistics & Data Analysis, 61, 148-157.
#'
#' @seealso \link{dHYPERPO2}, \link{HYPERPO}.
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
#' Returns a \code{gamlss.family} object which can be used
#' to fit a hyper-Poisson distribution version 2
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_HYPERPO2.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
HYPERPO2 <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "HYPERPO2",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "HYPERPO2",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("HYPERPO2", "Hyper-Poisson-2"),
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
                   dm   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Second derivatives

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dHYPERPO2(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pHYPERPO2", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(mean(y), length(y))),
                 sigma.initial = expression(sigma <- rep(1, length(y))),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}
