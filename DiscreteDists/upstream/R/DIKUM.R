#' The discrete Inverted Kumaraswamy family
#'
#' @author Daniel Felipe Villa Rengifo, \email{dvilla@unal.edu.co}
#'
#' @description
#' The function \code{DIKUM()} defines the
#' discrete Inverted Kumaraswamy distribution
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' El-Helbawy, A. A., Hegazy, M. A., Al-Dayian, G. R., & Abd EL-Kader, R. E. (2022). A discrete analog of the inverted Kumaraswamy distribution: properties and estimation with application to COVID-19 data. Pakistan Journal of Statistics and Operation Research, 18(1), 297-328.
#'
#' @seealso \link{dDIKUM}.
#'
#' @details
#' The discrete Inverted Kumaraswamy distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = (1-(2+x)^{-\mu})^{\sigma}-(1-(1+x)^{-\mu})^{\sigma}}
#'
#' with \eqn{\mu > 0} and \eqn{\sigma > 0}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\alpha} and \eqn{\beta}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a discrete Inverted Kumaraswamy distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_DIKUM.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DIKUM <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "DIKUM",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "DIKUM",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("DIKUM", "discrete-Inverted-Kumaraswamy"),
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

                 # Primeras derivadas, por ahora son computacionales

                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDIKUM(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDIKUM(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Segundas derivadas, por ahora son computacionales

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDIKUM(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDIKUM(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dDIKUM(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDIKUM(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dDIKUM(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pDIKUM", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_DIKUM(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_DIKUM(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu, sigma) {
                   y = 0:999
                   p1 <- (1-(2+y)^(-1*mu))^(sigma)
                   p2 <- (1-(1+y)^(-1*mu))^sigma
                   p <-  p1 - p2
                   return(sum(y*p))
                 },
                 variance = function(mu, sigma) {
                   y = 0:999
                   p1 <- (1-(2+y)^(-1*mu))^(sigma)
                   p2 <- (1-(1+y)^(-1*mu))^sigma
                   p <-  p1 - p2
                   var = sum((y^2)*p)-(sum(y*p))^2
                   return(var)
                 }

  ),
  class=c("gamlss.family", "family"))
}

