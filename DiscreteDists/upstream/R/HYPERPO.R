#' The hyper Poisson family
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' The function \code{HYPERPO()} defines the
#' hyper Poisson distribution
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
#' @seealso \link{dHYPERPO}.
#'
#' @details
#' The hyper-Poisson distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu^x}{_1F_1(1;\mu;\sigma)}\frac{\Gamma(\sigma)}{\Gamma(x+\sigma)}}
#'
#' where the function \eqn{_1F_1(a;c;z)} is defined as
#'
#' \eqn{_1F_1(a;c;z) = \sum_{r=0}^{\infty}\frac{(a)_r}{(c)_r}\frac{z^r}{r!}}
#'
#' and \eqn{(a)_r = \frac{\gamma(a+r)}{\gamma(a)}} for \eqn{a>0} and \eqn{r} positive integer.
#'
#' Note: in this implementation we changed the original parameters \eqn{\lambda} and \eqn{\gamma}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a hyper-Poisson distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_HYPERPO.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
HYPERPO <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "HYPERPO",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "HYPERPO",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("HYPERPO", "Hyper-Poisson"),
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
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))

                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   dldd
                 },

                 # Second derivatives

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))

                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))

                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dHYPERPO(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pHYPERPO", type="Discrete",
                                               ymin=0, y=y, mu=mu, sigma=sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_HYPERPO(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_HYPERPO(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}
