#' The GGEO family
#'
#' @author Valentina Hurtado Sepúlveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' The function \code{GGEO()} defines the
#' Generalized Geometric distribution
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "logit" link as the default for the sigma. Other links are "probit" and "cloglog"'(complementary log-log)
#'
#' @references
#' Gómez-Déniz, E. (2010). Another generalization of the geometric distribution. Test, 19, 399-415.
#'
#' @seealso \link{dGGEO}.
#'
#' @details
#' The GGEO distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\sigma \mu^x (1-\mu)}{(1-(1-\sigma) \mu^{x+1})(1-(1-\sigma) \mu^{x})}}
#'
#' with \eqn{0 < \mu < 1} and \eqn{\sigma > 0}. If \eqn{\sigma=1}, the GGEO distribution
#' reduces to the geometric distribution with success probability \eqn{1-\mu}.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a GGEO distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_GGEO.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GGEO <- function (mu.link="logit", sigma.link="log") {

  mstats <- checklink("mu.link", "GGEO",
                      substitute(mu.link), c("logit", "probit", "cloglog", "cauchit"))
  dstats <- checklink("sigma.link", "GGEO",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("GGEO", "Generalized geometric"),
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

                 # Primeras derivadas manuales

                 dldm = function(y, mu, sigma) {
                   part1 <- (y+1)*(sigma-1)*mu^y / ((sigma-1)*mu^(y+1)+1)
                   part2 <-     y*(sigma-1)*mu^(y-1) / ((sigma-1)*mu^y+1)
                   dldm <- y/mu - 1/(1-mu) - part1 - part2
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dldd <- 1/sigma-(mu^(y+1))/(1-(mu^(y+1))*(1-sigma))-(mu^y)/(1-(mu^y)*(1-sigma))
                   dldd
                 },

                 # Segundas derivadas manuales

                 d2ldm2 = function(y, mu, sigma) {
                   part1 <- (y+1)*(sigma-1)*mu^y / ((sigma-1)*mu^(y+1)+1)
                   part2 <-     y*(sigma-1)*mu^(y-1) / ((sigma-1)*mu^y+1)
                   dldm <- y/mu - 1/(1-mu) - part1 - part2
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },
                 d2ldmdd = function(y, mu, sigma) {
                   part1 <- (y+1)*(sigma-1)*mu^y / ((sigma-1)*mu^(y+1)+1)
                   part2 <-     y*(sigma-1)*mu^(y-1) / ((sigma-1)*mu^y+1)
                   dldm <- y/mu - 1/(1-mu) - part1 - part2

                   dldd <- 1/sigma-(mu^(y+1))/(1-(mu^(y+1))*(1-sigma))-(mu^y)/(1-(mu^y)*(1-sigma))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },
                 d2ldd2  = function(y, mu, sigma) {
                   dldd <- 1/sigma-(mu^(y+1))/(1-(mu^(y+1))*(1-sigma))-(mu^y)/(1-(mu^y)*(1-sigma))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dGGEO(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pGGEO", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_GGEO(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_GGEO(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(0 < mu & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}

