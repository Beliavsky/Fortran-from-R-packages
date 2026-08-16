#' The Bernoulli-geometric distribution
#'
#' @author Hermes Marques, \email{hermes.marques@ufrn.br}
#'
#' @description
#' The function \code{BerG()} defines the
#' Bernoulli-geometric distribution,
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Bourguignon, M., & de Medeiros, R. M. (2022). A simple and useful regression model for fitting count data. Test, 31(3), 790-827.
#'
#' @seealso \link{dBerG}.
#'
#' @details
#' The BerG distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{(1-\mu+\sigma)}{(1+\mu+\sigma)}} if \eqn{x=0},
#'
#' \eqn{f(x | \mu, \sigma) = 4 \mu \frac{(\mu+\sigma-1)^{x-1}}{(\mu+\sigma+1)^{x+1}}} if \eqn{x=1, 2, ...},
#'
#' with \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\sigma>|\mu-1|}.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a BerG distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_BerG.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BerG <- function (mu.link="log", sigma.link="log") {
  mstats <- checklink("mu.link", "BerG",
                      substitute(mu.link),
                      c("sqrt", "log", "identity"))
  dstats <- checklink("sigma.link", "BerG",
                      substitute(sigma.link),
                      c("sqrt", "log", "identity"))
  structure(list(family = c("BerG", "Bernoulli-geometric (BerG) distribution"),
                 parameters = list(mu=TRUE, sigma=TRUE),
                 nopar = 2,
                 type = "Discrete",

                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),

                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,

                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,

                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,

                 # First derivates
                 dldm = function(y, mu, sigma){
                   res <- ifelse(y == 0,
                                 2*(sigma + 1)/((mu - sigma - 1)*(mu + sigma + 1)),
                                 1/mu + (y - 1)/(mu + sigma - 1) - (y + 1)/(mu + sigma + 1)
                   )
                   res
                 },

                 dldd = function(y, mu, sigma){
                   res <- ifelse(y == 0,
                                 2*mu/((sigma - mu + 1)*(sigma + mu + 1)),
                                 (y - 1)/(mu + sigma - 1) - (y + 1)/(mu + sigma + 1)
                   )
                   res
                 },

                 # Second derivates
                 d2ldm2 = function(y, mu, sigma){
                   res <- ifelse(y == 0,
                                 -4*(sigma + 1)*mu/((mu - sigma - 1)^2*(mu + sigma + 1)^2),
                                 (y + 1)/(mu + sigma + 1)^2 + (1 - y)/(mu + sigma - 1)^2 - 1/mu^2
                   )
                   res
                 },

                 d2ldd2 = function(y, mu, sigma){
                   res <- ifelse(y == 0,
                                 -4*mu*(sigma + 1)/((sigma - mu + 1)^2 * (sigma + mu + 1)^2),
                                 (y + 1)/(mu + sigma + 1)^2 + (1 - y)/(mu + sigma - 1)^2
                   )
                   res
                 },
                 d2ldmdd = function(y, mu, sigma){
                   res <- ifelse(y == 0,
                                 2*((sigma + 1)^2 + mu^2)/((sigma - mu + 1)^2*(mu + sigma + 1)^2),
                                 2*((mu + sigma)*(mu + sigma - 2*y) + 1)/
                                   ((mu + sigma + 1)^2*(mu + sigma - 1)^2)
                   )
                   res
                 },

                 G.dev.incr = function(y, mu, sigma, pw=1, ...) -2*dBerG(y, mu, sigma, log=TRUE),

                 rqres = expression(rqres(pfun="pBerG", type="Discrete",
                                          ymin=0, y=y, mu=mu, sigma=sigma)),

                 mu.initial = expression({ mu <- mean(y) }),
                 sigma.initial = expression({ sigma <- var(y)/mean(y)} ),

                 mu.valid = function(mu) all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y)  all(y >= 0),
                 mean = function(mu, sigma) mu,
                 variance = function(mu, sigma) mu * sigma
  ),
  class = c("gamlss.family","family"))
}

