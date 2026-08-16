#' The Discrete Lindley family
#'
#' @author Yojan Andrés Alcaraz Pérez, \email{yalcaraz@unal.edu.co}
#'
#' @description
#' The function \code{DLD()} defines the
#' Discrete Lindley distribution
#' a single parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#'
#' @references
#' Bakouch, H. S., Jazi, M. A., & Nadarajah, S. (2014). A new discrete distribution. Statistics, 48(1), 200-240.
#'
#' @seealso \link{dDLD}.
#'
#' @details
#' The Discrete Lindley distribution with parameters \eqn{\mu > 0} has a support
#' 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu) = \frac{e^{-\mu x}}{1 + \mu}  (\mu(1 - 2e^{-\mu}) + (1- e^{-\mu})(1+\mu x))}
#'
#' The parameter \eqn{\mu} can be interpreted as a strict upper bound on the failure rate function
#'
#' The conventional discrete distributions (such as geometric, Poisson, etc.) are not
#' suitable for various scenarios like reliability, failure times, and counts. Consequently,
#' alternative discrete distributions have been created by adapting well-known continuous
#' models for reliability and failure times. Among these, the discrete Weibull distribution
#' stands out as the most widely used. But models like these require two parameters and not many
#' of the known discrete distributions can provide accurate models for both times and counts,
#' which the Discrete Lindley distribution does.
#'
#' Note: in this implementation we changed the original parameters \eqn{\theta} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a Discrete Lindley distribution
#' in the \code{gamlss()} function.
#'
#' @example  examples/examples_DLD.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DLD <- function (mu.link="log") {
  mstats <- checklink("mu.link", "DLD",
                      substitute(mu.link), c("log", "identity"))

  structure(list(family=c("DLD", "Discrete Lindley"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Discrete",
                 mu.link    = as.character(substitute(mu.link)),
                 mu.linkfun = mstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 mu.dr      = mstats$mu.eta,

                  # Primera derivada manual:
                 dldm_manual= function(y, mu) {
                   p1 <- -y
                   p2 <- -(1/(1+mu))
                   p3num <- exp(mu)-1+2*mu+y*exp(mu)-y+mu*y
                   p3den <- mu*exp(mu)-2*mu+exp(mu)+y*mu*exp(mu)-1-mu*y
                   p3 <- p3num/p3den
                   dldm <- p1 + p2 + p3
                   dldm
                 },

                 # Segunda derivada manual
                 d2ldm2_manual = function(y, mu) {
                   p1 <- -y
                   p2 <- -(1/(1+mu))
                   p3num <- exp(mu)-1+2*mu+y*exp(mu)-y+mu*y
                   p3den <- mu*exp(mu)-2*mu+exp(mu)+y*mu*exp(mu)-1-mu*y
                   p3 <- p3num/p3den
                   dldm <- p1 + p2 + p3
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },

                 G.dev.incr = function(y, mu, pw = 1, ...) -2*dDLD(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pDLD", type="Discrete",
                                               ymin = 0, y = y, mu = mu)),


                 mu.initial = expression(mu <- rep(estim_mu_DLD(y), length(y)) ),

                 mu.valid = function(mu) all(0 < mu),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu) {
                   num_mean <- exp(-mu)*(2*mu-mu*exp(-mu)+1-exp(-mu))
                   den_mean <- (1+mu)*(1-exp(-mu))^2
                   media <- num_mean/den_mean
                   return(media)
                 },

                 variance = function(mu) {
                  num_var <- exp(-mu)*(mu*exp(-2*mu)-3*mu^2*exp(-mu)+3*mu+2*mu^2-2*exp(-mu)+1+exp(-2*mu)-4*mu*exp(-mu))
                  den_var <- (1+mu)^2*(1-exp(-mu))^4
                  varianza <- num_var/den_var
                  return(varianza)
                 }
  ),
  class=c("gamlss.family", "family"))
}
