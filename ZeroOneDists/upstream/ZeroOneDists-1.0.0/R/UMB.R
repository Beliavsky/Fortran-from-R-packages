#' Unit Maxwell-Boltzmann family
#'
#' @author David Villegas Ceballos, \email{david.villegas1@udea.edu.co}
#'
#' @description
#' The function \code{UMB()} defines the Unit Maxwell-Boltzmann distribution, a one parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu.
#'
#' @seealso \link{dUMB}
#'
#' @details
#' The Unit Maxwell-Boltzmann distribution with parameter \eqn{\mu}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu) = \frac{\sqrt(2/\pi) \log^2(1/x) \exp(-\frac{\log^2(1/x)}{2\mu^2})}{\mu^3 x} }
#'
#' for \eqn{0 < x < 1} and \eqn{\mu > 0}.
#'
#' @returns
#' Returns a gamlss.family object which can be used to fit a
#' UMB distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_UMB.R
#'
#' @references
#' Biçer, C., Bakouch, H. S., Biçer, H. D., Alomair, G.,
#' Hussain, T., y Almohisen, A. (2024).
#' Unit Maxwell-Boltzmann Distribution and Its Application to
#' Concentrations Pollutant Data.
#' Axioms, 13(4), 226.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @importFrom stats dnorm pnorm
#' @export
UMB <- function (mu.link="log") {

  mstats <- checklink("mu.link", "Unit Maxwell-Boltzmann",
                      substitute(mu.link), c("log", "inverse", "own"))

  structure(list(family=c("UMB", "Unit Maxwell-Boltzmann"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Continuous",

                 mu.link    = as.character(substitute(mu.link)),
                 mu.linkfun    = mstats$linkfun,
                 mu.linkinv    = mstats$linkinv,
                 mu.dr    = mstats$mu.eta,

                 # Primeras derivadas
                 dldm = function(y, mu) {
                   L <- log(1 / y)
                   dldd <- (-3 / mu + (L^2) / mu^3)
                   dldd
                 },

                 # Segundas derivadas
                 d2ldm2 = function(y, mu) {
                   L <- log(1 / y)
                   d2ldm2 <- 3 / mu^2 - (3 * L^2) / mu^4
                   return(d2ldm2)
                 },

                 G.dev.incr = function(y, mu, ...) -2*dUMB(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pUMB", type="Continuous", y=y, mu=mu)),

                 mu.initial = expression(mu <- rep(sqrt(sum(log(y)^2)/(3*length(y))), length(y)) ),
                 mu.valid    = function(mu)    all(mu > 0),
                 y.valid = function(y) all(y > 0) && all(y < 1)

  ),
  class=c("gamlss.family", "family"))
}

