#' The Discrete Burr Hatke family
#'
#' @author Valentina Hurtado Sepulveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' The function \code{DBH()} defines the
#' Discrete Burr Hatke distribution
#' a single parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "logit" link as the default for the mu parameter. Other links are "probit" and "cloglog"'(complementary log-log)
#'
#' @references
#' El-Morshedy, M., Eliwa, M. S., & Altun, E. (2020). Discrete Burr-Hatke distribution with properties, estimation methods and regression model. IEEE access, 8, 74359-74370.
#'
#' @seealso \link{dDBH}.
#'
#' @details
#' The Discrete Burr-Hatke distribution with parameter \eqn{\mu} has a support
#' 0, 1, 2, ... and its probability mass function (pmf) is given by
#'
#' \eqn{f(x | \mu) = (\frac{1}{x+1}-\frac{\mu}{x+2})\mu^{x}}
#'
#' The pmf is log-convex for all values of \eqn{0 < \mu < 1}, where \eqn{\frac{f(x+1;\mu)}{f(x;\mu)}}
#' is an increasing function in \eqn{x} for all values of the parameter \eqn{\mu}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\lambda} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a Discrete Burr-Hatke distribution
#' in the \code{gamlss()} function.
#'
#' @example  examples/examples_DBH.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DBH <- function (mu.link="logit") {
  mstats <- checklink("mu.link", "DBH",
                      substitute(mu.link), c("logit", "probit", "cloglog", "cauchit"))

  structure(list(family=c("DBH", "Burr Hatke"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Discrete",
                 mu.link    = as.character(substitute(mu.link)),
                 mu.linkfun = mstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 mu.dr      = mstats$mu.eta,

                 # Primeras derivadas, se usaron las deriv manuales

                 dldm = function(y, mu) {
                   dldm <- (-y^2*mu+y^2-2*y*mu+2*y-mu)/(mu*(y+2-mu*(y+1)))
                   dldm
                 },

                 # Segundas derivadas, se usaron las deriva manuales

                 d2ldm2 = function(y, mu) {
                   dldm <- (-y^2*mu+y^2-2*y*mu+2*y-mu)/(mu*(y+2-mu*(y+1)))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },

                 G.dev.incr = function(y, mu, pw = 1, ...) -2*dDBH(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pDBH", type="Discrete",
                                               ymin = 0, y = y, mu = mu)),


                 mu.initial = expression(mu <- rep(estim_mu_DBH(y), length(y)) ),

                 mu.valid = function(mu) all(0 < mu & mu < 1),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu) {
                   -log(1-mu)/mu - 1
                 },

                 variance = function(mu) {
                   media <- -log(1-mu)/mu - 1
                   varianza <- mu/2*((2*(mu-3))/(mu*(mu-1))+(6*log(1-mu))/mu^2) - media^2
                   return(varianza)
                 }

  ),
  class=c("gamlss.family", "family"))
}
