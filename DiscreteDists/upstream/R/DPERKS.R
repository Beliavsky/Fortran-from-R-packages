#' The Discrete Perks family
#'
#' @author Veronica Seguro Varela, \email{vseguro@unal.edu.co}
#'
#' @description
#' The function \code{DPERKS()} defines the
#' Discrete Perks distribution
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Tyagi, A., Choudhary, N., & Singh, B. (2020). A new discrete
#' distribution: Theory and applications to discrete failure
#' lifetime and count data. J. Appl. Probab. Statist, 15, 117-143.
#'
#' @seealso \link{dDPERKS}.
#'
#' @details
#' The discrete Perks distribution with parameters \eqn{\mu > 0} and \eqn{\sigma > 0}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu(1+\mu)(e^\sigma-1)e^{\sigma x}}{(1+\mu e^{\sigma x})(1+\mu e^{\sigma(x+1)})} }
#'
#' Note: in this implementation we changed the original parameters
#' \eqn{\lambda} for \eqn{\mu} and \eqn{\beta} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#'
#' @example examples/examples_DPERKS.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DPERKS <- function(mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "DPERKS",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "DPERKS",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("DPERKS", "discrete-Perks"),
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
                   dldm <- 1/mu + 1/(mu+1) - exp(sigma*y)/(1+mu*exp(sigma*y)) - exp(sigma*(y+1))/(1+mu*exp(sigma*(y+1)))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dldd <- exp(sigma)/(exp(sigma)-1) + y - mu*y*exp(sigma*y)/(1+mu*exp(sigma*y)) - mu*(y+1)*exp(sigma*(y+1))/(1+mu*exp(sigma*(y+1)))
                   dldd
                 },

                 # Second derivatives

                 d2ldm2 = function(y, mu, sigma) {
                   dldm <- 1/mu + 1/(mu+1) - exp(sigma*y)/(1+mu*exp(sigma*y)) - exp(sigma*(y+1))/(1+mu*exp(sigma*(y+1)))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },
                 d2ldmdd = function(y, mu, sigma) {
                   dldm <- 1/mu + 1/(mu+1) - exp(sigma*y)/(1+mu*exp(sigma*y)) - exp(sigma*(y+1))/(1+mu*exp(sigma*(y+1)))

                   dldd <- exp(sigma)/(exp(sigma)-1) + y - mu*y*exp(sigma*y)/(1+mu*exp(sigma*y)) - mu*(y+1)*exp(sigma*(y+1))/(1+mu*exp(sigma*(y+1)))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },
                 d2ldd2  = function(y, mu, sigma) {
                   dldd <- exp(sigma)/(exp(sigma)-1) + y - mu*y*exp(sigma*y)/(1+mu*exp(sigma*y)) - mu*(y+1)*exp(sigma*(y+1))/(1+mu*exp(sigma*(y+1)))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dDPERKS(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pDPERKS", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 #mu.initial    = expression(mu    <- rep(estim_mu_sigma_DPERKS(y)[1], length(y)) ),
                 #sigma.initial = expression(sigma <- rep(estim_mu_sigma_DPERKS(y)[2], length(y)) ),

                 mu.initial    = expression(mu    <- rep(1, length(y)) ),
                 sigma.initial = expression(sigma <- rep(1, length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}
#'
#' Initial values for DPERKS
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_DPERKS <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_DPERKS,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' logLik function for DPERKS
#' @description Calculates logLik for DPERKS distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @export
logLik_DPERKS <- function(logparam=c(0, 0), x){
  return(sum(dDPERKS(x,
                     mu    = exp(logparam[1]),
                     sigma = exp(logparam[2]),
                     log=TRUE)))
}

