#' The DMOLBE family
#'
#' @author Olga Usuga, \email{olga.usuga@udea.edu.co}
#'
#' @description
#' The function \code{DMOLBE()} defines the
#' Discrete Marshall-Olkin Length Biased Exponential distribution
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Aljohani, H. M., Ahsan-ul-Haq, M., Zafar, J., Almetwally, E. M.,
#' Alghamdi, A. S., Hussam, E., & Muse, A. H. (2023). Analysis of
#' Covid-19 data using discrete Marshall–Olkinin length biased exponential:
#' Bayesian and frequentist approach. Scientific Reports, 13(1), 12243.
#'
#' @seealso \link{dDMOLBE}.
#'
#' @details
#' The DMOLBE distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\sigma ((1+x/\mu)\exp(-x/\mu)-(1+(x+1)/\mu)\exp(-(x+1)/\mu))}{(1-(1-\sigma)(1+x/\mu)\exp(-x/\mu)) ((1-(1-\sigma)(1+(x+1)/\mu)\exp(-(x+1)/\mu))}}
#'
#' with \eqn{\mu > 0} and \eqn{\sigma > 0}
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a DMOLBE distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_DMOLBE.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DMOLBE <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "DMOLBE",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "DMOLBE",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("DMOLBE", "Discrete Marshall-Olkin Length Biased Exponential"),
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
                   t1 <- (1 + y/mu) * exp(-y/mu)
                   t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
                   A  <- t1 - t2
                   B1 <- 1 - (1 - sigma) * t1
                   B2 <- 1 - (1 - sigma) * t2
                   dA  <- (y^2 / mu^3) * exp(-y/mu) - ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
                   dB1 <- -(1 - sigma) * (y^2 / mu^3) * exp(-y/mu)
                   dB2 <- -(1 - sigma) * ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
                   dldm <- (dA / A) - (dB1 / B1) - (dB2 / B2)

                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   t1 <- (1 + y/mu) * exp(-y/mu)
                   t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
                   B1 <- 1 - (1 - sigma) * t1
                   B2 <- 1 - (1 - sigma) * t2
                   dldd <- (1 / sigma) - (t1 / B1) - (t2 / B2)

                   dldd
                 },

                 # Second derivatives

                 d2ldm2 = function(y, mu, sigma) {
                   t1 <- (1 + y/mu) * exp(-y/mu)
                   t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
                   A  <- t1 - t2
                   B1 <- 1 - (1 - sigma) * t1
                   B2 <- 1 - (1 - sigma) * t2
                   dA  <- (y^2 / mu^3) * exp(-y/mu) - ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
                   dB1 <- -(1 - sigma) * (y^2 / mu^3) * exp(-y/mu)
                   dB2 <- -(1 - sigma) * ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
                   dldm <- (dA / A) - (dB1 / B1) - (dB2 / B2)

                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   t1 <- (1 + y/mu) * exp(-y/mu)
                   t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
                   A  <- t1 - t2
                   B1 <- 1 - (1 - sigma) * t1
                   B2 <- 1 - (1 - sigma) * t2
                   dA  <- (y^2 / mu^3) * exp(-y/mu) - ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
                   dB1 <- -(1 - sigma) * (y^2 / mu^3) * exp(-y/mu)
                   dB2 <- -(1 - sigma) * ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
                   dldm <- (dA / A) - (dB1 / B1) - (dB2 / B2)

                   t1 <- (1 + y/mu) * exp(-y/mu)
                   t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
                   B1 <- 1 - (1 - sigma) * t1
                   B2 <- 1 - (1 - sigma) * t2
                   dldd <- (1 / sigma) - (t1 / B1) - (t2 / B2)

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   t1 <- (1 + y/mu) * exp(-y/mu)
                   t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
                   B1 <- 1 - (1 - sigma) * t1
                   B2 <- 1 - (1 - sigma) * t2
                   dldd <- (1 / sigma) - (t1 / B1) - (t2 / B2)

                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dDMOLBE(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pDMOLBE", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_DMOLBE(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_DMOLBE(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu, sigma, x_max=100) {
                   x <- 1:x_max
                   num <- (1+x/mu) * exp(-x/mu)
                   den <- 1 - (1-sigma) *(1+x/mu) * exp(-x/mu)
                   the_mean <- sigma * sum(num/den)
                   return(the_mean)
                 },

                 variance = function(mu, sigma, x_max=100) {
                   x <- 1:x_max
                   num <- (1+x/mu) * exp(-x/mu)
                   den <- 1 - (1-sigma) *(1+x/mu) * exp(-x/mu)
                   the_mean <- sigma * sum(num/den)
                   the_var <- sigma * sum((2*x+1) * num/den) - the_mean^2
                   return(the_var)
                 }

  ),
  class=c("gamlss.family", "family"))
}
#' logLik function for DMOLBE
#' @description Calculates logLik for DMOLBE distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @export
logLik_DMOLBE <- function(logparam=c(0, 0), x){
  return(sum(dDMOLBE(x     = x,
                     mu    = exp(logparam[1]),
                     sigma = exp(logparam[2]),
                     log=TRUE)))
}
#' Initial values for DMOLBE
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_DMOLBE <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_DMOLBE,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
