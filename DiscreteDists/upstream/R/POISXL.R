#' The Discrete Poisson XLindley
#'
#' @author Mariana Blandon Mejia, \email{mblandonm@unal.edu.co}
#'
#' @description
#' The function \code{POISXL()} defines the
#' Discrete Poisson XLindley distribution
#' a single parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#'
#' @references
#' Ahsan-ul-Haq, M., Al-Bossly, A., El-Morshedy, M., & Eliwa, M. S. (2022). Poisson XLindley distribution for count data: statistical and reliability properties with estimation techniques and inference. Computational Intelligence and neuroscience, 2022(1), 6503670.
#'
#' @seealso \link{dPOISXL}.
#'
#' @details
#' The Discrete Poisson XLindley distribution with parameters \eqn{\mu} has a support
#' 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu) = \frac{\mu^2(x+\mu^2+3(1+\mu))}{(1+\mu)^{4+x}}}; with \eqn{\mu>0}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\alpha} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a Discrete Poisson XLindley distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_POISXL.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
POISXL <- function (mu.link="log") {

  mstats <- checklink("mu.link", "POISXL",
                      substitute(mu.link), c("log"))

  structure(list(family=c("POISXL", "Poisson-XLindley"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Discrete",
                 mu.link    = as.character(substitute(mu.link)),
                 mu.linkfun = mstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 mu.dr      = mstats$mu.eta,

                 # Primeras derivadas, por ahora son computacionales
                 dldm = function(y, mu) {
                   dm   <- gamlss::numeric.deriv(dPOISXL(y, mu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 # Segundas derivadas, por ahora son computacionales
                 d2ldm2 = function(y, mu) {
                   dm   <- gamlss::numeric.deriv(dPOISXL(y, mu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 G.dev.incr = function(y, mu, pw = 1, ...) -2*dPOISXL(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pPOISXL", type="Discrete",
                                               ymin = 0, y = y, mu = mu)),

                 mu.initial = expression(mu <- rep(estim_mu_POISXL(y), length(y))),

                 mu.valid = function(mu) all(mu > 0),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu) {
                   num <- mu^2 + 2*mu + 2
                   den <- mu*(1+mu)^2
                   return(num/den)
                 },

                 variance = function(mu) {
                   numerator <- mu^5 + 5*mu^4 + 11*mu^3 + 14*mu^2 + 10*mu + 2
                   denominator <- mu^2 * (1 + mu)^4
                   return(numerator / denominator)
                 }


  ),
  class=c("gamlss.family", "family"))
}

