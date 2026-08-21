#' Unit Half Logistic-Geometry distribution
#'
#' @author Juan Diego Suarez Hernandez, \email{jsuarezhe@unal.edu.co}
#'
#' @description
#' The Unit Half Logistic-Geometry family
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#'
#' @seealso \link{dUHLG}
#'
#' @details
#' The Unit Half Logistic-Geometry distribution with parameter \code{mu},
#' has density given by
#'
#' \eqn{f(x| \mu) = \frac{2 \mu}{(\mu+(2-\mu)x)^2} }
#'
#' for \eqn{0 < x < 1} and \eqn{\mu > 0}.
#'
#' @returns
#' Returns a gamlss.family object which can be used to fit a
#' UHLG distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_UHLG.R
#'
#' @references
#' Ramadan, A. T., Tolba, A. H., & El-Desouky, B. S. (2022). A unit half-logistic geometric distribution and its application in insurance. Axioms, 11(12), 676.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
UHLG <- function (mu.link="log"){
  mstats <- checklink("mu.link", "Unit Half Logistic-Geometry",
                      substitute(mu.link), c("log", "own"))

  structure(list(family=c("UHLG", "Unit Half Logistic-Geometry"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Continuous",

                 mu.link    = as.character(substitute(mu.link)),

                 mu.linkfun    = mstats$linkfun,

                 mu.linkinv    = mstats$linkinv,

                 mu.dr    = mstats$mu.eta,

                 # First derivates
                 dldm = function(y, mu) {
                   dldm   <- 1/mu - 2*(1-y)/(mu+(2-mu)*y)
                   dldm
                 },

                 # Second derivates

                 d2ldm2 = function(y, mu) {
                   dm   <- gamlss::numeric.deriv(dUHLG(y, mu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 G.dev.incr = function(y, mu, ...) -2*dUHLG(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pUHLG", type="Continuous", y=y, mu=mu)),

                 mu.initial    = expression(mu    <- rep(0.5, length(y))),

                 mu.valid    = function(mu)    all(mu > 0),

                 y.valid = function(y) all(y > 0 & y < 1)


  ),
  class=c("gamlss.family", "family"))
}
