# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exp Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Exp",
         contains = "Distribution",
         slots = c(rate = "numeric"),
         prototype = list(rate = 1))

#' @title Exponential Distribution
#' @name Exp
#'
#' @description
#' The Exponential distribution is a continuous probability distribution often
#' used to model the time between independent events that occur at a constant
#' average rate. It is defined by the rate parameter \eqn{\lambda > 0}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Exp`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Exp`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param rate numeric. The distribution parameter.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Exponential distribution is
#' given by: \deqn{ f(x; \lambda) = \lambda e^{-\lambda x}, \quad x \geq 0 .}
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dexp()], [pexp()], [qexp()], [rexp()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Exp Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' rate <- 5
#' D <- Exp(rate)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, c(0.3, 2, 10)) # density function
#' p(D, c(0.3, 2, 10)) # distribution function
#' qn(D, c(0.4, 0.8)) # inverse distribution function
#' x <- r(D, 100) # random generator function
#'
#' # alternative way to use the function
#' df <- d(D) ; df(x) # df is a function itself
#'
#' # ------------------
#' # Moments
#' # ------------------
#'
#' mean(D) # Expectation
#' median(D) # Median
#' mode(D) # Mode
#' var(D) # Variance
#' sd(D) # Standard Deviation
#' skew(D) # Skewness
#' kurt(D) # Excess Kurtosis
#' entro(D) # Entropy
#' finf(D) # Fisher Information Matrix
#'
#' # List of all available moments
#' mom <- moments(D)
#' mom$mean # expectation
#'
#' # ------------------
#' # Point Estimation
#' # ------------------
#'
#' ll(D, x)
#' llexp(x, rate)
#'
#' eexp(x, type = "mle")
#' eexp(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("exp", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vexp(rate, type = "mle")
#' vexp(rate, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Exp <- function(rate = 1) {
  new("Exp", rate = rate)
}

setValidity("Exp", function(object) {
  if(length(object@rate) != 1) {
    stop("rate has to be a numeric of length 1")
  }
  if(object@rate <= 0) {
    stop("rate has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
setMethod("d", signature = c(distr = "Exp", x = "numeric"),
          function(distr, x, log = FALSE) {
            dexp(x, rate = distr@rate, log = log)
          })

#' @rdname Exp
setMethod("p", signature = c(distr = "Exp", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pexp(q, rate = distr@rate, lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Exp
setMethod("qn", signature = c(distr = "Exp", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qexp(p, rate = distr@rate, lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Exp
setMethod("r", signature = c(distr = "Exp", n = "numeric"),
          function(distr, n) {
            rexp(n, rate = distr@rate)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
setMethod("mean",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 / x@rate

})

#' @rdname Exp
setMethod("median",
          signature  = c(x = "Exp"),
          definition = function(x) {

  log(2) / x@rate

})

#' @rdname Exp
setMethod("mode",
          signature  = c(x = "Exp"),
          definition = function(x) {

  0

})

#' @rdname Exp
setMethod("var",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 / x@rate ^ 2

})

#' @rdname Exp
setMethod("sd",
          signature  = c(x = "Exp"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Exp
setMethod("skew",
          signature  = c(x = "Exp"),
          definition = function(x) {

  2

})

#' @rdname Exp
setMethod("kurt",
          signature  = c(x = "Exp"),
          definition = function(x) {

  6

})

#' @rdname Exp
setMethod("entro",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 - log(x@rate)

})

#' @rdname Exp
setMethod("finf",
          signature  = c(x = "Exp"),
          definition = function(x) {

  1 / x@rate ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
llexp <- function(x, rate) {
  ll(Exp(rate), x)
}

#' @rdname Exp
setMethod("ll",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x) {

  rate <- distr@rate
  length(x) * log(rate) - rate * sum(x)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
eexp <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Exp()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Exp
setMethod("mle",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  list(rate = 1 / mean(x))

})

#' @rdname Exp
setMethod("me",
          signature  = c(distr = "Exp", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  mle(distr, x, na.rm = na.rm)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Exp
#' @export
vexp <- function(rate, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Exp(rate)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Exp
setMethod("avar_mle",
          signature  = c(distr = "Exp"),
          definition = function(distr) {

  rate <- distr@rate
  c(rate = rate ^ 2)

})

#' @rdname Exp
setMethod("avar_me",
          signature  = c(distr = "Exp"),
          definition = function(distr) {

  avar_mle(distr)

})
