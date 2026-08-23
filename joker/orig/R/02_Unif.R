# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unif Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Unif",
         contains = "Distribution",
         slots = c(min = "numeric", max = "numeric"),
         prototype = list(min = 0, max = 1))

#' @title Uniform Distribution
#' @name Unif
#'
#' @description
#' The Uniform distribution is an absolute continuous probability distribution
#' where all intervals of the same length within the distribution's support are
#' equally probable. It is defined by two parameters: the lower bound \eqn{a}
#' and the upper bound \eqn{b}, with \eqn{a < b}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Unif`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Unif`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param min,max numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Uniform distribution is:
#' \deqn{ f(x; a, b) = \frac{1}{b - a}, \quad a \le x \le b .}
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dunif()], [punif()], [qunif()],
#' [runif()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Uniform Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- 3 ; b <- 5
#' D <- Unif(a, b)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, c(0.3, 0.8, 0.5)) # density function
#' p(D, c(0.3, 0.8, 0.5)) # distribution function
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
#' var(D) # Variance
#' sd(D) # Standard Deviation
#' skew(D) # Skewness
#' kurt(D) # Excess Kurtosis
#' entro(D) # Entropy
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
#' llunif(x, a, b)
#'
#' eunif(x, type = "mle")
#' eunif(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("unif", x) # the distr argument can be a character
Unif <- function(min = 0, max = 1) {
  new("Unif", min = min, max = max)
}

setValidity("Unif", function(object) {
  if(length(object@min) != 1) {
    stop("min has to be a numeric of length 1")
  }
  if(length(object@max) != 1) {
    stop("max has to be a numeric of length 1")
  }
  if(object@min >= object@max) {
    stop("min must be less than max")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
setMethod("d", signature = c(distr = "Unif", x = "numeric"),
          function(distr, x, log = FALSE) {
            dunif(x, min = distr@min, max = distr@max, log = log)
          })

#' @rdname Unif
setMethod("p", signature = c(distr = "Unif", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            punif(q, min = distr@min, max = distr@max,
                  lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Unif
setMethod("qn", signature = c(distr = "Unif", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qunif(p, min = distr@min, max = distr@max,
                  lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Unif
setMethod("r", signature = c(distr = "Unif", n = "numeric"),
          function(distr, n) {
            runif(n, min = distr@min, max = distr@max)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
setMethod("mean",
          signature  = c(x = "Unif"),
          definition = function(x) {

  (x@max + x@min) / 2

})


#' @rdname Unif
setMethod("median",
          signature  = c(x = "Unif"),
          definition = function(x) {

            (x@max + x@min) / 2

          })

#' @rdname Unif
setMethod("mode",
          signature  = c(x = "Unif"),
          definition = function(x) {

            warning("The mode is any element in the support (or its interior) of
            a Uniform distribution. The mean is returned by default.")
            return((x@max + x@min) / 2)

          })

#' @rdname Unif
setMethod("var",
          signature  = c(x = "Unif"),
          definition = function(x) {

  (x@max - x@min) ^ 2 / 12

})

#' @rdname Unif
setMethod("sd",
          signature  = c(x = "Unif"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Unif
setMethod("skew",
          signature  = c(x = "Unif"),
          definition = function(x) {

  0

})

#' @rdname Unif
setMethod("kurt",
          signature  = c(x = "Unif"),
          definition = function(x) {

  - 1.2

})

#' @rdname Unif
setMethod("entro",
          signature  = c(x = "Unif"),
          definition = function(x) {

  log(x@max - x@min)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
#' @export
llunif <- function(x, min, max) {
  ll(Unif(min, max), x)
}

#' @rdname Unif
setMethod("ll",
          signature  = c(distr = "Unif", x = "numeric"),
          definition = function(distr, x) {

  m <- distr@min
  M <- distr@max
  if (max(x) > M || min(x) < m) {
    return(0)
  } else {
    return(- length(x) * log(M - m))
  }

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Unif
#' @export
eunif <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Unif()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Unif
setMethod("mle",
          signature  = c(distr = "Unif", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  list(min = min(x), max = max(x))

})

#' @rdname Unif
setMethod("me",
          signature  = c(distr = "Unif", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  m <- mean(x)
  s <- sqrt(3) * bsd(x)

  list(min = m - s, max = m + s)

})
