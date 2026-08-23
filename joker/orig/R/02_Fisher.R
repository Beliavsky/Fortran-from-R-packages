# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fisher Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Fisher",
         contains = "Distribution",
         slots = c(df1 = "numeric", df2 = "numeric"),
         prototype = list(df1 = 1, df2 = 1))

#' @title Fisher Distribution
#' @name Fisher
#'
#' @description
#' The Fisher (F) distribution is an absolute continuous probability
#' distribution that arises frequently in the analysis of variance (ANOVA) and
#' in hypothesis testing. It is defined by two degrees of freedom parameters
#' \eqn{d_1 > 0} and \eqn{d_2 > 0}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Fisher`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Fisher`. For the
#' log-likelihood functions, `x` is the sample of observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param df1,df2 numeric. The distribution degrees of freedom parameters.
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#'
#' @details
#' The probability density function (PDF) of the F-distribution is given by:
#' \deqn{ f(x; d_1, d_2) = \frac{\sqrt{\left(\frac{d_1 x}{d_1 x +
#' d_2}\right)^{d_1} \left(\frac{d_2}{d_1 x + d_2}\right)^{d_2}}}{x B(d_1/2,
#' d_2/2)}, \quad x > 0 .}
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [df()], [pf()], [qf()], [rf()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Fisher Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' df1 <- 14 ; df2 <- 20
#' D <- Fisher(df1, df2)
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
#' llf(x, df1, df2)
#'
Fisher <- function(df1 = 1, df2 = 1) {
  new("Fisher", df1 = df1, df2 = df2)
}

setValidity("Fisher", function(object) {
  if(length(object@df1) != 1) {
    stop("df1 has to be a numeric of length 1")
  }
  if(length(object@df2) != 1) {
    stop("df2 has to be a numeric of length 1")
  }
  if(object@df1 <= 0) {
    stop("df1 has to be positive")
  }
  if(object@df2 <= 0) {
    stop("df2 has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Fisher
setMethod("d", signature = c(distr = "Fisher", x = "numeric"),
          function(distr, x, log = FALSE) {
            df(x, df1 = distr@df1, df2 = distr@df2, ncp = 0, log = log)
          })

#' @rdname Fisher
setMethod("p", signature = c(distr = "Fisher", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pf(q, df1 = distr@df1, df2 = distr@df2, ncp = 0,
               lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Fisher
setMethod("qn", signature = c(distr = "Fisher", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qf(p, df1 = distr@df1, df2 = distr@df2, ncp = 0,
               lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Fisher
setMethod("r", signature = c(distr = "Fisher", n = "numeric"),
          function(distr, n) {
            rf(n, df1 = distr@df1, df2 = distr@df2, ncp = 0)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Fisher
setMethod("mean",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  if (x@df2 > 2) {
    return(x@df2 * x@df1 / (x@df1 * (x@df2 - 2)))
  } else {
    stop("Expectation is undefined for F distribution with df2 <= 2.")
  }

})

#' @rdname Fisher
setMethod("median",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  qf(0.5, df1 = x@df1, df2 = x@df2)

})

#' @rdname Fisher
setMethod("mode",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  if (x@df1 > 2) {
    return(x@df1 - 2) * x@df2 / (x@df1 * (x@df2 + 2))
  } else {
    stop("Expectation is undefined for F distribution with df1 <= 2.")
  }

})

#' @rdname Fisher
setMethod("var",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (x@df2 > 4) {
    return(2 * (n1 ^ 2 + n1 * (n2 - 2)) *
             (n2 / n1) ^ 2 / ((n2 - 2) ^ 2 * (n2 - 4)))
  } else {
    stop("Variance is undefined for Fdistribution with df2 <= 4.")
  }

})

#' @rdname Fisher
setMethod("sd",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Fisher
setMethod("skew",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (n2 > 6) {
    ((2 * n1 + n2 - 2) * sqrt(8 * (n2 - 4))) /
      ((n2 - 6) * sqrt(n1 * (n1 + n2 - 2)))
  }  else {
    stop("Skewness is undefined for F distribution with df2 <= 6.")
  }

})

#' @rdname Fisher
setMethod("kurt",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  if (n2 > 8) {
    12 * (n1 * (5 * n2 - 22) * (n1 + n2 -2) + (n2 - 4) * (n2 - 2) ^ 2) /
      n1 * (n2 - 6) * (n2 - 8) * (n1 + n2 - 2)
  }  else {
    stop("Kurtosis is undefined for F distribution with df2 <= 8.")
  }

})

#' @rdname Fisher
setMethod("entro",
          signature  = c(x = "Fisher"),
          definition = function(x) {

  n1 <- x@df1
  n2 <- x@df2

  lgamma(n1 / 2) + lgamma(n2 / 2) - lgamma((n1 + n2) / 2) +
    (1 - n1 / 2) * digamma(1 + n1 / 2) - (1 + n2 / 2) * digamma(1 + n2 / 2) +
    ((n1 + n2) / 2) * digamma((n1 + n2) / 2) + log(n2) - log(n1)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Fisher
#' @export
llf <- function(x, df1, df2) {
  ll(Fisher(df1, df2), x)
}

#' @rdname Fisher
setMethod("ll",
          signature  = c(distr = "Fisher", x = "numeric"),
          definition = function(distr, x) {

  d1 <- distr@df1
  d2 <- distr@df2
  n <- length(x)
  s <- sum(log(x))
  t <- sum(log(d1 * x + d2))

  (n * d1 * log(d1) + n * d2 * log(d2) + d1 * s - (d1 + d2) * t) / 2 -
    s - n * lbeta(d1 / 2, d2 / 2)

})
