# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chisq Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Chisq",
         contains = "Distribution",
         slots = c(df = "numeric"),
         prototype = list(df = 1))

#' @title Chi-Square Distribution
#' @name Chisq
#'
#' @description
#' The Chi-Square distribution is a continuous probability distribution commonly
#' used in statistical inference, particularly in hypothesis testing and
#' confidence interval estimation. It is defined by the degrees of freedom
#' parameter \eqn{k > 0}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Chisq`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Chisq`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param df numeric. The distribution degrees of freedom parameter.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Chi-Square distribution is
#' given by: \deqn{ f(x; k) = \frac{1}{2^{k/2}\Gamma(k/2)} x^{k/2 - 1} e^{-x/2},
#' \quad x > 0.}
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dchisq()], [pchisq()], [qchisq()],
#' [rchisq()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Chi-Square Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' df <- 4
#' D <- Chisq(df)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, c(0.3, 2, 20)) # density function
#' p(D, c(0.3, 2, 20)) # distribution function
#' qn(D, c(0.4, 0.8)) # inverse distribution function
#' x <- r(D, 100) # random generator function
#'
#' # alternative way to use the function
#' den <- d(D) ; den(x) # den is a function itself
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
#' llchisq(x, df)
#'
#' echisq(x, type = "mle")
#' echisq(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("chisq", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vchisq(df, type = "mle")
#' vchisq(df, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Chisq <- function(df = 1) {
  new("Chisq", df = df)
}

setValidity("Chisq", function(object) {
  if(length(object@df) != 1) {
    stop("df has to be a numeric of length 1")
  }
  if(object@df < 0) {
    stop("df has to be non-negative")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
setMethod("d", signature = c(distr = "Chisq", x = "numeric"),
          function(distr, x, log = FALSE) {
            dchisq(x, df = distr@df, ncp = 0, log = log)
          })

#' @rdname Chisq
setMethod("p", signature = c(distr = "Chisq", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pchisq(q, df = distr@df, ncp = 0,
                   lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Chisq
setMethod("qn", signature = c(distr = "Chisq", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qchisq(p, df = distr@df, ncp = 0,
                   lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Chisq
setMethod("r", signature = c(distr = "Chisq", n = "numeric"),
          function(distr, n) {
            rchisq(n, df = distr@df, ncp = 0)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
setMethod("mean",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  x@df

})

#' @rdname Chisq
setMethod("median",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  warning("The median of the Chi-Squared Distribution is not
          available in closed-form. An approximation is provided.")

  k <- x@df
  k * (1 - 2 / (9 * k)) ^ 3

})

#' @rdname Chisq
setMethod("mode",
          signature  = c(x = "Chisq"),
          definition = function(x) {

            max(x@df - 2, 0)

          })

#' @rdname Chisq
setMethod("var",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  2 * x@df

})

#' @rdname Chisq
setMethod("sd",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Chisq
setMethod("skew",
          signature  = c(x = "Chisq"),
          definition = function(x) {

 2 ^ 1.5 / sqrt(x@df)

})

#' @rdname Chisq
setMethod("kurt",
          signature  = c(x = "Chisq"),
          definition = function(x) {

12 / x@df

})

#' @rdname Chisq
setMethod("entro",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  k2 <- x@df / 2

  k2 + log(2 * gamma(k2)) + (1 - k2) * digamma(k2)

})

#' @rdname Chisq
setMethod("finf",
          signature  = c(x = "Chisq"),
          definition = function(x) {

  c(df = 0.25 * trigamma(x@df / 2))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
#' @export
llchisq <- function(x, df) {
  ll(Chisq(df), x)
}

#' @rdname Chisq
setMethod("ll",
          signature  = c(distr = "Chisq", x = "numeric"),
          definition = function(distr, x) {

    df <- distr@df
  - length(x) * (lgamma(df / 2) + df * log(2) / 2) +
    (df / 2 - 1) * sum(log(x)) - sum(x) / 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
#' @export
echisq <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Chisq()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Chisq
setMethod("mle",
          signature  = c(distr = "Chisq", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  list(df = 2 * idigamma(mean(log(x)) - log(2)))

})

#' @rdname Chisq
setMethod("me",
          signature  = c(distr = "Chisq", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  list(df = mean(x))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Chisq
#' @export
vchisq <- function(df, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Chisq(df)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Chisq
setMethod("avar_mle",
          signature  = c(distr = "Chisq"),
          definition = function(distr) {

  1 / finf(distr)

})

#' @rdname Chisq
setMethod("avar_me",
          signature  = c(distr = "Chisq"),
          definition = function(distr) {

  c(df = 2 * distr@df)

})
