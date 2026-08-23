# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Norm Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Norm",
         contains = "Distribution",
         slots = c(mean = "numeric", sd = "numeric"),
         prototype = list(mean = 0, sd = 1))

#' @title Normal Distribution
#' @name Norm
#'
#' @description
#' The Normal or Gaussian distribution, is an absolute continuous probability
#' distribution characterized by two parameters: the mean \eqn{\mu} and the
#' standard deviation \eqn{\sigma > 0}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Norm`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Norm`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param mean,sd numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Normal distribution is:
#' \deqn{ f(x; \mu, \sigma) = \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{1}{2}
#' \left(\frac{x - \mu}{\sigma}\right)^2} .}
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dnorm()], [pnorm()], [qnorm()],
#' [rnorm()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Normal Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' m <- 3 ; s <- 5
#' D <- Norm(m, s)
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
#' enorm(x, type = "mle")
#' enorm(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("norm", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vnorm(m, s, type = "mle")
#' vnorm(m, s, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Norm <- function(mean = 0, sd = 1) {
  new("Norm", mean = mean, sd = sd)
}

setValidity("Norm", function(object) {
  if(length(object@mean) != 1) {
    stop("mean has to be a numeric of length 1")
  }
  if(length(object@sd) != 1) {
    stop("sd has to be a numeric of length 1")
  }
  if(object@sd <= 0) {
    stop("sd has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
setMethod("d", signature = c(distr = "Norm", x = "numeric"),
          function(distr, x, log = FALSE) {
            dnorm(x, mean = distr@mean, sd = distr@sd, log = log)
          })

#' @rdname Norm
setMethod("p", signature = c(distr = "Norm", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pnorm(q, mean = distr@mean, sd = distr@sd,
                  lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Norm
setMethod("qn", signature = c(distr = "Norm", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qnorm(p, mean = distr@mean, sd = distr@sd,
                  lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Norm
setMethod("r", signature = c(distr = "Norm", n = "numeric"),
          function(distr, n) {
            rnorm(n, mean = distr@mean, sd = distr@sd)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
setMethod("mean",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("median",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("mode",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@mean

})

#' @rdname Norm
setMethod("var",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@sd ^ 2

})

#' @rdname Norm
setMethod("sd",
          signature  = c(x = "Norm"),
          definition = function(x) {

  x@sd

})

#' @rdname Norm
setMethod("skew",
          signature  = c(x = "Norm"),
          definition = function(x) {

  0

})

#' @rdname Norm
setMethod("kurt",
          signature  = c(x = "Norm"),
          definition = function(x) {

  0

})

#' @rdname Norm
setMethod("entro",
          signature  = c(x = "Norm"),
          definition = function(x) {

  log(2 * pi * exp(1) * x@sd ^ 2) / 2

})

#' @rdname Norm
setMethod("finf",
          signature  = c(x = "Norm"),
          definition = function(x) {

  sd <- x@sd

  mat <- matrix(c(1 / sd ^ 2, 0, 0, 2 / sd ^ 2), 2, 2)
  prm_names <- c("mean", "sd")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
#' @export
llnorm <- function(x, mean, sd) {
  ll(Norm(mean, sd), x)
}

#' @rdname Norm
setMethod("ll",
          signature  = c(distr = "Norm", x = "numeric"),
          definition = function(distr, x) {

  m <- distr@mean
  s <- distr@sd

  - 0.5 * length(x) * log(2 * pi * s ^ 2) -
              0.5 * sum((x - m) ^ 2) / s ^ 2

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
#' @export
enorm <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Norm()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Norm
setMethod("mle",
          signature  = c(distr = "Norm", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  list(mean = mean(x), sd = bsd(x))

})

#' @rdname Norm
setMethod("me",
          signature  = c(distr = "Norm", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  mle(distr, x, na.rm = na.rm)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Norm
#' @export
vnorm <- function(mean, sd, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Norm(mean, sd)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Norm
setMethod("avar_mle",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Norm
setMethod("avar_me",
          signature  = c(distr = "Norm"),
          definition = function(distr) {

  avar_mle(distr)

})
