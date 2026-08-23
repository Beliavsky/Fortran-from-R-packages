# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Binom Distribution                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Binom",
         contains = "Distribution",
         slots = c(size = "numeric", prob = "numeric"),
         prototype = list(size = 1, prob = 0.5))

#' @title Binom Distribution
#' @name Binom
#'
#' @description
#' The binomial distribution is a discrete probability distribution which models
#' the probability of having x successes in n independent Bernoulli trials with
#' success probability p.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Binom`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Binom`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param size number of trials (zero or more).
#' @param prob numeric. Probability of success on each trial.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability mass function (PMF) of the binomial distribution is given
#' by: \deqn{ f(x; n, p) = \binom{n}{x} p^x (1 - p)^{n - x}, \quad N \in
#' \mathbb{N}, \quad p \in (0, 1),} with \eqn{x \in \{0, 1, \dots, N\}}.
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dbinom()], [pbinom()], [qbinom()],
#' [rbinom()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Binomial Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' N <- 10 ; p <- 0.7
#' D <- Binom(N, p)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, 0:N) # density function
#' p(D, 0:N) # distribution function
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
#' llbinom(x, N, p)
#'
#' ebinom(x, size = N, type = "mle")
#' ebinom(x, size = N, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vbinom(N, p, type = "mle")
#' vbinom(N, p, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Binom <- function(size = 1, prob = 0.5) {
  new("Binom", size = size, prob = prob)
}

setValidity("Binom", function(object) {
  if(length(object@size) != 1) {
    stop("size has to be a numeric of length 1")
  }
  if(!is_natural(object@size)) {
    stop("size has to be a natural number")
  }
  if(length(object@prob) != 1) {
    stop("prob has to be a numeric of length 1")
  }
  if(object@prob <= 0 || object@prob >= 1) {
    stop("prob has to be between 0 and 1")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
setMethod("d", signature = c(distr = "Binom", x = "numeric"),
          function(distr, x, log = FALSE) {
            dbinom(x, size = distr@size, prob = distr@prob, log = log)
          })

#' @rdname Binom
setMethod("p", signature = c(distr = "Binom", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pbinom(q, size = distr@size, prob = distr@prob,
                   lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Binom
setMethod("qn", signature = c(distr = "Binom", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qbinom(p, size = distr@size, prob = distr@prob,
                   lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Binom
setMethod("r", signature = c(distr = "Binom", n = "numeric"),
          function(distr, n) {
            rbinom(n, size = distr@size, prob = distr@prob)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
setMethod("mean",
          signature  = c(x = "Binom"),
          definition = function(x) {

  x@size * x@prob

})

#' @rdname Binom
setMethod("var",
          signature  = c(x = "Binom"),
          definition = function(x) {

  x@size * x@prob * (1 - x@prob)

})

#' @rdname Binom
setMethod("sd",
          signature  = c(x = "Binom"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Binom
setMethod("skew",
          signature  = c(x = "Binom"),
          definition = function(x) {

  p <- x@prob
  q <- 1 - p
  (q - p) / sqrt(x@size * p * q)

})

#' @rdname Binom
setMethod("kurt",
          signature  = c(x = "Binom"),
          definition = function(x) {

  p <- x@prob
  q <- 1 - p

  (1 - 6 * p * q) / (x@size * p * q)

})

#' @rdname Binom
setMethod("entro",
          signature  = c(x = "Binom"),
          definition = function(x) {

  warning("The entropy given is an approximation in the O(1 / n) order.")
  p <- x@prob
  0.5 * log(2 * pi  * exp(1) * x@size * p * (1 - p), base = 2)

})

#' @rdname Binom
setMethod("finf",
          signature  = c(x = "Binom"),
          definition = function(x) {

  x@size / (x@prob * (1 - x@prob))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
#' @export
llbinom <- function(x, size, prob) {
  ll(distr = Binom(size, prob), x)
}

#' @rdname Binom
setMethod("ll",
          signature  = c(distr = "Binom", x = "numeric"),
          definition = function(distr, x) {

  N <- distr@size
  p <- distr@prob
  n <- length(x)
  s <- sum(x)
  y <- sum(unlist(lapply(x, FUN = function(x) { lchoose(N, x) })))

  log(p) * s + log(1 - p) * (n * N - s) + y

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
#' @export
ebinom <- function(x, size, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Binom(size)
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Binom
setMethod("mle",
          signature  = c(distr = "Binom", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  p <- mean(x) / distr@size

  if (p > 1) {
    stop("Success probability ", p, ", greater than 1.
          Did you forget to specify the size of the Binomial?")
  }

  list(prob = p)

})

#' @rdname Binom
setMethod("me",
          signature  = c(distr = "Binom", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  mle(distr, x, na.rm = na.rm)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Binom
#' @export
vbinom <- function(size, prob, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Binom(size, prob)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Binom
setMethod("avar_mle",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  prob <- distr@prob
  size <- distr@size
  c(prob = prob * (1 - prob) / size)

})

#' @rdname Binom
setMethod("avar_me",
          signature  = c(distr = "Binom"),
          definition = function(distr) {

  avar_mle(distr)

})
