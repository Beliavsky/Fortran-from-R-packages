# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Weib Distribution                                                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Weib",
         contains = "Distribution",
         slots = c(shape = "numeric", scale = "numeric"),
         prototype = list(shape = 1, scale = 1))

#' @title Weibull Distribution
#' @name Weib
#'
#' @description
#' The Weibull distribution is an absolute continuous probability distribution,
#' parameterized by a shape parameter \eqn{k > 0} and a scale parameter
#' \eqn{\lambda > 0}.
#'
#' @srrstats {G1.0} A list of publications is provided in references.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Weib`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Weib`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param shape,scale numeric. The non-negative distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me or lme).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle and me
#' optimization. See Details.
#'
#' @details
#' The probability density function (PDF) of the Weibull distribution is:
#' \deqn{ f(x; k, \lambda) = \frac{k}{\lambda}\left(\frac{x}{\lambda}
#' \right)^{k - 1} \exp\left[-\left(\frac{x}{\lambda}\right)^k\right],
#' \quad x \geq 0 .}
#'
#' For the parameter estimation, both the MLE and the ME cannot be explicitly
#' derived. However, the L-moment estimator (`type = "lme"`) is available, and
#' is used as initialization for the numerical approximation of the MLE and the
#' ME.
#'
#' The MLE and ME of the Weibull distribution parameters is not available in
#' closed form and has to be approximated numerically. The optimization can be
#' performed on the shape parameter \eqn{k\in(0,+\infty)}.
#'
#' For the MLE, this is done with `optim()`. The default method used is the
#' L-BFGS-B method with lower bound `1e-5` and upper bound `Inf`. The `par0`
#' argument can either be a numeric (satisfying `lower <= par0 <= upper`) or a
#' character specifying the closed-form estimator to be used as initialization
#' for the algorithm (`"lme"` - the default value).
#'
#' For the ME, this is done with `uniroot()`. Again, the `par0` argument can
#' either be a numeric (satisfying `lower <= par0 <= upper`) or a character
#' specifying the closed-form estimator to be used as initialization for the
#' algorithm (`"mle"` or `"lme"` - the default value). The lower and upper
#' bounds are set by default to `0.5` and `Inf`, respectively. Note that the
#' ME equations involve the \eqn{\Gamma(1 + 1 \ k)}, which can become unreliable
#' for small values of `k`, hence the `0.5` lower bound. Specifying a lower
#' bound below `0.5` will result in a warning and be ignored.
#'
#' @inherit distributions return
#'
#' @importFrom stats uniroot
#'
#' @seealso
#' Functions from the `stats` package: [dweibull()], [pweibull()], [qweibull()],
#' [rweibull()]
#'
#' @references Kim, H. M., Jang, Y. H., Arnold, B. C., & Zhao, J. (2024).
#' New efficient estimators for the Weibull distribution. Communications in
#' Statistics-Theory and Methods, 53(13), 4576-4601.
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Weibull Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- 3 ; b <- 5
#' D <- Weib(a, b)
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
#' llweibull(x, a, b)
#'
#' eweibull(x, type = "mle")
#' eweibull(x, type = "me")
#' eweibull(x, type = "lme")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("weib", x) # the distr argument can be a character
Weib <- function(shape = 1, scale = 1) {
  new("Weib", shape = shape, scale = scale)
}

setValidity("Weib", function(object) {
  if(length(object@shape) != 1) {
    stop("shape has to be a numeric of length 1")
  }
  if(object@shape <= 0) {
    stop("shape has to be positive")
  }
  if(length(object@scale) != 1) {
    stop("scale has to be a numeric of length 1")
  }
  if(object@scale <= 0) {
    stop("scale has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
setMethod("d", signature = c(distr = "Weib", x = "numeric"),
          function(distr, x, log = FALSE) {
            dweibull(x, shape = distr@shape, scale = distr@scale, log = log)
          })

#' @rdname Weib
setMethod("p", signature = c(distr = "Weib", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pweibull(q, shape = distr@shape, scale = distr@scale,
                     lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Weib
setMethod("qn", signature = c(distr = "Weib", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qweibull(p, shape = distr@shape, scale = distr@scale,
                     lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Weib
setMethod("r", signature = c(distr = "Weib", n = "numeric"),
          function(distr, n) {
            rweibull(n, shape = distr@shape, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
setMethod("mean",
          signature  = c(x = "Weib"),
          definition = function(x) {

  x@scale * gamma(1 + 1 / x@shape)

})

#' @rdname Weib
setMethod("median",
          signature  = c(x = "Weib"),
          definition = function(x) {

  x@scale * log(2) ^ (1 / x@shape)

})

#' @rdname Weib
setMethod("mode",
          signature  = c(x = "Weib"),
          definition = function(x) {

  if (x@shape > 1) {
    return(x@scale * (1 - 1 / x@shape) ^ (1 / x@shape))
  } else {
    return(0)
  }

})

#' @rdname Weib
setMethod("var",
          signature  = c(x = "Weib"),
          definition = function(x) {

  (gamma(1 + 2 / x@shape) - gamma(1 + 1 / x@shape) ^ 2) * x@scale ^ 2

})

#' @rdname Weib
setMethod("sd",
          signature  = c(x = "Weib"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Weib
setMethod("skew",
          signature  = c(x = "Weib"),
          definition = function(x) {

  m <- mean(x)
  s2 <- var(x)

  (gamma(1 + 3 / x@shape) * x@scale ^ 3 - 3 * m * s2 - m ^ 3) / s2 ^ 1.5

})

#' @rdname Weib
setMethod("kurt",
          signature  = c(x = "Weib"),
          definition = function(x) {

  g1 <- gamma(1 + 1 / x@shape)
  g2 <- gamma(1 + 2 / x@shape)
  g3 <- gamma(1 + 3 / x@shape)
  g4 <- gamma(1 + 4 / x@shape)

  (- 6 * g1 ^ 4 + 12 * g1 ^ 2 * g2 - 3 * g2 ^ 2 - 4 * g1 * g3 + g4) /
    (g2 - g1 ^ 2) ^ 2

})

#' @rdname Weib
setMethod("entro",
          signature  = c(x = "Weib"),
          definition = function(x) {

  x@scale * (1 - 1 / x@shape) + log(x@scale / x@shape) + 1

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
#' @export
llweibull <- function(x, shape, scale) {
  ll(distr = Weib(shape, scale), x)
}

#' @rdname Weib
setMethod("ll",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x) {

  k <- distr@shape
  l <- distr@scale
  n <- length(x)

  n * log(k) - n * k * log(l) + (k - 1) * sum(log(x)) - sum((x / l) ^ k)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Weib"),
          definition = function(par, tx, distr) {

  log(par) - log(mean(tx ^ par)) + (par - 1) * mean(log(tx)) - 1

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Weib"),
          definition = function(par, tx, distr) {

  par ^ (-1) + mean(log(tx)) - par * mean(tx ^ (par - 1)) / mean(tx ^ par)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Weib
#' @export
eweibull <- function(x, type = "mle", ...) {

  type <- match.arg(tolower(type), choices = c("mle", "me", "lme"))

  if (type %in% c("mle", "me")) {

    return(do.call(type, list(distr = Weib(), x = x, ...)))

  } else if (type == "lme") {

    x <- sort(x, decreasing = TRUE)
    n <- length(x)

    m1 <- mean(x)
    m2 <- 2 * sum(((n-1):0) * x) / (n * (n - 1)) - m1

    a <- - log(2) / log(1 - m2 / m1)
    b <- m1 / gamma(1 + 1 / a)

    return(list(shape = a, scale = b))

  }
}

#' @rdname Weib
setMethod("mle",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x,
                                par0 = "lme",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf,
                                na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  par0 <- check_optim(par0, method, lower, upper,
                      choices = c("lme"), len = 1)

  if (is.character(par0)) {
    par0 <- eweibull(x, type = par0)$shape
  }

  par <- optim(par = par0,
               fn = lloptim,
               gr = dlloptim,
               tx = x,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  list(shape = par, scale = mean(x ^ par) ^ (1 / par))

})

#' @rdname Weib
setMethod("me",
          signature  = c(distr = "Weib", x = "numeric"),
          definition = function(distr, x,
                                par0 = "lme",
                                lower = 0.5,
                                upper = Inf,
                                na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  if (lower < 0.5) {
    warning("lower bound provided was l < 0.5.",
            "Changed to l = 0.5 for numerical stability")
    lower <- 0.5
  }
  par0 <- check_optim(par0, method = "Brent", lower, upper,
                      choices = c("mle", "lme"), len = 1)

  if (is.character(par0)) {
    par0 <- eweibull(x, type = par0)$shape
  }

  shape <- uniroot(f = function(k) {
                    log(gamma(1 + 2 / k) - gamma(1 + 1 / k) ^ 2) -
                      2 * lgamma(1 + 1 / k) - 2 * (log(sd(x)) - log(mean(x)))
                  },
                  interval = c(max(1, par0 - 1), par0 + 1),
                  extendInt = "yes")$root

  list(shape = shape, scale = mean(x) / gamma(1 + 1 / shape))

})
