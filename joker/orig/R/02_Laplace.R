# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Laplace Distribution                                                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Laplace",
         contains = "Distribution",
         slots = c(mu = "numeric", sigma = "numeric"),
         prototype = list(mu = 0, sigma = 1))

#' @title Laplace Distribution
#' @name Laplace
#'
#' @description
#' The Laplace distribution, also known as the double exponential distribution,
#' is a continuous probability distribution that is often used to model data
#' with sharp peaks and heavy tails. It is parameterized by a location parameter
#' \eqn{\mu} and a scale parameter \eqn{b > 0}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Laplace`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Laplace`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param mu,sigma numeric. The distribution parameters.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability density function (PDF) of the Laplace distribution is:
#' \deqn{ f(x; \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right) .}
#'
#' @inherit distributions return
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Laplace Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' m <- 3 ; s <- 5
#' D <- Laplace(m, s)
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
#' elaplace(x, type = "mle")
#' elaplace(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("laplace", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vlaplace(m, s, type = "mle")
#' vlaplace(m, s, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Laplace <- function(mu = 0, sigma = 1) {
  new("Laplace", mu = mu, sigma = sigma)
}

setValidity("Laplace", function(object) {
  if(length(object@mu) != 1) {
    stop("mu has to be a numeric of length 1")
  }
  if(length(object@mu) != 1) {
    stop("mu has to be a numeric of length 1")
  }
  if(object@sigma <= 0) {
    stop("sigma has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
dlaplace <- function(x, mu, sigma, log = FALSE) {
  y <- 0.5 * dexp(abs(x - mu), rate = 1 / sigma)
  if (log) {
    return(log(y))
  } else {
    return(y)
  }
}

#' @rdname Laplace
#' @export
plaplace <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  y <- unlist(lapply(q, function(x) {
    if (x >= mu) {
      return(1 - 0.5 * exp((mu - x) / sigma))
    } else {
      return(0.5 * exp((x - mu) / sigma))
    }
  }))

  if (!lower.tail) {
    y <- 1 - y
  }
  if (log.p) {
    return(log(y))
  } else {
    return(y)
  }

}

#' @rdname Laplace
#' @export
qlaplace <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE) {

  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }

  unlist(lapply(p, function(x) {
    if (x >= 0.5) {
      return(mu + qexp(2 * (x - 0.5), rate = 1 / sigma))
    } else {
      return(mu - qexp(2 * x, rate = 1 / sigma, lower.tail = FALSE))
    }
  }))

}

#' @rdname Laplace
#' @export
rlaplace <- function(n, mu, sigma) {
  (2 * rbern(n, 0.5) - 1) * rexp(n, rate = 1 / sigma) + mu
}

#' @rdname Laplace
setMethod("d", signature = c(distr = "Laplace", x = "numeric"),
          function(distr, x, log = FALSE) {
            dlaplace(x, distr@mu, distr@sigma, log = log)
          })

#' @rdname Laplace
setMethod("p", signature = c(distr = "Laplace", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            plaplace(q, distr@mu, distr@sigma,
                     lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Laplace
setMethod("qn", signature = c(distr = "Laplace", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qlaplace(p, distr@mu, distr@sigma,
                     lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Laplace
setMethod("r", signature = c(distr = "Laplace", n = "numeric"),
          function(distr, n) {
            rlaplace(n, distr@mu, distr@sigma)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
setMethod("mean",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  x@mu

})

#' @rdname Laplace
setMethod("median",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  x@mu

})

#' @rdname Laplace
setMethod("mode",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  x@mu

})

#' @rdname Laplace
setMethod("var",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  2 * x@sigma ^ 2

})

#' @rdname Laplace
setMethod("sd",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Laplace
setMethod("skew",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  0

})

#' @rdname Laplace
setMethod("kurt",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  6

})

#' @rdname Laplace
setMethod("entro",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  log(2 * x@sigma * exp(1))

})

#' @rdname Laplace
setMethod("finf",
          signature  = c(x = "Laplace"),
          definition = function(x) {

  mat <- matrix(c(1, 0, 0, 1 / x@sigma), 2, 2)
  prm_names <- c("mu", "sigma")
  dimnames(mat) <- list(prm_names, prm_names)

  mat

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
lllaplace <- function(x, mu, sigma) {
  ll(distr = Laplace(mu, sigma), x)
}

#' @rdname Laplace
setMethod("ll",
          signature  = c(distr = "Laplace", x = "numeric"),
          definition = function(distr, x) {

  - length(x) * log(2 * distr@sigma) - sum(abs(x - distr@mu)) / distr@sigma

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
elaplace <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Laplace()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Laplace
setMethod("mle",
          signature  = c(distr = "Laplace", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  m <- median(x)
  list(mu = m, sigma = mean(abs(x - m)))

})

#' @rdname Laplace
setMethod("me",
          signature  = c(distr = "Laplace", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  mle(distr, x, na.rm = na.rm)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Laplace
#' @export
vlaplace <- function(mu, sigma, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Laplace(mu, sigma)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Laplace
setMethod("avar_mle",
          signature  = c(distr = "Laplace"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Laplace
setMethod("avar_me",
          signature  = c(distr = "Laplace"),
          definition = function(distr) {

  avar_mle(distr)

})
