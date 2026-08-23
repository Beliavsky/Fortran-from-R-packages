# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cat Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Cat",
         contains = "Distribution",
         slots = c(prob = "numeric"),
         prototype = list(prob = c(0.5, 0.5)))

#' @title Cat Distribution
#' @name Cat
#'
#' @description
#' The Categorical distribution is a discrete probability distribution that
#' describes the probability of a single trial resulting in one of \eqn{k}
#' possible categories. It is a generalization of the Bernoulli distribution
#' and a special case of the multinomial distribution with \eqn{n = 1}.
#'
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a, G2.2} Assertions on the length and type
#' of input is implemented.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Cat`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Cat`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param prob numeric. Probability vector of success for each category.
#' @param dim numeric. The probability vector dimension. See Details.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log logical. Should the logarithm of the probability be
#' returned?
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#'
#' @details
#' The probability mass function (PMF) of the categorical distribution is given
#' by: \deqn{ f(x; p) = \prod_{i=1}^k p_i^{x_i},}
#' subject to \eqn{ \sum_{i=1}^{k} x_i = n }.
#'
#' The estimation of `prob` from a sample would by default return a vector of
#' probabilities corresponding to the categories that appeared in the sample and
#' 0 for the rest. However, the parameter dimension cannot be uncovered by the
#' sample, it has to be provided separately. This can be done with the argument
#' `dim`. If `dim` is not supplied, the dimension will be retrieved from the
#' `distr` argument. Categories that did not appear in the sample will have 0
#' probabilities appended to the end of the prob vector.
#'
#' Note that the actual dimension of the probability parameter vector is `k-1`,
#' therefore the Fisher information matrix and the asymptotic variance -
#' covariance matrix of the estimators is of dimension `(k-1)x(k-1)`.
#'
#' @inherit distributions return
#'
#' @seealso
#' [dmultinom()], [rmultinom()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Categorical Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' p <- c(0.1, 0.2, 0.7)
#' D <- Cat(p)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, 2) # density function
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
#' mode(D) # Mode
#' var(D) # Variance
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
#' llcat(x, p)
#'
#' ecat(x, dim = 3, type = "mle")
#' ecat(x, dim = 3, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("cat", dim = 3, x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vcat(p, type = "mle")
#' vcat(p, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Cat <- function(prob = c(0.5, 0.5)) {
  new("Cat", prob = prob)
}

setValidity("Cat", function(object) {
  if(length(object@prob) <= 1) {
    stop("prob has to be a numeric of length at least 2")
  }
  if(any(object@prob <= 0) || any(object@prob >= 1)) {
    stop("prob has to be between 0 and 1")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
#' @export
dcat <- function(x, prob, log = FALSE) {

  if (any(prob < 0) || abs(sum(prob) - 1) > 1e-8) {
    stop("prob must be a valid probability vector")
  }

  if (any(x %% 1 != 0)) {
    warning("non-integer x")
  }

  y <- unlist(lapply(x, function(x) {
                if (x %in% seq_along(prob)) {
                  return(prob[x])
                } else {
                  return(0)
                }
              }))

  if (log) {
    return(log(y))
  } else {
    return(y)
  }

}

#' @rdname Cat
#' @export
rcat <- function(n, prob) {
  if (length(n) > 1) {
    warning("n has length > 1. The object's length will be used as the sample
            size")
    n <- length(n)
  } else if (!is.numeric(n) || n < 0) {
    stop("n must be a positive numeric (which will be converted to integer)")
  }
  n <- as.integer(n)
  sample(seq_along(prob), n, prob = prob, replace = TRUE)
}

#' @rdname Cat
setMethod("d", signature = c(distr = "Cat", x = "numeric"),
          function(distr, x, log = FALSE) {
            dcat(x, prob = distr@prob, log = log)
          })

#' @rdname Cat
setMethod("r", signature = c(distr = "Cat", n = "numeric"),
          function(distr, n) {
            rcat(n, prob = distr@prob)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
setMethod("mean",
          signature  = c(x = "Cat"),
          definition = function(x) {

  x@prob

})

#' @rdname Cat
setMethod("mode",
          signature  = c(x = "Cat"),
          definition = function(x) {

  which(x@prob == max(x@prob))

})

#' @rdname Cat
setMethod("var",
          signature  = c(x = "Cat"),
          definition = function(x) {

  k <- length(x@prob)

  diag(x@prob) - matrix(x@prob, k, 1) %*% matrix(x@prob, 1, k)

})

#' @rdname Cat
setMethod("entro",
          signature  = c(x = "Cat"),
          definition = function(x) {

  p <- x@prob

  - p * log(p) - (1 - p) * log(1 - p)

})

#' @rdname Cat
setMethod("finf",
          signature  = c(x = "Cat"),
          definition = function(x) {

  k <- length(x@prob)

  if (k == 2) {
    y <- 1 / x@prob[-k]
  } else {
    y <- diag(1 / x@prob[-k])
  }
  D <- y + matrix(1, k - 1, 1) %*% matrix(1, 1, k - 1) /
    x@prob[k]

  rownames(D) <- paste0("prob", seq_along(x@prob[-k]))
  colnames(D) <- paste0("prob", seq_along(x@prob[-k]))
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
#' @export
llcat <- function(x, prob) {
  ll(Cat(prob), x)
}

#' @rdname Cat
setMethod("ll",
          signature  = c(distr = "Cat", x = "numeric"),
          definition = function(distr, x) {

  sum(log(distr@prob[x]))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
#' @export
ecat <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Cat()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Cat
setMethod("mle",
          signature  = c(distr = "Cat", x = "numeric"),
          definition = function(distr, x, dim = NULL, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  if (is.null(dim)) {
    dim <- length(distr@prob)
  }

  p <- unname(table(x) / length(x))

  if (dim < length(p)) {
    stop("Dimension of Cat distribution supplied was ", dim, ", but ",
         length(p), " categories found in the sample.")
  }

  p <- c(p, rep(0, length = dim - length(p)))

  list(prob = p)

})

#' @rdname Cat
setMethod("me",
          signature  = c(distr = "Cat", x = "numeric"),
          definition = function(distr, x, dim = NULL, na.rm = FALSE) {

  mle(distr, x, dim, na.rm = na.rm)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cat
#' @export
vcat <- function(prob, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Cat(prob)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Cat
setMethod("avar_mle",
          signature  = c(distr = "Cat"),
          definition = function(distr) {

  as.matrix(nearPD(solve(finf(distr))))

})

#' @rdname Cat
setMethod("avar_me",
          signature  = c(distr = "Cat"),
          definition = function(distr) {

  avar_mle(distr)

})
