# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dir Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Dir",
         contains = "Distribution",
         slots = c(alpha = "numeric"),
         prototype = list(alpha = c(1, 1)))

#' @title Dirichlet Distribution
#' @name Dir
#'
#' @description
#' The Dirichlet distribution is an absolute continuous probability,
#' specifically a multivariate generalization of the beta distribution,
#' parameterized by a vector \eqn{\boldsymbol{\alpha} =
#' (\alpha_1, \alpha_2, ..., \alpha_k)} with \eqn{\alpha_i > 0}.
#'
#' @srrstats {G1.0} A list of publications is provided in references.
#' @srrstats {G1.1} This is the first implementation the Dirichlet SAME in
#' general. The MLE algorithm is considerably improved.
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a, G2.2} Assertions on the length and type
#' of input is implemented in the `ddir` and `rdir` functions.
#' @srrstats {G2.4, G2.4a} Explicit conversion to the appropriate data type is
#' implemented in the `rdir` function.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Dir`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Dir`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param alpha numeric. The non-negative distribution parameter vector.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param log logical. Should the logarithm of the probability be
#' returned?
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization.
#'
#' @details
#' The probability density function (PDF) of the Dirichlet distribution is given
#' by:
#' \deqn{ f(x_1, ..., x_k; \alpha_1, ..., \alpha_k) =
#' \frac{1}{B(\boldsymbol{\alpha})} \prod_{i=1}^k x_i^{\alpha_i - 1}, }
#' where \eqn{B(\boldsymbol{\alpha})} is the multivariate Beta function:
#' \deqn{ B(\boldsymbol{\alpha}) = \frac{\prod_{i=1}^k
#' \Gamma(\alpha_i)}{\Gamma\left(\sum_{i=1}^k \alpha_i\right)} }
#' and \eqn{\sum_{i=1}^k x_i = 1}, \eqn{x_i > 0}.
#'
#' @inherit distributions return
#'
#' @references
#'
#' - Oikonomidis, I. & Trevezas, S. (2025), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Dir Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- c(0.5, 2, 5)
#' D <- Dir(a)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, c(0.3, 0.2, 0.5)) # density function
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
#' lldir(x, a)
#'
#' edir(x, type = "mle")
#' edir(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("dir", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vdir(a, type = "mle")
#' vdir(a, type = "me")
#'
#' avar_mle(D)
#' avar_me(D)
#'
#' v(D, type = "mle")
Dir <- function(alpha = c(1, 1)) {
  new("Dir", alpha = alpha)
}

setValidity("Dir", function(object) {
  if(length(object@alpha) <= 1) {
    stop("alpha has to be a numeric of length 2 or more")
  }
  if(any(object@alpha <= 0)) {
    stop("alpha has to be positive")
  }
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
ddir <- function(x, alpha, log = FALSE) {

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (ncol(x) != length(alpha)) {
    stop("the columns of x must be equal to the length of alpha")
  }
  if (any(alpha <= 0)) {
    stop("alpha must be positive")
  }

  y <- apply(x, 1,
        FUN = function(x) {
          if (any(x <= 0) || abs(sum(x) - 1) > 1e-10) {
            -Inf
          } else {
            lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
          }
        })

  if (!log) {
    return(exp(y))
  } else {
    return(y)
  }

}

#' @rdname Dir
#' @export
rdir <- function(n, alpha) {

  if (length(n) > 1) {
    warning("n has length > 1. The object's length will be used as the sample
            size")
    n <- length(n)
  } else if (!is.numeric(n) || n < 0) {
    stop("n must be a positive numeric (which will be converted to integer)")
  }
  if (any(alpha <= 0)) {
    stop("alpha must be positive")
  }

  k <- length(alpha)
  n <- as.integer(n)
  m <- matrix(rgamma(n * k, shape = alpha), nrow = n, ncol = k, byrow = TRUE)

  m / rowSums(m)

}

#' @rdname Dir
setMethod("d", signature = c(distr = "Dir", x = "numeric"),
          function(distr, x, log = FALSE) {
            ddir(x, alpha = distr@alpha, log = log)
          })

#' @rdname Dir
setMethod("d", signature = c(distr = "Dir", x = "matrix"),
          function(distr, x) {
            ddir(x, alpha = distr@alpha)
          })

#' @rdname Dir
setMethod("r", signature = c(distr = "Dir", n = "numeric"),
          function(distr, n) {
            rdir(n, alpha = distr@alpha)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
setMethod("mean",
          signature  = c(x = "Dir"),
          definition = function(x) {

  x@alpha / sum(x@alpha)

})

#' @rdname Dir
setMethod("mode",
          signature  = c(x = "Dir"),
          definition = function(x) {

  (x@alpha - 1) / (sum(x@alpha) - length(x@alpha))

})

#' @rdname Dir
setMethod("var",
          signature  = c(x = "Dir"),
          definition = function(x) {

  # Required variables
  a <- x@alpha
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)
  Ik <- diag(k)

  y <- - Matrix(a, k, 1) %*% Matrix(a, 1, k)
  diag(y) <- a * b
  y <- y / (a0 ^ 2 * (a0 + 1))
  as.matrix(nearPD(y))

})

#' @rdname Dir
setMethod("entro",
          signature  = c(x = "Dir"),
          definition = function(x) {

  a <- x@alpha
  a0 <- sum(a)
  ba <- sum(lgamma(a)) - lgamma(a0)

  ba + (a0 - length(a)) * digamma(a0) - sum((a - 1) * digamma(a))

})

#' @rdname Dir
setMethod("finf",
          signature  = c(x = "Dir"),
          definition = function(x) {

  a <- x@alpha
  k <- length(a)

  D <- diag(trigamma(a)) - matrix(trigamma(sum(a)), k, k)

  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
lldir <- function(x, alpha) {
  ll(distr = Dir(alpha), x)
}

#' @rdname Dir
setMethod("ll",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x) {

  a <- distr@alpha
  nrow(x) * (lgamma(sum(a)) - sum(lgamma(a))) + sum(log(x) %*% diag(a - 1))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Dir"),
          definition = function(par, tx, distr) {

  a <- idigamma(digamma(par) + tx)

  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * tx)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Dir"),
          definition = function(par, tx, distr) {

  # Shape parameters (a_i) as a function of a0
  a <- idigamma(digamma(par) + tx)

  # a_i derivative wrt a0
  da <- trigamma(par) / trigamma(a)

  # lloptim derivative wrt a0 (par)
  digamma(sum(a)) * sum(da) - sum(digamma(a) * da) + sum(tx * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
edir <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me", "same"))
  distr <- Dir()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Dir
setMethod("mle",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  par0 <- check_optim(par0, method, lower, upper,
                      choices = c("me", "same"), len = 1)

  if (is.character(par0)) {
    par0 <- sum(edir(x, type = par0)$alpha)
  }

  tx  <- colMeans(log(x))

  par <- optim(par = par0,
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  list(alpha = idigamma(digamma(par) + tx))

})

#' @rdname Dir
setMethod("me",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  m  <- colMeans(x)
  m2  <- colMeans(x ^ 2)

  a0 <- (1 - sum(m2)) / (sum(m2) - sum(m ^ 2))

  list(alpha = a0 * m)

})

#' @rdname Dir
setMethod("same",
          signature  = c(distr = "Dir", x = "matrix"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  m  <- colMeans(x)
  logm  <- colMeans(log(x))
  mlogm <- colMeans(x * log(x))

  list(alpha  = (length(m) - 1) * m / sum(mlogm - m * logm))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Dir
#' @export
vdir <- function(alpha, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me", "same"))
  distr <- Dir(alpha)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Dir
setMethod("avar_mle",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  a <- distr@alpha
  k <- length(a)
  a0 <- sum(a)
  trig <- 1 / trigamma(a)
  cons <- trigamma(a0) / (1 - trigamma(a0) * sum(trig))

  D <- diag(trig) + cons * Matrix(trig, k, 1) %*% Matrix(trig, 1, k)

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

#' @rdname Dir
setMethod("avar_me",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@alpha
  a0 <- sum(a)
  k <- length(a)
  dn <- a0 ^ 2 - sum(a ^ 2)
  a1 <- a0 + 1
  a2 <- a0 + 2
  a3 <- a0 + 3
  s2 <- sum(a^2)
  s3 <- sum(a^3)

  c0 <- ((- 4 * a0 * (a0 - 1) * a1 ^ 2 * s3 +
         (2 * a0 ^ 3 + a0 ^ 2 + a0) * s2 ^ 2 +
         (2 * a0 ^ 5 + 2 * a0 ^ 4 - 6 * a0 ^ 3 - 4 * a0 ^ 2 - 2 * a0) * s2 +
         a0 ^ 6 + a0 ^ 5 + 2 * a0 ^ 3) / (dn ^ 2 * a1 * a2 * a3))

  D <- (a0 / a1) * diag(a) + 2*a0 / (dn*a2) * (a %*% t(a^2) + a^2 %*% t(a)) +
    c0 * a %*% t(a)

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})

#' @rdname Dir
setMethod("avar_same",
          signature  = c(distr = "Dir"),
          definition = function(distr) {

  # Required variables
  a <- distr@alpha
  a0 <- sum(a)
  b <- a0 - a
  k <- length(a)
  Ik <- diag(k)
  Amat <- Matrix(a, k, 1) %*% Matrix(a, 1, k)
  mat2 <- Matrix(1/a, k, 1) %*% Matrix(1, 1, k)

  par1 <- 1 / ((k - 1) * (a0 + 1))
  par2 <- 1 / ((k - 1) ^ 2 * (a0 + 1))

  c <- - par2 *  sum(a^2*trigamma(a)) + a0 * par2 * sum(a*trigamma(a)) +
    (a0+2) * par1

  D <- Amat * (c - (mat2 + Matrix::t(mat2)) * a0 * par1  +
                 a0 * diag(1/ a) / (a0+1))

  D <- as.matrix(nearPD(D))
  rownames(D) <- paste0("alpha", seq_along(a))
  colnames(D) <- paste0("alpha", seq_along(a))
  D

})
