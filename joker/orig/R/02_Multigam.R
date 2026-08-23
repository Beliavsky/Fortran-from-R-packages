# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multigam Distribution                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Multigam",
         contains = "Distribution",
         slots = c(shape = "numeric", scale = "numeric"),
         prototype = list(shape = 1, scale = 1))

#' @title Multivariate Gamma Distribution
#' @name Multigam
#'
#' @description
#' The multivariate gamma distribution is a multivariate absolute continuous
#' probability distribution, defined as the cumulative sum of independent
#' gamma random variables with possibly different shape parameters
#' \eqn{\alpha_i > 0, i\in\{1, \dots, k\}} and the same scale \eqn{\beta > 0}.
#'
#' @srrstats {G1.0} A list of publications is provided in references.
#' @srrstats {G1.1} This is the first implementation the Multigam SAME in R. The
#' MLE algorithm is considerably improved.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Multigam`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Multigam`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param shape,scale numeric. The non-negative distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param log logical. Should the logarithm of the probability be
#' returned?
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization. See Details.
#'
#' @details
#' The probability density function (PDF) of the multivariate gamma distribution
#' is given by:
#' \deqn{ f(x; \alpha, \beta) =
#' \frac{\beta^{-\alpha_0}}{\prod_{i=1}^k\Gamma(\alpha_i)}, e^{-x_k/\beta}
#' x_1^{\alpha_1-1}\prod_{i=1}^k (x_i - x_{i-1})^{(\alpha_i-1)} \quad x > 0. }
#'
#' The MLE of the multigamma distribution parameters is not available in closed
#' form and has to be approximated numerically. This is done with `optim()`.
#' Specifically, instead of solving a \eqn{(k+1)} optimization problem w.r.t
#' \eqn{\alpha, \beta}, the optimization can be performed on the shape parameter
#' sum \eqn{\alpha_0:=\sum_{i=1}^k\alpha \in(0,+\infty)^k}. The default method
#' used is the L-BFGS-B method with lower bound `1e-5` and upper bound `Inf`.
#' The `par0` argument can either be a numeric (satisfying
#' `lower <= par0 <= upper`) or a character specifying the closed-form estimator
#' to be used as initialization for the algorithm (`"me"` or `"same"` - the
#' default value).
#'
#' @inherit distributions return
#'
#' @references
#'
#' - Mathal, A. M., & Moschopoulos, P. G. (1992). A form of multivariate gamma
#' distribution. Annals of the Institute of Statistical Mathematics, 44, 97-106.
#'
#' - Oikonomidis, I. & Trevezas, S. (2025), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Multivariate Gamma Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- c(0.5, 3, 5) ; b <- 5
#' D <- Multigam(a, b)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, c(0.3, 2, 10)) # density function
#'
#' # alternative way to use the function
#' df <- d(D) ; df(c(0.3, 2, 10)) # df is a function itself
#'
#' x <- r(D, 100) # random generator function
#'
#' # ------------------
#' # Moments
#' # ------------------
#'
#' mean(D) # Expectation
#' var(D) # Variance
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
#' llmultigam(x, a, b)
#'
#' emultigam(x, type = "mle")
#' emultigam(x, type = "me")
#' emultigam(x, type = "same")
#'
#' mle(D, x)
#' me(D, x)
#' same(D, x)
#' e(D, x, type = "mle")
#'
#' mle("multigam", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vmultigam(a, b, type = "mle")
#' vmultigam(a, b, type = "me")
#' vmultigam(a, b, type = "same")
#'
#' avar_mle(D)
#' avar_me(D)
#' avar_same(D)
#'
#' v(D, type = "mle")
Multigam <- function(shape = 1, scale = 1) {
  new("Multigam", shape = shape, scale = scale)
}

setValidity("Multigam", function(object) {
  if(any(object@shape <= 0)) {
    stop("shape has to be a vector with positive elements")
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

#' @rdname Multigam
#' @export
dmultigam <- function(x, shape, scale, log = FALSE) {

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (ncol(x) != length(shape)) {
    stop("the columns of x must be equal to the length of shape")
  }
  if (any(shape <= 0) || scale <= 0) {
    stop("shape and scale must be positive")
  }

  y <- apply(x, 1,
             FUN = function(x) {
               if (any(x <= 0)) {
                 -Inf
               } else {
                 sum(shape * log(fd(x))) - sum(shape) * log(scale) -
                   sum(lgamma(shape)) - x[length(x)] / scale
               }
             })

  if (!log) {
    return(exp(y))
  } else {
    return(y)
  }

}

#' @rdname Multigam
#' @export
rmultigam <- function(n, shape, scale) {

  k <- length(shape)
  x <- matrix(nrow = n, ncol = k)
  for (j in 1:k) {
    x[, j] <- stats::rgamma(n, shape[j], scale = scale)
  }

  t(apply(x, 1, cumsum))

}

#' @rdname Multigam
setMethod("d", signature = c(distr = "Multigam", x = "numeric"),
          function(distr, x, log = FALSE) {
            dmultigam(x, shape = distr@shape, scale = distr@scale, log = log)
          })

#' @rdname Multigam
setMethod("d", signature = c(distr = "Multigam", x = "matrix"),
          function(distr, x, log = FALSE) {
            dmultigam(x, shape = distr@shape, scale = distr@scale, log = log)
          })

#' @rdname Multigam
setMethod("r", signature = c(distr = "Multigam", n = "numeric"),
          function(distr, n) {
            rmultigam(n, shape = distr@shape, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
setMethod("mean",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  cumsum(x@shape * x@scale)

})

#' @rdname Multigam
setMethod("var",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  cumsum(x@shape * x@scale ^ 2)

})

#' @rdname Multigam
setMethod("finf",
          signature  = c(x = "Multigam"),
          definition = function(x) {

  # Preliminaries
  a <- x@shape
  b <- x@scale
  k <- length(a)
  a0 <- sum(a)

  D <- rbind(cbind(diag(trigamma(a)),
                   matrix(1 / b, nrow = k, ncol = 1)),
             c(rep(1 / b, k), a0 / b ^ 2))

  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
#' @export
llmultigam <- function(x, shape, scale) {
  ll(Multigam(shape, scale), x)
}

#' @rdname Multigam
setMethod("ll",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x) {

  k <- length(distr@shape)
  sum(apply(x, MARGIN = 1, FUN = dmultigam,
            shape = distr@shape, scale = distr@scale, log = TRUE))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Multigam"),
          definition = function(par, tx, distr) {

  k <- length(tx) - 1
  logz <- tx[1:k]
  xk <- tx[k + 1]

  b <- xk / par
  a <- idigamma(logz - log(b))

  - sum(a) * log(b) - sum(lgamma(a)) - xk / b + sum((a - 1) * logz)

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Multigam"),
          definition = function(par, tx, distr) {

  k <- length(tx) - 1
  logz <- tx[1:k]
  xk <- tx[k + 1]

  b <- xk / par
  a <- idigamma(logz - log(b))

  db <- - xk / par ^ 2
  da <- 1 / (par * trigamma(a))

  - sum(da) * log(b) - sum(a) * db / b - sum(digamma(a) * da) +
    xk * db / b ^ 2 + sum(logz * da)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
#' @export
emultigam <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me", "same"))
  distr <- Multigam()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Multigam
setMethod("mle",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  par0 <- check_optim(par0, method, lower, upper,
                      choices = c("me", "same"), len = 1)

  if (is.character(par0)) {
    par0 <- sum(unlist(emultigam(x, type = par0)$shape))
  }

  k <- ncol(x)
  logz <- colMeans(log(fd(x)))
  xk <- mean(x[, k])
  tx <- c(logz, xk)

  par <- optim(par = par0,
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  b <- xk / par
  a <- idigamma(logz - log(b))

  list(shape = a, scale = b)

})

#' @rdname Multigam
setMethod("me",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  z <- fd(x)
  mz <- colMeans(z)
  scale <- mean(colVar(z) / mz)
  shape <- mz / scale

  list(shape = shape, scale = scale)

})

#' @rdname Multigam
setMethod("same",
          signature  = c(distr = "Multigam", x = "matrix"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  z <- fd(x)
  scale <- mean(diag(stats::cov(z, log(z))))
  shape <- colMeans(z) / scale

  list(shape = shape, scale = scale)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Multigam
#' @export
vmultigam <- function(shape, scale, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me", "same"))
  distr <- Multigam(shape, scale)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Multigam
setMethod("avar_mle",
          signature  = c(distr = "Multigam"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  trinv <- 1 / trigamma(a)
  cons <- a0 - sum(trinv)

  D <- diag(cons * trinv) + Matrix(trinv, k, 1) %*% Matrix(trinv, 1, k)
  D <- cbind(D, - b * trinv)
  D <- rbind(D, t(c(- b * trinv, b ^ 2))) / cons

  D <- as.matrix(nearPD(D))
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

#' @rdname Multigam
setMethod("avar_me",
          signature  = c(distr = "Multigam"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  # Matrix A
  A11 <- (Matrix(a, k, 1) %*% Matrix(2 + 1 / a, 1, k) / k + diag(1, k, k)) / b
  A12 <- - Matrix(a, k, 1) %*% Matrix(1 / a, 1, k) / (k * b ^ 2)
  A21 <- - (2 + 1 / a) / k
  A22 <- 1 / (a * k * b)
  A <- rbind(cbind(A11, A12), Matrix(c(A21, A22), 1, 2 * k))

  # Matrix B
  B11 <- a * b ^ 2
  B22 <- 2 * a * (a + 1) * b ^ 4 * ( 2 * a + 3)
  B12 <- 2 * a * (a + 1) * b ^ 3
  B <- rbind(cbind(diag(B11), diag(B12)),
             cbind(diag(B12), diag(B22)))
  B <- nearPD(B)

  # Matrix D
  D <- nearPD(A %*% B %*% Matrix::t(A))
  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})

#' @rdname Multigam
setMethod("avar_same",
          signature  = c(distr = "Multigam"),
          definition = function(distr) {

  # Preliminaries
  a <- distr@shape
  b <- distr@scale
  k <- length(a)
  a0 <- sum(a)

  # Matrix A
  A12 <- Matrix(a, k, 1) %*% Matrix(a, 1, k) / k
  A13 <- - Matrix(a, k, 1) %*% Matrix(1, 1, k) / (k * b)
  A21 <- - (digamma(a) + log(b)) / k
  A22 <- - a * b / k
  A23 <- rep(1 / k, k)
  A11 <- (- Matrix(a, k, 1) %*% Matrix(A21, 1, k) + diag(1, k, k)) / b
  A <- rbind(cbind(A11, A12, A13), Matrix(c(A21, A22, A23), 1, 3 * k))

  # Matrix B
  B11 <- a * b ^ 2
  B22 <- trigamma(a)
  B33 <- a * (a + 1) * b ^ 2 *
    (trigamma(a + 2) + (digamma(a + 2) + log(b)) ^ 2) -
    (a * b) ^ 2 * (digamma(a + 1) + log(b)) ^ 2
  B12 <- rep(b, k)
  B13 <- a * (a + 1) * b ^ 2 * (digamma(a + 2) + log(b)) -
    (a * b) ^ 2 * (digamma(a + 1) + log(b))
  B23 <- a * b * (trigamma(a + 1) + (digamma(a + 1) + log(b)) ^ 2) -
    a * b * (digamma(a) + log(b)) * (digamma(a + 1) + log(b))
  B <- rbind(cbind(diag(B11), diag(B12), diag(B13)),
             cbind(diag(B12), diag(B22), diag(B23)),
             cbind(diag(B13), diag(B23), diag(B33)))
  B <- nearPD(B)

  # Matrix D
  D <- nearPD(A %*% B %*% Matrix::t(A))
  D <- as.matrix(D)
  rownames(D) <- c(paste0("shape", seq_along(a)), "scale")
  colnames(D) <- c(paste0("shape", seq_along(a)), "scale")
  D

})
