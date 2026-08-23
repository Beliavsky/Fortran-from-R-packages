# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gam Distribution                                                          ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Gam",
         contains = "Distribution",
         slots = c(shape = "numeric", scale = "numeric"),
         prototype = list(shape = 1, scale = 1))

#' @title Gamma Distribution
#' @name Gam
#'
#' @description
#' The Gamma distribution is an absolute continuous probability distribution
#' with two parameters: shape \eqn{\alpha > 0} and scale \eqn{\beta > 0}.
#'
#' @srrstats {G1.0} A list of publications is provided in references.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Gam`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Gam`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param shape,scale numeric. The non-negative distribution parameters.
#' @param type character, case ignored. The estimator type (mle, me, or same).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization. See Details.
#'
#' @details
#' The probability density function (PDF) of the Gamma distribution is given by:
#' \deqn{ f(x; \alpha, \beta) = \frac{\beta^{-\alpha} x^{\alpha-1}
#' e^{-x/\beta}}{\Gamma(\alpha)}, \quad x > 0. }
#'
#' The MLE of the gamma distribution parameters is not available in closed form
#' and has to be approximated numerically. This is done with `optim()`. The
#' optimization can be performed on the shape parameter
#' \eqn{\alpha\in(0,+\infty)}. The default method used is the L-BFGS-B method
#' with lower bound `1e-5` and upper bound `Inf`. The `par0` argument can either
#' be a numeric (satisfying `lower <= par0 <= upper`) or a character specifying
#' the closed-form estimator to be used as initialization for the algorithm
#' (`"me"` or `"same"` - the default value).
#'
#' @inherit distributions return
#'
#' @references
#'
#' - Wiens, D. P., Cheng, J., & Beaulieu, N. C. (2003). A class of method of
#' moments estimators for the two-parameter gamma family. Pakistan Journal of
#' Statistics, 19(1), 129-141.
#'
#' - Ye, Z. S., & Chen, N. (2017). Closed-form estimators for the gamma
#' distribution derived from likelihood equations. The American Statistician,
#' 71(2), 177-181.
#'
#' - Tamae, H., Irie, K. & Kubokawa, T. (2020), A score-adjusted approach to
#' closed-form estimators for the gamma and beta distributions, Japanese Journal
#' of Statistics and Data Science 3, 543–561.
#'
#' - Papadatos, N. (2022), On point estimators for gamma and beta distributions,
#' arXiv preprint arXiv:2205.10799.
#'
#' @seealso
#' Functions from the `stats` package: [dgamma()], [pgamma()], [qgamma()],
#' [rgamma()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Gamma Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' a <- 3 ; b <- 5
#' D <- Gam(a, b)
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
#' ll(D, x)
#' llgamma(x, a, b)
#'
#' egamma(x, type = "mle")
#' egamma(x, type = "me")
#' egamma(x, type = "same")
#'
#' mle(D, x)
#' me(D, x)
#' same(D, x)
#' e(D, x, type = "mle")
#'
#' mle("gam", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vgamma(a, b, type = "mle")
#' vgamma(a, b, type = "me")
#' vgamma(a, b, type = "same")
#'
#' avar_mle(D)
#' avar_me(D)
#' avar_same(D)
#'
#' v(D, type = "mle")
Gam <- function(shape = 1, scale = 1) {
  new("Gam", shape = shape, scale = scale)
}

setValidity("Gam", function(object) {
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

#' @rdname Gam
setMethod("d", signature = c(distr = "Gam", x = "numeric"),
          function(distr, x, log = FALSE) {
            dgamma(x, shape = distr@shape, scale = distr@scale, log = log)
          })

#' @rdname Gam
setMethod("p", signature = c(distr = "Gam", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pgamma(q, shape = distr@shape, scale = distr@scale,
                   lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Gam
setMethod("qn", signature = c(distr = "Gam", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qgamma(p, shape = distr@shape, scale = distr@scale,
                   lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Gam
setMethod("r", signature = c(distr = "Gam", n = "numeric"),
          function(distr, n) {
            rgamma(n, shape = distr@shape, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
setMethod("mean",
          signature  = c(x = "Gam"),
          definition = function(x) {

  x@shape * x@scale

})

#' @rdname Gam
#' @export
setMethod("median",
          signature  = c(x = "Gam"),
          definition = function(x) {

  qgamma(0.5, shape = x@shape, scale = x@scale)

})

#' @rdname Gam
setMethod("mode",
          signature  = c(x = "Gam"),
          definition = function(x) {

  a <- x@shape
  b <- x@scale
  if(a >= 1) {
    return((a - 1) / b)
  } else {
    return(0)
  }

})

#' @rdname Gam
setMethod("var",
          signature  = c(x = "Gam"),
          definition = function(x) {

  x@shape * x@scale ^ 2

})

#' @rdname Gam
setMethod("sd",
          signature  = c(x = "Gam"),
          definition = function(x) {

  sqrt(var(x))

})

#' @rdname Gam
setMethod("skew",
          signature  = c(x = "Gam"),
          definition = function(x) {

  2 / sqrt(x@shape)

})

#' @rdname Gam
setMethod("kurt",
          signature  = c(x = "Gam"),
          definition = function(x) {

  6 / x@shape

})

#' @rdname Gam
setMethod("entro",
          signature  = c(x = "Gam"),
          definition = function(x) {

  a <- x@shape
  a + log(x@scale) + lgamma(a) + (1 - a) * digamma(a)

})

#' @rdname Gam
setMethod("finf",
          signature  = c(x = "Gam"),
          definition = function(x) {

  a <- x@shape
  b <- x@scale

  D <- matrix(c(trigamma(a), 1 / b, 1 / b, a / b ^ 2),
              nrow = 2, ncol = 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
#' @export
llgamma <- function(x, shape, scale) {
  ll(Gam(shape, scale), x)
}

#' @rdname Gam
setMethod("ll",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x) {

  a <- distr@shape
  b <- distr@scale
  - length(x) * (lgamma(a) + a * log(b)) + (a - 1) * sum(log(x)) - sum(x) / b

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Gam"),
          definition = function(par, tx, distr) {

  par * log(par) - lgamma(par) - (tx[1] + 1) * par + (par - 1) * tx[2]

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Gam"),
          definition = function(par, tx, distr) {

  log(par) - digamma(par) - tx[1] + tx[2]

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
#' @export
egamma <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me", "same"))
  distr <- Gam()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Gam
setMethod("mle",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x,
                                par0 = "same",
                                method = "L-BFGS-B",
                                lower = 1e-5,
                                upper = Inf, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  par0 <- check_optim(par0, method, lower, upper,
                      choices = c("me", "same"), len = 1)

  if (is.character(par0)) {
    par0 <- egamma(x, type = par0)$shape
  }

  tx <- c(log(mean(x)), mean(log(x)))

  par <- optim(par = par0,
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = distr,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  list(shape = par, scale = mean(x) / par)

})

#' @rdname Gam
setMethod("me",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  m  <- mean(x)
  m2 <- mean(x ^ 2)
  s2 <- m2 - m ^ 2

  list(shape = m ^ 2 / s2, scale = s2 / m)

})

#' @rdname Gam
setMethod("same",
          signature  = c(distr = "Gam", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  mx  <- mean(x)
  mlx <- mean(log(x))
  mxlx <- mean(x * log(x))
  cxlx <- mxlx - mx * mlx

  list(shape = mx / cxlx, scale = cxlx)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Gam
#' @export
vgamma <- function(shape, scale, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me", "same"))
  distr <- Gam(shape, scale)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Gam
setMethod("avar_mle",
          signature  = c(distr = "Gam"),
          definition = function(distr) {

  inv2x2(finf(distr))

})

#' @rdname Gam
setMethod("avar_me",
          signature  = c(distr = "Gam"),
          definition = function(distr) {

  a <- distr@shape
  b <- distr@scale

  s11 <- 2 * a * (a + 1)
  s22 <- b ^ 2 * (2 * a + 3) / a
  s12 <- - 2 * b * (a + 1)
  D <- matrix(c(s11, s12, s12, s22), nrow = 2, ncol = 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})

#' @rdname Gam
setMethod("avar_same",
          signature  = c(distr = "Gam"),
          definition = function(distr) {

  a <- distr@shape
  b <- distr@scale

  c1 <- 1 + a * trigamma(a + 1)
  c2 <- 1 + a * trigamma(a)

  v11 <- a ^ 2 * c1
  v21 <- - a * b * c1
  v22 <- b ^ 2 * c2

  D <- matrix(c(v11, v21, v21, v22), 2, 2)

  rownames(D) <- c("shape", "scale")
  colnames(D) <- c("shape", "scale")

  D

})
