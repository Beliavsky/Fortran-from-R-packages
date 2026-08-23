# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cauchy Distribution                                                       ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Distribution           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("Cauchy",
         contains = "Distribution",
         slots = c(location = "numeric", scale = "numeric"),
         prototype = list(location = 0, scale = 1))

#' @title Cauchy Distribution
#' @name Cauchy
#'
#' @description
#' The Cauchy distribution is an absolute continuous probability distribution
#' characterized by its location parameter \eqn{x_0} and scale parameter
#' \eqn{\gamma > 0}.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Cauchy`.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Cauchy`. For the
#' log-likelihood and the estimation functions, `x` is the sample of
#' observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param location,scale numeric. Location and scale parameters.
#' @param type character, case ignored. The estimator type (mle or me).
#' @param log,log.p logical. Should the logarithm of the probability be
#' returned?
#' @param lower.tail logical. If TRUE (default), probabilities are
#' \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.
#' @param na.rm logical. Should the `NA` values be removed?
#' @param ... extra arguments.
#' @param par0,method,lower,upper arguments passed to optim for the mle
#' optimization.
#'
#' @details
#' The probability density function (PDF) of the Cauchy distribution is given
#' by: \deqn{ f(x; x_0, \gamma) = \frac{1}{\pi \gamma \left[1 + \left(\frac{x -
#' x_0}{\gamma}\right)^2\right]}.}
#'
#' The MLE of the Cauchy distribution parameters is not available in closed form
#' and has to be approximated numerically. This is done with `optim()`.
#' The default method used is the L-BFGS-B method with lower bounds
#' `c(-Inf, 1e-5)` and upper bounds `c(Inf, Inf)`. The `par0` argument can
#' either be a numeric (both elements satisfying `lower <= par0 <= upper`)
#' or a character specifying the closed-form estimator to be used as
#' initialization for the algorithm (`"me"` - the default value).
#'
#' Note that the `me()` estimator for the Cauchy distribution is not a
#' \emph{moment} estimator; it utilizes the sample median instead of the sample
#' mean.
#'
#' @inherit distributions return
#'
#' @seealso
#' Functions from the `stats` package: [dcauchy()], [pcauchy()], [qcauchy()],
#' [rcauchy()]
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Cauchy Distribution Example
#' # -----------------------------------------------------
#'
#' # Create the distribution
#' x0 <- 3 ; scale <- 5
#' D <- Cauchy(x0, scale)
#'
#' # ------------------
#' # dpqr Functions
#' # ------------------
#'
#' d(D, c(-5, 3, 10)) # density function
#' p(D, c(-5, 3, 10)) # distribution function
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
#' median(D) # Median
#' mode(D) # Mode
#' entro(D) # Entropy
#' finf(D) # Fisher Information Matrix
#'
#' # ------------------
#' # Point Estimation
#' # ------------------
#'
#' ll(D, x)
#' llcauchy(x, x0, scale)
#'
#' ecauchy(x, type = "mle")
#' ecauchy(x, type = "me")
#'
#' mle(D, x)
#' me(D, x)
#' e(D, x, type = "mle")
#'
#' mle("cauchy", x) # the distr argument can be a character
#'
#' # ------------------
#' # Estimator Variance
#' # ------------------
#'
#' vcauchy(x0, scale, type = "mle")
#' avar_mle(D)
#' v(D, type = "mle")
Cauchy <- function(location = 0, scale = 1) {
  new("Cauchy", location = location, scale = scale)
}

setValidity("Cauchy", function(object) {
  if(length(object@location) != 1) {
    stop("location has to be a numeric of length 1")
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

#' @rdname Cauchy
setMethod("d", signature = c(distr = "Cauchy", x = "numeric"),
          function(distr, x, log = FALSE) {
          dcauchy(x, location = distr@location, scale = distr@scale, log = log)
          })

#' @rdname Cauchy
setMethod("p", signature = c(distr = "Cauchy", q = "numeric"),
          function(distr, q, lower.tail = TRUE, log.p = FALSE) {
            pcauchy(q, location = distr@location, scale = distr@scale,
                    lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Cauchy
setMethod("qn", signature = c(distr = "Cauchy", p = "numeric"),
          function(distr, p, lower.tail = TRUE, log.p = FALSE) {
            qcauchy(p, location = distr@location, scale = distr@scale,
                    lower.tail = lower.tail, log.p = log.p)
          })

#' @rdname Cauchy
setMethod("r", signature = c(distr = "Cauchy", n = "numeric"),
          function(distr, n) {
            rcauchy(n, location = distr@location, scale = distr@scale)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
setMethod("mean",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The mean of the Cauchy distribution is not defined. NaN is returned")
  NaN

})


#' @rdname Cauchy
setMethod("median",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  x@location

})

#' @rdname Cauchy
setMethod("mode",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  x@location

})

#' @rdname Cauchy
setMethod("var",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The variance of the Cauchy distribution is not defined.
          NaN is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("sd",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The standard deviation of the Cauchy distribution is not
          defined. NaN is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("skew",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The skewness of the Cauchy distribution is not defined. NaN
          is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("kurt",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  warning("The kurtosis of the Cauchy distribution is not defined. NaN
          is returned.")
  NaN

})

#' @rdname Cauchy
setMethod("entro",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  log(4 * pi * x@scale)

})

#' @rdname Cauchy
setMethod("finf",
          signature  = c(x = "Cauchy"),
          definition = function(x) {

  D <- diag(2) / (2 * x@scale ^ 2)

  rownames(D) <- c("location", "scale")
  colnames(D) <- c("location", "scale")
  D

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
#' @export
llcauchy <- function(x, location, scale) {
  ll(distr = Cauchy(location, scale), x)
}

#' @rdname Cauchy
setMethod("ll",
          signature  = c(distr = "Cauchy", x = "numeric"),
          definition = function(distr, x) {


  - length(x) * log(pi * distr@scale) -
      sum(log (1 + ((x - distr@location) / distr@scale) ^ 2))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Cauchy"),
          definition = function(par, tx, distr) {

  log(par[2]) - mean(log((tx - par[1]) ^ 2 + par[2] ^ 2))

})

setMethod("dlloptim",
          signature  = c(par = "numeric", tx = "numeric", distr = "Cauchy"),
          definition = function(par, tx, distr) {

  c(2 * mean((tx - par[1]) / ((tx - par[1]) ^ 2 + par[2] ^ 2)),
    1 / par[2] - 2 * par[2] * mean(1 / ((tx - par[1]) ^ 2 + par[2] ^ 2)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
#' @export
ecauchy <- function(x, type = "mle", ...) {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Cauchy()
  do.call(type, list(distr = distr, x = x, ...))
}

#' @rdname Cauchy
setMethod("mle",
          signature  = c(distr = "Cauchy", x = "numeric"),
          definition = function(distr, x,
                                par0 = "me",
                                method = "L-BFGS-B",
                                lower = c(-Inf, 1e-5),
                                upper = c(Inf, Inf), na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)
  par0 <- check_optim(par0, method, lower, upper,
                      choices = c("me"), len = 2)

  if (is.character(par0)) {
    par0 <- unlist(ecauchy(x, type = par0))
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

  par <- unname(par)
  list(location = par[1], scale = par[2])

})

#' @rdname Cauchy
#' @export
setMethod("me",
          signature  = c(distr = "Cauchy", x = "numeric"),
          definition = function(distr, x, na.rm = FALSE) {

  x <- check_data(x, na.rm = na.rm)

  x0 <- median(x)

  list(location = x0, scale = median(abs(x - x0)))

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @rdname Cauchy
#' @export
vcauchy <- function(location, scale, type = "mle") {
  type <- match.arg(tolower(type), choices = c("mle", "me"))
  distr <- Cauchy(location, scale)
  do.call(paste0("avar_", type), list(distr = distr))
}

#' @rdname Cauchy
setMethod("avar_mle",
          signature  = c(distr = "Cauchy"),
          definition = function(distr) {

  inv2x2(finf(distr))

})
