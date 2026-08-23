# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generics & Classes                                                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Distribution S4 Classes
#' @name distributions
#' @aliases d p q r
#'
#' @description
#' A collection of S4 classes that provide a flexible and structured way to work
#' with probability distributions.
#'
#' @param n number of observations. If `length(n) > 1`, the length is taken to
#' be the number required.
#' @param distr an object of class `Distribution` or one of its subclasses.
#' @param x For the density function, `x` is a numeric vector of quantiles. For
#' the moments functions, `x` is an object of class `Distribution`  or one of
#' its subclasses. For the log-likelihood and the estimation functions, `x` is
#' the sample of observations.
#' @param p numeric. Vector of probabilities.
#' @param q numeric. Vector of quantiles.
#' @param ... extra arguments.
#'
#' @details
#' These S4 generic methods can work both as functions and as functionals
#' (functions that return functions). The available distribution families are
#' coded as S4 classes, specifically subclasses of the `Distribution`
#' superclass. The methods can be used in two ways:
#'
#' Option 1: If both the `distr` argument and `x` or `n` are supplied, then the
#' function is evaluated directly, as usual.
#'
#' Option 2: If only the `distr` argument is supplied, the method returns a
#' function that takes as input the missing argument `x` or `n`, allowing the
#' user to work with the function object itself. See examples.
#'
#' Looking for a specific distribution family?
#' This help page is general. Use the help page of each distribution to see the
#' available methods for the class, details, and examples. Check the See Also
#' section.
#'
#' @return
#' Each type of function returns a different type of object:
#'
#' - Distribution Functions: When supplied with one argument (`distr`), the
#' `d()`, `p()`, `q()`, `r()`, `ll()` functions return the density, cumulative
#' probability, quantile, random sample generator, and log-likelihood functions,
#' respectively. When supplied with both arguments (`distr` and `x`), they
#' evaluate the aforementioned functions directly.
#'
#' - Moments: Returns a numeric, either vector or matrix depending on the moment
#' and the distribution. The `moments()` function returns a list with all the
#' available methods.
#'
#' - Estimation: Returns a list, the estimators of the unknown parameters. Note
#' that in distribution families like the binomial, multinomial, and negative
#' binomial, the size is not returned, since it is considered known.
#'
#' - Variance: Returns a named matrix. The asymptotic covariance matrix of the
#' estimator.
#'
#' @seealso [moments], [loglikelihood], [estimation], [Bern],
#' [Beta], [Binom], [Cat], [Cauchy], [Chisq], [Dir], [Exp],
#' [Fisher], [Gam], [Geom], [Laplace], [Lnorm], [Multigam], [Multinom],
#' [Nbinom], [Norm], [Pois], [Stud], [Unif], [Weib]
#'
#' @inherit Beta examples
NULL

setClass("Distribution")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## d, p, q, r             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @describeIn distributions density function
setGeneric("d", function(distr, x, ...) {
  standardGeneric("d")
})

#' @describeIn distributions cumulative distribution function
setGeneric("p", function(distr, q, ...) {
  standardGeneric("p")
})

#' @describeIn distributions generalized inverse distribution function
setGeneric("qn", function(distr, p, ...)
  standardGeneric("qn"))

#' @describeIn distributions random sample generator function
setGeneric("r", function(distr, n, ...) {
  standardGeneric("r")
})

#' @title Distribution Functionals
#' @name functionals
#'
#' @param distr a `Distribution` object.
#' @param x,q,p,n missing. Arguments not supplied.
#' @param ... extra arguments.
#'
#' @details
#' When `x`, `q`, `p`, or `n` are missing, the methods return a function that
#' takes as input the missing argument, allowing the user to work with the
#' function object itself. See examples.
#'
#' @return
#' When supplied with one argument, the `d()`, `p()`, `q()`, `r()` `ll()`
#' functions return the density, cumulative probability, quantile, random sample
#' generator, and log-likelihood functions, respectively.
#'
#' @inherit distributions description seealso examples
NULL

#' @rdname functionals
setMethod("d", signature = c(distr = "Distribution", x = "missing"),
          function(distr, x, ...) {
            function(x) {d(distr, x, ...)}
          })

#' @rdname functionals
setMethod("p", signature = c(distr = "Distribution", q = "missing"),
          function(distr, q, ...) {
            function(q) {p(distr, q, ...)}
          })

#' @rdname functionals
setMethod("qn", signature = c(distr = "Distribution", p = "missing"),
          function(distr, p, ...) {
            function(p) {qn(distr, p, ...)}
          })

#' @rdname functionals
setMethod("r", signature = c(distr = "Distribution", n = "missing"),
          function(distr, n, ...) {
            function(n) {r(distr, n, ...)}
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Moments                ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Moments - Parametric Quantities of Interest
#' @name moments
#'
#' @description A set of functions that calculate the theoretical moments
#' (expectation, variance, skewness, excess kurtosis) and other important
#' parametric functions (median, mode, entropy, Fisher information) of a
#' distribution.
#'
#' @param x a `Distribution` object.
#' @param y,use,na.rm arguments in `mean` and `var` standard methods from the
#' `stats` package not used here.
#' @param ... extra arguments.
#'
#' @details
#' Given a distribution, these functions calculate the theoretical moments and
#' other parametric quantities of interest. Some distribution-function
#' combinations are not available; for example, the `sd()` function is
#' available only for univariate distributions.
#'
#' The `moments()` function automatically finds the available methods for a
#' given distribution and results all of the results in a list.
#'
#' Technical Note:
#' The argument of the moment functions does not follow the naming convention of
#' the package, i.e. the `Distribution` object is names `x` rather than `distr`.
#' This is due to the fact that most of the generics are already defined in the
#' `stats` package (`mean`, `median`, `mode`, `var`, `sd`), therefore the first
#' argument was already named `x` and could not change.
#'
#' @return Numeric, either vector or matrix depending on the moment and the
#' distribution. The `moments()` function returns a list with all the available
#' methods.
#'
#' @export
#'
#' @seealso [distributions], [loglikelihood], [estimation]
#'
#' @inherit distributions examples
moments <- function(x) {
  mom <- get_moment_methods(x)
  y <- lapply(mom, FUN = function(m) { do.call(m, list(x = x)) })
  names(y) <- mom
  y
}

#' @rdname moments
#' @name mean
#' @usage mean(x, ...)
NULL

#' @describeIn moments Median
setGeneric("median")

#' @describeIn moments Mode
setGeneric("mode")

#' @describeIn moments Variance
setGeneric("var")

#' @describeIn moments Standard Deviation
setGeneric("sd")

#' @describeIn moments Skewness
setGeneric("skew", function(x, ...) {
  standardGeneric("skew")
})

#' @describeIn moments Kurtosis
setGeneric("kurt", function(x, ...) {
  standardGeneric("kurt")
})

#' @describeIn moments Entropy
setGeneric("entro", function(x, ...) {
  standardGeneric("entro")
})

#' @describeIn moments Fisher Information (numeric or matrix)
setGeneric("finf", function(x, ...) {
  standardGeneric("finf")
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Likelihood             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Log-Likelihood Function
#' @name loglikelihood
#' @aliases ll
#'
#' @description
#' This function calculates the log-likelihood of an independent and
#' identically distributed (iid) sample from a distribution. See Details.
#'
#' @param distr A `Distribution` object
#' @param x numeric. A sample under estimation.
#' @param ... extra arguments.
#'
#' @details
#' The log-likelihood functions are provided in two forms: the `ll<name>`
#' distribution-specific version that follows the stats package conventions, and
#' the S4 generic `ll`. Examples for the `ll<name>` version are included in the
#' distribution-specific help pages, e.g. `?Beta` (all distributions can be
#' found in the See Also section of the `Distributions` help page).
#'
#' As with the `d()`, `p()`, `q()`, `r()` methods, `ll()` can be supplied only
#' with `distr` to return the log-likelihood function (i.e. it can be used as a
#' functional), or with both `distr` and `x` to be evaluated directly.
#'
#' In some distribution families like beta and gamma, the MLE cannot be
#' explicitly derived and numerical optimization algorithms have to be employed.
#' Even in ``good" scenarios, with plenty of observations and a smooth
#' optimization function, extra care should be taken to ensure a fast and right
#' convergence if possible. Two important steps are taken in package in this
#' direction:
#'
#' 1. The log-likelihood function is analytically calculated for each
#' distribution family, so that constant terms with respect to the parameters
#' can be removed, leaving only the sufficient statistics as a requirement for
#' the function evaluation.
#'
#' 2. Multidimensional problems are reduced to unidimensional ones by utilizing
#' the score equations.
#'
#' The resulting function that is inserted in the optimization algorithm is
#' called `lloptim()`, and is not to be confused with the actual log-likelihood
#' function `ll()`. The corresponding derivative is called `dlloptim()`. These
#' functions are used internally and are not exported.
#'
#' Therefore, whenever numerical computation of the MLE is required, `optim()`
#' is called to optimize `lloptim()`, using the ME or SAME as the starting point
#' (user's choice), and the L-BFGS-U optimization algorithm, with lower and
#' upper limits defined by default as the parameter space boundary. Illustrative
#' examples can be found in the package vignette.
#'
#' @return If only the `distr` argument is supplied, `ll()` returns a function.
#' If both `distr` and `x` are supplied, `ll()` returns a numeric, the value of
#' the log-likelihood function.
#'
#' @export
#'
#' @seealso [distributions], [moments], [estimation]
#'
#' @inherit distributions examples
setGeneric("ll", signature = c("distr", "x"),
           function(distr, x, ...) { standardGeneric("ll") })

#' @rdname functionals
setMethod("ll",
          signature  = c(distr = "Distribution", x = "missing"),
          definition = function(distr, x, ...) {
            function(x) {
              ll(distr, x, ...)
            }
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Score                  ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

setGeneric("lloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("lloptim") })

setGeneric("dlloptim", signature = c("par", "tx", "distr"),
           function(par, tx, distr, ...) { standardGeneric("dlloptim") })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Estimation             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Parameter Estimation
#' @name estimation
#' @aliases e mle me same
#'
#' @srrstats {G1.0} A list of publications is provided in references.
#' @srrstats {G1.3} The statistical terminology used in the package is
#' elementary and well-established in the literature. The abbreviations used for
#' the estimation methods (MLE, ME, SAME) are described in the documentation.
#'
#' @description
#' This set of functions estimates the parameters of a random sample according
#' to a specified family of distributions. See details.
#'
#' @param distr A `Distribution` object or a `character`. The distribution
#' family assumed.
#' @param x numeric. A sample under estimation.
#' @param type character, case ignored. The estimator type.
#' @param ... extra arguments.
#'
#' @details
#' The package covers three major estimation methods: maximum likelihood
#' estimation (MLE), moment estimation (ME), and score-adjusted estimation
#' (SAME).
#'
#' In order to perform parameter estimation, a new `e<name>()` member is added
#' to the `d()`, `p()`, `q()`, `r()` family, following the standard `stats` name
#' convention. These functions take two arguments, the observations `x` (an
#' atomic vector for univariate or a matrix for multivariate distributions) and
#' the `type` of estimation method to use (a character with possible values
#' `"mle"`, `"me"`, and `"same"`.)
#'
#' Point estimation functions are available in two versions, the distribution
#' specific one, e.g. `ebeta()`, and the S4 generic ones, namely `mle()`,
#' `me()`, and `same()`. A general function called `e()` is also implemented,
#' covering all distributions and estimators.
#'
#' @return list. The estimator of the unknown parameters. Note that in
#' distribution families like the binomial, multinomial, and negative binomial,
#' the size is not returned, since it is considered known.
#'
#' @importFrom Matrix Matrix nearPD Cholesky
#' @export
#'
#' @references
#'
#' General Textbooks
#'
#' - Van der Vaart, A. W. (2000), Asymptotic statistics, Vol. 3,
#' Cambridge university press.
#'
#' Beta and gamma distribution families
#'
#' - Ye, Z.-S. & Chen, N. (2017), Closed-form estimators for the gamma
#' distribution derived from likelihood equations, The American Statistician
#' 71(2), 177–181.
#'
#' - Tamae, H., Irie, K. & Kubokawa, T. (2020), A score-adjusted approach to
#' closed-form estimators for the gamma and beta distributions, Japanese Journal
#' of Statistics and Data Science 3, 543–561.
#'
#' - Mathal, A. & Moschopoulos, P. (1992), A form of multivariate gamma
#' distribution, Annals of the Institute of Statistical Mathematics 44, 97–106.
#'
#' - Oikonomidis, I. & Trevezas, S. (2023), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @seealso [mle], [me], [same]
#'
#' @inherit distributions examples
e <- function(distr, x, type = "mle", ...) {
  abbr <- get_class_abbr(distr)
  if (abbr %in% c("binom", "nbinom")) {
    args <- list(x = x, type = type, size = distr@size, ...)
  } else if (abbr == "cat") {
    args <- list(x = x, type = type, dim = length(distr@prob), ...)
  } else {
    args <- list(x = x, type = type, ...)
  }
  do.call(paste0("e", abbr), args)
}

#' @describeIn estimation Maximum Likelihood Estimator
setGeneric("mle", signature = c("distr", "x"),
           function(distr, x, ...) { standardGeneric("mle") })

#' @rdname functionals
setMethod("mle",
          signature  = c(distr = "Distribution", x = "missing"),
          definition = function(distr, x, ...) {
            function(x) {
              mle(distr, x, ...)
            }
          })

#' @rdname estimation
setMethod("mle",
          signature  = c(distr = "character", x = "ANY"),
          definition = function(distr, x, ...) {
            mle(get_distr_class(distr), x, ...)
          })

#' @describeIn estimation Moment Estimator
setGeneric("me", signature = c("distr", "x"),
           function(distr, x, ...) { standardGeneric("me") })

#' @rdname functionals
setMethod("me",
          signature  = c(distr = "Distribution", x = "missing"),
          definition = function(distr, x, ...) {
            function(x) {
              me(distr, x, ...)
            }
          })

#' @rdname estimation
setMethod("me",
          signature  = c(distr = "character", x = "ANY"),
          definition = function(distr, x, ...) {
            me(get_distr_class(distr), x, ...)
          })

#' @describeIn estimation Score - Adjusted Moment Estimation
setGeneric("same", signature = c("distr", "x"),
           function(distr, x, ...) { standardGeneric("same") })

#' @rdname functionals
setMethod("same",
          signature  = c(distr = "Distribution", x = "missing"),
          definition = function(distr, x, ...) {
            function(x) {same(distr, x, ...)}
          })

#' @rdname estimation
setMethod("same",
          signature  = c(distr = "character", x = "ANY"),
          definition = function(distr, x, ...) {
            same(get_distr_class(distr), x, ...)
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variance               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Estimator Variance
#' @name variance
#' @aliases v avar_mle avar_me avar_same
#'
#' @description
#' These functions calculate the variance (or variance - covariance matrix in
#' the multidimensional case) of an estimator, given a specified family of
#' distributions and the true parameter values.
#'
#' @param distr A `Distribution` object.
#' @param type character, case ignored. The estimator type.
#' @param ... extra arguments.
#'
#' @return numeric, or matrix for multidimensional cases.
#'
#' @export
#'
#' @seealso [avar_mle], [avar_me], [avar_same]
#'
#' @inherit estimation references examples
v <- function(distr, type, ...) {
  do.call(paste0("v", get_class_abbr(distr)),
          c(get_params(distr), list(type = type, ...)))
}

#' @describeIn variance Asymptotic Variance of the Maximum Likelihood Estimator
setGeneric("avar_mle", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_mle") })

#' @describeIn variance Asymptotic Variance of the Moment Estimator
setGeneric("avar_me", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_me") })

#' @describeIn variance Asymptotic Variance of the Score-Adjusted Moment
#' Estimator
setGeneric("avar_same", signature = c("distr"),
           function(distr, ...) { standardGeneric("avar_same") })
