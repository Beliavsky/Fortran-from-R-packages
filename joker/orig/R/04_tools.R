# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tools                                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Checks                 ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Check the Data
#'
#' @description
#' This function checks that the data argument supplied by the user for
#' parameter estimation are of the appropriate type.
#'
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a, G2.2} Assertions on the length and type
#' of input is implemented via the S4 methods along with additional checks made
#' by internal functions.
#' @srrstats {G2.13, G2.14, G2.14a, G2.14b, G2.14c} Estimation functions allow
#' users to specify whether missing (`NA`) data are to be ignored or produce an
#' error.
#' @srrstats {G2.15, G2.16} Estimation functions always check for `NaN` and
#' `Inf` values, which produce an error.
#' @srrstats {G5.3} The `check_data()` function checks the data for missing and
#' undefined values, and tests verifying that these values are appropriately
#' handled are included in the test files.
#'
#' @param x an object to be checked. Should be vector, matrix, or array.
#' @param na.rm logical. Should the NA values be removed?
#'
#' @returns The object `x`, possibly without `NA` values if `x` is a vector
#' containing `NA` values and `na.rm = TRUE`.
#'
#' @keywords internal
check_data <- function(x, na.rm = FALSE) {

  if (!is.logical(na.rm) || length(na.rm) > 1) {
    stop("na.rm must be a logical of length 1.")
  }
  if (!is_numatvec(x) && !is_nummat(x) && !is_numarray(x)) {
    stop("x must be a numeric vector, matrix, or array")
  }
  if (any(is.infinite(x))) {
    stop("x cannot include Inf values")
  }
  if (any(is.nan(x))) {
    stop("x cannot include NaN values")
  }
  if (any(is.na(x))) {
    if (!na.rm) {
      warning("Found NA values with na.rm = FALSE")
    } else if (!is.vector(x)) {
      stop("Found matrix/array with NA values. Cannot use na.rm = TRUE")
    } else {
      x <- x[!is.na(x)]
    }
  }
  x
}

#' @title Check Optim Arguments
#'
#' @description
#' Checks that the arguments supplied by the user are appropriate to be passed
#' to `optim()`. Used internally in parameter estimation.
#'
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a, G2.2} Assertions on the length and type
#' of input is implemented via the S4 methods along with additional checks made
#' by internal functions.
#'
#' @param par0 numeric or character. If numeric, it is passed to optim as the
#' initial estimation, i.e. the `par` argument. If character, the corresponding
#' estimation method is called and the result is passed to optim instead.
#' @param method,lower,upper arguments passed to optim.
#' @param choices character. A vector of available estimation methods for the
#' `par0` argument
#' @param len integer. The appropriate length of the `lower` and `upper`
#' argument, as well as `par` if it is numeric.
#'
#' @returns `par0`, possibly altered via `match.arg()` if it is a character.
#' @keywords internal
check_optim <- function(par0, method, lower, upper, choices = NULL, len = 1) {

  if (length(lower) != len || length(upper) != len) {
    stop("lower and upper must have length ", len)
  }
  if (!is.character(par0) && !is.numeric(par0)) {
    stop("par0 must be a character or a numeric")
  }
  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop("lower and upper must be numeric")
  }

  if (is.numeric(par0)) {
    if (any(par0 < lower) || any(par0 > upper)) {
      stop("par0 must be within the lower and upper bounds")
    }
    if (length(par0) != len) {
      stop("par0 must have length ", len)
    }
  } else {
    par0 <- match.arg(tolower(par0), choices = choices)
  }

  method <- match.arg(tolower(method),
                      choices = c("nelder-mead", "bfgs", "cg", "l-bfgs-b",
                                  "sann", "brent"))

  par0
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Progress Bar           ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Progress Bar
#' @name progress
#' @aliases progress_bar
#'
#' @description
#' Create a progress bar to be used with for loops that can possibly take a lot
#' of time.
#'
#' @param iter integer. The current iteration step of the for loop.
#' @param total integer. The total number of iterations.
#' @param start POSIXct. The start time, as returned by `Sys.time()`.
#' @param message character. A message appearing before the progress bar.
#' @param width integer. The length of the progress bar.
#' @param seconds integer. Seconds to be converted into hh:mm:ss format.
#'
#' @returns
#' `format_hms()` returns a character in the hh:mm:ss format. `progress_bar()`
#' prints the progress bar on the console, calling `cat()`, therefore as returns
#' an invisible `NULL`.
#'
#' @importFrom utils flush.console
#'
#' @keywords internal
progress_bar <- function(iter, total, start, message = NULL, width = 20) {

  percent <- iter / total
  hashes <- floor(percent * width)
  spaces <- width - hashes
  bar <- paste0(paste(rep("=", hashes), collapse = ""),
                paste(rep(" ", spaces), collapse = ""))

  elapsed_sec <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  eta_sec <- elapsed_sec / iter * (total - iter)

  cat("\r", message, sprintf("[%s] %3.0f%% | Elapsed: %s | ETA: %s",
                             bar, percent * 100, format_hms(elapsed_sec),
                             format_hms(eta_sec)))
  flush.console()

  if (iter == total) {
    cat("\n")
  }
}

#' @rdname progress
format_hms <- function(seconds) {
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  secs <- round(seconds %% 60)
  sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Structures             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Distribution Handling Helpers
#' @name distrhelpers
#'
#' @description
#' This set of functions help handle the distribution classes. See Details.
#'
#' @param x,distr an object of class `Distribution`.
#' @param list logical. Should a list be returned? If `FALSE`, the return object
#' is wrapped in `unlist()`.
#' @param prm,i A list containing three elements (`name`, `pos`, `val`) and the
#' i-th element of `val` to be updated as the new parameter. See
#' `small_metrics()`.
#'
#' @returns Depends on the function. See Details.
#'
#' @seealso [small_metrics()], [large_metrics()]
#'
#' @keywords internal
NULL

#' @describeIn distrhelpers Returns a character vector with the available moment
#' methods for the distribution.
get_moment_methods <- function(x) {

  # All available moments
  mom <- c("mean", "median", "mode", "var", "sd", "skew", "kurt",
           "entro", "finf")

  # Get class methods
  df_meth <- attr(methods(class = class(x)), "info")
  meth <- df_meth[df_meth$from == "joker", ]$generic

  mom[mom %in% meth]

}

#' @describeIn distrhelpers Turns the S4 class in the name (character) used in
#' the usual `stats` dpqr syntax.
get_class_abbr <- function(distr) {

  y <- tolower(class(distr)[1])

  switch(y,
         "gam" = "gamma",
         "stud" = "t",
         "fisher" = "f",
         "weib" = "weibull",
         y)

}

#' @describeIn distrhelpers Turns the distribution name from a character to an
#' S4 class.
get_distr_class <- function(distr) {

  distr <- paste(toupper(substr(distr, 1, 1)),
                 substr(tolower(distr), 2, nchar(distr)), sep = "")

  if (distr == "Gamma") {
    distr <- "Gam"
  }

  new(Class = distr)

}

#' @describeIn distrhelpers Turns an S4 distr object to a list.
s4_to_list <- function(distr) {

  # Get the slot names
  names <- methods::slotNames(class(distr))

  # Initialize an empty list to store slot values
  y <- list()

  # Loop through the slot names and extract slot values
  for (name in names) {
    y[[name]] <- methods::slot(distr, name)
  }

  y

}

#' @describeIn distrhelpers Get the parameters of a distribution as a list.
get_params <- function(distr, list = TRUE) {
  params <- s4_to_list(distr)
  params["name"] <- NULL
  params["ncp"] <- NULL

  if (!list) {
    params <- unlist(params)
  }

  params
}

#' @describeIn distrhelpers Get the unknown parameters of a distribution as a
#' list.
get_unknown_params <- function(distr, list = TRUE) {

  params <- get_params(distr)

  if (is(distr, "Binom") || is(distr, "Nbinom") || is(distr, "Multinom")) {
    params$size <- NULL
  }

  if (!list) {
    params <- unlist(params)
  }

  params
}

#' @describeIn distrhelpers Update the distribution parameters. Returns the
#' distribution object. Used inside the `small_metrics()` and `large_metrics()`
#' functions.
update_params <- function(distr, prm, i) {

  # Position of parameter (e.g. the third element of the shape parameter vector)
  if (is.null(prm$pos)) {
    prm$pos <- 1
  }

  params <- s4_to_list(distr)
  params[[prm$name]][prm$pos] <- prm$val[i]

  do.call("new", c(params, Class = class(distr)))

}

#' Turn an array to a data.frame
#'
#' @description
#' This function turns an array to a data.frame. It is used by the
#' `small_metrics()` and `large_metrics()` functions.
#'
#' @param x array.
#'
#' @returns data.frame.
#'
#' @keywords internal
array_to_df <- function(x) {

  dn <- dimnames(x)
  names(dn) <- names(dimnames(x))

  df <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE)
  df$Value <- as.vector(x)

  df

}

#' Indexing Functions
#'
#' @description
#' A set of functions that allow to index a matrix or array. These functions are
#' used internally when the dimension length of an object can vary.
#'
#' @param x atomic vector, matrix, or array. An object to be indexed.
#' @param i integer. The index.
#'
#' @returns A vector or matrix, subset of the original object.
#'
#' @keywords internal
set1of1 <- function(x, i) {
  x
}

#' @rdname set1of1
set1of2 <- function(x, i) {
  x[i, ]
}

#' @rdname set1of1
set1of3 <- function(x, i) {
  x[i, , ]
}

#' @rdname set1of1
set2of3 <- function(x, i) {
  x[, i, ]
}

#' Sequence Functions for Matrices
#'
#' @description
#' This set of functions extent the `seq_along()` functions for `matrix`
#' objects.
#'
#' @param x matrix.
#'
#' @returns A sequence of integers from 1 to the number of rows or columns of
#' the matrix.
#'
#' @keywords internal
seqcol <- function(x) {
  if (!is.matrix(x)) {
    stop("x is not a matrix.")
  }
  seq_along(x[1, ])
}

#' @rdname seqcol
seqrow <- function(x) {
  if (!is.matrix(x)) {
    stop("x is not a matrix.")
  }
  seq_along(x[ , 1])
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Statistics             ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Statistical Functions
#' @name stats
#'
#' @description
#' A set of statistics that extend the ones made available in the `stats`
#' package.
#'
#' @param x numeric for `bvar()` and `bsd()`, matrix for `rowVar()` and
#' `colVar()`.
#'
#' @returns `bvar()` and `bsd()` return a single numeric, `rowVar()` and
#' `colVar()` return a numeric vector.
#'
#' @keywords internal
NULL

#' @describeIn stats Biased sample variance
bvar <- function(x) {
  if (!is_numatvec(x)) {
    stop("x is not a numeric matrix.")
  }
  ((length(x) - 1) / length(x)) * var(x)
}

#' @describeIn stats Biased sample standard deviation
bsd <- function(x) {
  if (!is_numatvec(x)) {
    stop("x is not a numeric matrix.")
  }
  sqrt(bvar(x))
}

#' @describeIn stats Biased sample variance by matrix row
rowVar <- function(x) {
  if (!is_nummat(x)) {
    stop("x is not a numeric matrix.")
  }
  rowMeans(x ^ 2) - rowMeans(x) ^ 2
}

#' @describeIn stats Biased sample variance by matrix column
colVar <- function(x) {
  if (!is_nummat(x)) {
    stop("x is not a numeric matrix.")
  }
  colMeans(x ^ 2) - colMeans(x) ^ 2
}

#' @describeIn stats Sample skewness
setMethod("skew",
          signature  = c(x = "numeric"),
          definition = function(x) {

  mean((x - mean(x)) ^ 3)

})

#' @describeIn stats Sample kurtosis
setMethod("kurt",
          signature  = c(x = "numeric"),
          definition = function(x) {

  mean((x - mean(x)) ^ 4)

})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Gamma Function         ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Inverse Digamma Function
#'
#' @description
#' The inverse of the digamma function, i.e. the derivative of the log-gamma
#' function.
#'
#' @param x numeric. The point to evaluate the function.
#' @param ... extra arguments passed to `optim()`.
#'
#' @return numeric. The evaluated function.
#'
#' @export
#'
#' @details
#' The `idigamma()` function implements the inverse of the digamma function
#' \eqn{\psi}. It is a numerical approximation based on the Brent optimization
#' algorithm. Specifically, `idigamma()` makes a call to `optim()` in order to
#' solve the equation \eqn{\psi(x) = y}; more accurately, to find the minimum of
#' \eqn{f(x) = \log\Gamma(x) - xy}, whose derivative is
#' \eqn{f'(x) = \psi(x) - y}. The optimization is restricted within the tight
#' bounds derived by Batir (2017). The function is vectorized.
#'
#' @references
#' Necdet Batir (2017), INEQUALITIES FOR THE INVERSES OF THE POLYGAMMA FUNCTIONS
#' https://arxiv.org/pdf/1705.06547
#'
#' Oikonomidis, I. & Trevezas, S. (2023), Moment-Type Estimators for the
#' Dirichlet and the Multivariate Gamma Distributions, arXiv,
#' https://arxiv.org/abs/2311.15025
#'
#' @seealso [optim()]
#'
#' @examples
#' \donttest{
#' idigamma(2)
#' }
idigamma <- function(x, ...) {

  unlist(lapply(x, FUN = function(x) {

    # Tight bounds derived in https://arxiv.org/pdf/1705.06547
    l <- 1 / log(1 + exp(- x))
    u <- exp(x) + 0.5

    # Quasi-Newton optimization
    optim(par = (l + u) / 2,
          fn = function(y) {lgamma(y) - x * y},
          #gr = function(y) {digamma(y) - x},
          method = "Brent",
          lower = l,
          upper = u, ...)$par
  }))

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Matrix Algebra         ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Matrix Wrappers
#'
#' @description
#' Simple wrappers of functions from the `Matrix` package.
#'
#' @param x matrix
#' @param ... extra arguments passed to `Matrix::Matrix()`.
#'
#' @returns matrix
#'
#' @keywords internal
Matrix <- function(...) {
  Matrix::Matrix(...)
}

#' @rdname Matrix
nearPD <- function(x) {
  Matrix::nearPD(x)$mat
}

#' 2x2 Inverse Matrix
#'
#' @description
#' Calculates the inverse of a 2x2 matrix.
#'
#' @param x A 2x2 matrix.
#'
#' @returns A 2x2 matrix, the inverse of `x`
#'
#' @keywords internal
inv2x2 <- function(x) {

  if (!is.matrix(x) || any(dim(x) != 2)) {
    stop("x must be a 2x2 matrix")
  }

  det <- x[1, 1] * x[2, 2] - x[1, 2] * x[2, 1]

  if (det == 0) {
    stop("The matrix is singular, its inverse does not exist.")
  } else {
    inv <- matrix(c(x[2, 2], -x[2, 1], -x[1, 2], x[1, 1]), nrow = 2)
    inv <- inv / det
    dimnames(inv) <- dimnames(x)
    return(inv)
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Is it?                 ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Is it?
#' @name isit
#'
#' @description
#' A set of functions that check whether an object has the desired
#' characteristics.
#'
#' @param x numeric vector or matrix.
#' @param tol numeric. The tolerance for a numeric to be considered a whole
#' number.
#' @returns logical. TRUE or FALSE, depending on whether the object satisfies
#' the checks that define the characteristic.
#'
#' @keywords internal
NULL

#' @describeIn isit Is the object an integer in the mathematical sense?
is_whole <- function(x, tol = .Machine$double.eps^0.5) {
  is.numeric(x) && all(abs(x - round(x)) < tol)
}

#' @describeIn isit Is the object an atomic vector (not matrix or array)?
is_atvec <- function(x) {
  is.atomic(x) && is.vector(x)
}

#' @describeIn isit Is the object a numeric atomic vector?
is_numatvec <- function(x) {
  is.numeric(x) && is.atomic(x) && is.vector(x)
}

#' @describeIn isit Is the object a numeric atomic matrix?
is_nummat <- function(x) {
  is.numeric(x) && is.matrix(x)
}

#' @describeIn isit Is the object a numeric atomic array?
is_numarray <- function(x) {
  is.numeric(x) && is.array(x)
}

#' @describeIn isit Is the object a symmetric matrix?
is_symmetric <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x is not a numeric matrix.")
  }
  sum(x == t(x)) == (nrow(x) ^ 2)
}

#' @describeIn isit Is the object a positive definite matrix?
is_pd <- function(x) {

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x is not a numeric matrix.")
  }

  if (!is_symmetric(x)) {
    return(FALSE)
  }

  eigs <- eigen(x, symmetric = TRUE)$values

  if (any(is.complex(eigs)) || any(eigs < 0)) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}

#' @describeIn isit Are all the elements finite and positive?
is_pos <- function(x) {
  all(is.finite(x)) && all(x > 0)
}

#' @describeIn isit Is the object an integer in the mathematical sense?
is_integer <- function(x) {
  identical(x, round(x))
}

#' @describeIn isit Is the object a natural number in the mathematical sense?
is_natural <- function(x) {
  is_integer(x) && (x > 0)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multivariate Gamma     ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Forward Difference
#'
#' @description
#' Calculates the forward difference of a vector or matrix.
#'
#' @param x numeric vector or matrix.
#'
#' @details
#' The function is used internally in the `Multigam` distribution.
#'
#' @returns an atomic vector or matrix of the same dimensions as `x`.
#'
#' @keywords internal
fd <- function(x) {
  if (!is_numatvec(x) && !is_nummat(x)) {
    stop("x must be a numeric atomic vector or matrix.")
  } else if (is.matrix(x)) {
    return(cbind(x[, 1], t(diff(t(x)))))
  } else {
    return(c(x[1], diff(x)))
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generalized Gamma      ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

stacy_to_prentice <- function(par) {
  a <- par$a
  b <- par$b
  k <- par$k
  list(mu = log(a) + digamma(k) / b,
       sigma = 1 / (b * sqrt(k)),
       q = 1 / sqrt(k))
}

prentice_to_stacy <- function(par) {
  mu <- par$mu
  sigma <- par$sigma
  q <- par$q
  k <- 1 / q ^ 2
  b <- q / sigma
  a <- exp(mu - digamma(k) / b)
  list(a = a, b = b, k = k)
}
