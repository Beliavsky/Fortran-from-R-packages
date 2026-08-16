#' Sum of One-Dimensional Functions
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @param f an R function taking a numeric first argument and returning a numeric vector of the same length.
#' @param lower the lower limit of sum. Can be infinite.
#' @param upper the upper limit of sum. Can be infinite.
#' @param ... additional arguments to be passed to f.
#' @param abs.tol absolute accuracy requested.
#'
#' @return This function returns the sum value.
#'
#' @example  examples/examples_add.R
#'
#' @importFrom utils tail
#' @export
add <- function (f, lower, upper, ..., abs.tol = .Machine$double.eps) {

  f <- match.fun(f)
  ff <- function(x) f(x, ...)

  if (lower >= upper)
    stop("invalid parameter values")
  stopifnot(length(lower) == 1, length(upper) == 1)

  # My auxiliar functions ------------------------------------------------------

  # First function
  add_minusinf_to_inf <- function(ff, ..., abs.tol) {

    x <- seq(from=-100, to=100) #to ensure a sum with at least 201 values
    ans <- sum(ff(x))
    x <- tail(x, n=1L) + 1 # The next value

    while (TRUE) {
      next_term <- ff(x) + ff(-x)
      ans <- ans + next_term
      if (abs(next_term) < abs.tol) break
      x <- x + 1
    }
    list(value=ans, abs.error=abs(next_term))
  }

  # Second function
  add_lower_to_inf <- function(ff, lower, ..., abs.tol) {

    x <- seq(from=lower, to=lower+300) #to ensure a sum with at least 301 values
    ans <- sum(ff(x))
    x <- tail(x, n=1L) + 1 # The next value

    while (TRUE) {
      next_term <- ff(x)
      ans <- ans + next_term
      if (abs(next_term) < abs.tol) break
      x <- x + 1
    }
    list(value=ans, abs.error=abs(next_term))
  }

  # Third function
  add_minusinf_to_upper <- function(ff, upper, ..., abs.tol) {

    x <- seq(from=upper, to=upper-300) #to ensure a sum with at least 301 values
    ans <- sum(ff(x))
    x <- tail(x, n=1L) - 1 # The next value

    while (TRUE) {
      next_term <- ff(x)
      ans <- ans + next_term
      if (abs(next_term) < abs.tol) break
      x <- x - 1
    }
    list(value=ans, abs.error=abs(next_term))
  }

  # End my auxiliar functions --------------------------------------------------

  # Sum with finite lower and upper
  if (is.finite(lower) && is.finite(upper)) {
    wk <- sum(ff(seq(from=lower, to=upper, by=1)))
    wk <- list(value=wk, abs.error=0)
  }

  else {
    if (is.na(lower) || is.na(upper))
      stop("a limit is NA or NaN")

    if (is.finite(lower)) {
      wk <- add_lower_to_inf(ff, lower, ..., abs.tol=abs.tol)
    }
    else if (is.finite(upper)) {
      wk <- add_minusinf_to_upper(ff, upper, ..., abs.tol=abs.tol)
    }
    else {
      wk <- add_minusinf_to_inf(ff, ..., abs.tol=abs.tol)
    }
  }
  return(wk)
}
