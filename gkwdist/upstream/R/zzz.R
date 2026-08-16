utils::globalVariables(c(
  "packageVersion", "dev.interactive", "calculateMean", "index", "theoretical", "observed", "lower",
  "upper", "log", "log_scale", "::", ":::", "x", "Theoretical", "Empirical", "value", "loglik", "cook_dist",
  "fitted", "abs_resid", "leverage", "y_obs", "linpred", "resid", "model_label", "metric", "y", "Type", "Deviation",
  "object", "p_empirical", "p_theoretical", "statistic", "type", "Residual", "Family", "Value", "Criterion", "Parameter"
))

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


## usethis namespace: start
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom numDeriv grad hessian
#' @import RcppArmadillo
#' @import graphics
## usethis namespace: end
NULL

#' @useDynLib gkwdist, .registration = TRUE
NULL
