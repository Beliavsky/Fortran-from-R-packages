#' Fast Positive Definiteness Check
#'
#' @description
#' Efficiently checks if a matrix is positive definite using Cholesky decomposition.
#'
#' @param M Matrix. The matrix to verify.
#'
#' @return Logical. \code{TRUE} if positive definite, \code{FALSE} otherwise.
#' @export
is_pd_fast <- function(M) {
  if (is.null(M) || any(is.na(M))) return(FALSE)
  !inherits(try(chol(M), silent = TRUE), "try-error")
}