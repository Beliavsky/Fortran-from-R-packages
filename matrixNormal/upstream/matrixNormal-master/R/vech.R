#' Half-Vectorization of a matrix
#'
#' @description Stacks elements of the lower triangle of a numeric symmetric matrix \emph{A}.
#'
#' @details
#' For a symmetric matrix \emph{A}, the vectorization of \emph{A} contains more information than necessary. The half-vectorization, denoted \code{vech()}, of a symmetric square n by n matrix \emph{A} is the vectorization of the lower triangular portion.
#'
#' @note
#' Unlike other \code{vech()} functions available on CRAN, matrixNormal version may inherit names from matrices to their vectorized forms.
#'
#' @inheritParams vec
#' @inheritParams is.symmetric.matrix
#' @return A vector with n(n+1)/2 elements.
#'
#' @export
#'
#' @examples
#' x <- matrix(c(1, 2, 2, 4),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(1:2, c("Sex", "Smoker"))
#' )
#' print(x)
#'
#' # Example 1
#' vech(x)
#' # If you just want the vectorized form
#' vech(x, use.Names = FALSE)
#'
#' # Example 2: If one has NA's
#' x[1, 2] <- x[2, 1] <- NA
#' vech(x)
vech <- function(A, use.Names = TRUE, tol = .Machine$double.eps^0.5) {
  if (is.vector(A)) {
    if (length(A) == 1) {
      return(A)
    } else {
      stop("vech undefined for vectors")
    }
  }

  symm <- is.symmetric.matrix(A, tol)
  if (isFALSE(symm)) {
    stop(sprintf("% must be a numeric and symmetric matrix for half-vectorization."))
  }
  # } else { #if symm is TRUE or NA ...

  full <- vec(A, use.Names)
  stack <- A[lower.tri(A, diag = TRUE)]
  vech.A <- full[stack]

  return(vech.A)
}
