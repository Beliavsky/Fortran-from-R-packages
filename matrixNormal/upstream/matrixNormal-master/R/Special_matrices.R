#' Generating Special Matrices
#'
#' @family matrix
#' @keywords matrix identity ones
#'
#' @name special.matrix
#'
#' @description Creates an Identity Matrix \emph{I} and a Matrix of Ones \emph{J}. \itemize{
#' \item \code{I}(): Creates an identity matrix where the number of columns is n.  This is a diagonal matrix with all equal to one (1). An identity matrix is usually written as \emph{I}. Names of rows and columns (\code{dimnames}) are included.
#' \item \code{J}(): Creates a matrix of ones with any number of rows and columns. Names of rows and columns (\code{dimnames}) are included.
#' }
#'
#' @param n Number of rows in \emph{I} or \emph{J}.
#' @param m Number of columns in \emph{J}. Default: Same as number of rows.
#'
#' @examples
#' # To create an identity matrix of order 12
#' I(2)
#' # To make a matrix of 6 rows and 10 columns of all ones
#' J(6, 10)
#' # To make a matrix of unity, dimensions 6 x 6.
#' J(6)
#' @rdname  special.matrix
#' @export I
I <- function(n) {
  # Identity Matrix where number of columns is n.
  D <- diag(1, n)
  dimnames(D) <- list(1:n, 1:n)
  return(D)
}

#' @rdname  special.matrix
#' @export J
J <- function(n, m = n) {
  # A matrix of ones with n rows and m columns.
  matrix(c(1), nrow = n, ncol = m, dimnames = list(1:n, 1:m))
}
