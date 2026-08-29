#' Stacks a Matrix using matrix operator "vec"
#'
#' @family matrix
#' @keywords distribution matrix
#'
#' @description Returns a column vector that stacks the columns of \emph{A}, a m x n matrix.
#'
#' @note
#' \enumerate{
#' \item Unlike other `vec()` functions on CRAN, matrixNormal version inherits names from matrices to their vectorized forms.
#' \item vec() was adapted from Frederick Novomestky's \pkg{matrixcalc}. This function is edited so that it can take dimension names and return the matrix as a vector.
#' \item These functions were used as accessories for other matrixNormal functions.
#' }
#'
#' @references
#' Magnus, J. R. and H. Neudecker (1999). \emph{Matrix Differential Calculus with Applications in Statistics and Econometrics.} Second Edition, John Wiley, ed.
#'
#' @param A A matrix with m rows and n columns.
#' @param use.Names Logical. If TRUE, the names of A are taken to be names of the stacked matrix. Default: TRUE.
#' @return A vector with mn elements.
#'
#' @examples
#' M <- matrix(c(4, 5, 6, 7, 8, 9), nrow = 3)
#' M
#' vec(M)
#' if (!requireNamespace("matrixcalc", quietly = TRUE)) {
#' # Compare vec from \pkg{matrixcalc} and new function.
#' matrixcalc::vec(M)
#' # The names are rownames(M):colnames(M) in that order.
#' # Very similar to matrixcalc but dimension names are different.
#' } else {
#'   message("you need to install the package matrixcalc to compare this example.")
#' }
#' @export

vec <- function(A, use.Names = TRUE) {
  if (is.vector(A)) {
    return(A)
  }
  if (!is.matrix(A)) stop("argument A is not a matrix")
  if (!is.numeric(A)) stop("argument A is not a numeric matrix")
  vec.A <- as.vector(A)

  ## PAUL HARGARTEN ADDED THIS CODE:
  if (use.Names) {
    n <- nrow(A)
    c <- ncol(A)

    # if no dimnames, add # instead.
    if (is.null(dimnames(A)[[1]])) {
      dimnames(A)[[1]] <- 1:n
    }
    if (is.null(dimnames(A)[[2]])) {
      dimnames(A)[[2]] <- 1:c
    }
    names(vec.A) <- paste(rep(dimnames(A)[[1]], c),
      rep(dimnames(A)[[2]], each = n),
      sep = ":"
    )
  }
  return(vec.A)
}
