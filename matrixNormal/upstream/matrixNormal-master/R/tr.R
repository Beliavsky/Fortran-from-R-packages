#' Matrix Trace
#'
#' @family matrix
#' @keywords matrix
#'
#' @description Computes the trace of a square numeric matrix \emph{A}.
#'
#' @note If the argument is not a square numeric matrix, the function presents an error and
#'      terminates.
#'
#' @param A Square matrix.
#'
#' @examples
#' A <- matrix(seq(1, 16, 1), nrow = 4, byrow = TRUE)
#' A
#' tr(A)
#' tr(I(3))
#' @export

tr <- function(A) {
  # if( !is.numeric(A) & !is.matrix(A) )
  #      stop(paste( "A", "must be a numeric matrix."))
  if (nrow(A) != ncol(A)) {
    stop(paste("A", "is not a square matrix"))
  }
  stopifnot(nrow(A) == ncol(A))
  sum(diag(A))
}
