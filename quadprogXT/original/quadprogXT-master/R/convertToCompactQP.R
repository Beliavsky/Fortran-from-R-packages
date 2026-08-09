#' Make a Constraint Matrix Compact
#'
#' @param Amat A constraint matrix as defined in solve.QP.
#' @return A list containing the two elements to be passed to solve.QP.compact,
#' each named accordingly.
#' @seealso quadprog::solve.QP
#' @seealso quadprog::solve.QP.compact
#' @export
convertToCompact <- function(Amat){

  ## Get the dimensions of the input and test what is non-zero.
  d <- dim(Amat)
  Az <- Amat != 0

  ## Count the non-zero entries in each column and note the maximum. This will
  ## be the number of rows in the output matrices (minus one for Aind).
  cz <- colSums(Az)
  mc <- max(cz)

  ## Throw an error if there are columns of all zeros.
  if(any(cz == 0)) stop("Some columns of the constraint matrix are all zero.")

  ## Build the output matrices.
  Am <- matrix(data = 0, ncol = d[2L], nrow = mc)
  Ai <- matrix(data = 0, ncol = d[2L], nrow = mc + 1L)

  ## Note the indices for which the input matrix is non-zero as well as those
  ## for which the output matrices will be non-zero. The latter is determined
  ## by testing whether the row in the output matrix is less than the
  ## count of non-zero entries in the corresponding column of A given by cz.
  ii <- which(Az)
  io <- which(row(Am) <= cz[col(Am)])

  ## Write the non-zero entries of the input to Amat, and the non-zero counts
  ## to the first row of Aind.
  Am[io] <- Amat[ii]
  Ai[1L,] <- cz

  ## The entries in Aind (past the first row) correspond to the row index of
  ## the non-zero elements in the input matrix. The expression below is a bit
  ## nasty, however it's equivalent to row(A)[ii], just much faster when A is
  ## sparse.
  Ai[-1L, ][io] <- (ii - 1L) %% d[1L] + 1L

  ## Write the matrices to the output list and return.
  out <- list(Amat = Am, Aind = Ai)
  return(out)
}
