##
## copied from the Rsymphony package from
## Kurt Hornik and Stefan Theussl and Reinhard Harter
##
## Simple functions for converting "matrix" type objects into the
## sparse "column major order" (CSC, modulo offsets) format used by
## SYMPHONY.

## matind: vector of the row indices corresponding to each entry of
##   value
## values: vector of the values of nonzero entries of the constraint
##   matrix in column order.

#' Convert a plain matrix or simple triplet form matrix to a [Matrix::dgCMatrix-class] (implicit) form
#' @param x a matrix or a simple triplet form matrix
#' @return a list of row pointer, column pointer, and values corresponding to a [Matrix::dgCMatrix-class] object
make_csc_matrix <- function(x) UseMethod("make_csc_matrix")

#' @method make_csc_matrix default
make_csc_matrix.default <- function(x) stop("x must be a matrix/simple_triple_matrix/Matrix!")

#' @method make_csc_matrix matrix
make_csc_matrix.matrix <- function(x) {
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")

    ind <- which(x != 0, arr.ind = TRUE)
    list(matbeg = c(0L, cumsum(tabulate(ind[, 2L], ncol(x)))),
         matind = ind[, 1] - 1L,
         values = x[ind])
}

#' @method make_csc_matrix simple_triplet_matrix
make_csc_matrix.simple_triplet_matrix <- function(x) {
    if(!inherits(x, "simple_triplet_matrix"))
        stop("Argument 'x' must be of class 'simple_triplet_matrix'.")

    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
}

#' @method make_csc_matrix Matrix
#' @import Matrix
make_csc_matrix.Matrix <- function(x) {
  if (!inherits(x, "dgCMatrix")) {
    if (is(x, "ddiMatrix")) {
      indices <- seq_len(nrow(x))
      list(matbeg = c(0L, indices),
           matind = indices - 1L,
           values = x@x)
    } else {
      tmp <- as(as(x, "CsparseMatrix"), "generalMatrix", "dMatrix")
      list(matbeg = tmp@p, matind = tmp@i, values = tmp@x)
    }
  } else {
    list(matbeg = x@p, matind = x@i, values = x@x)
  }
}

## Added by @bnaras for sparse symmetric matrix.

#' Convert a plain matrix or simple triplet form matrix to a [Matrix::dsCMatrix-class] (implicit, upper) form
#' @param m a matrix or a simple triplet form matrix
#' @return a list of row pointer, column pointer, and values corresponding to a [Matrix::dsCMatrix-class] object
make_csc_symm_matrix <- function(m) UseMethod("make_csc_symm_matrix")

#' @method make_csc_symm_matrix matrix
make_csc_symm_matrix.matrix  <- function(m) {
  ind <- which(m!=0 ,arr.ind = TRUE)
  ind  <- ind[ind[, 1] <= ind[, 2], ] ## keep upper part only
  values  <- m[ind]
  x  <- list(i = ind[, 1], j = ind[, 2])  ## triplet form
  ind  <- order(x$j, x$i)  ## may not be needed
  list(matbeg = c(0L, cumsum(tabulate(x$j[ind], ncol(m)))),
       matind = x$i[ind] - 1L,
       values = values)
}

#' @method make_csc_symm_matrix Matrix
#' @import Matrix
make_csc_symm_matrix.simple_triplet_matrix  <- function(m) {
  ind <- which(m$i <= m$j)
  x <- list(i = m$i[ind] + 1L, j = m$j[ind] + 1L) ##make it 1-based
  values <- m$v[ind]
  ind  <- order(x$j, x$i)  ## may not be needed
  list(matbeg = c(0L, cumsum(tabulate(x$j[ind], m$ncol))),
       matind = x$i[ind] - 1L,
       values = values)
}

## Ensure that a matrix is of class "dgCMatrix"
## Note: in Matrix, a sparseMatrix need not be dgCMatrix;
## it could be ddiMatrix, dgTMatrix, etc. So coerce!
#' @importFrom methods as is
ensure_dgc_matrix <- function(mat, nrow, ncol) {
  if (inherits(mat, "simple_triplet_matrix")) {
    csc <- make_csc_matrix(mat)
    Matrix::sparseMatrix(i = csc[["matind"]], p = csc[["matbeg"]], x = csc[["values"]],
                         dims = c(nrow, ncol), index1 = FALSE)
  } else {
    if (!inherits(mat, "dgCMatrix")) {
      as(as(mat, "generalMatrix"), "CsparseMatrix")
    } else {
      mat
    }
  }
}
    
