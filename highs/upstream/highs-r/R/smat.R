# Test if Object is a Sparse Matrix
#
# This function tests if an object is a sparse matrix of class "dgCMatrix".
# The dgCMatrix class represents sparse numeric matrices in compressed sparse column format.
#
# @param x An R object to be tested.
#
# @return A logical value: TRUE if x is a dgCMatrix, FALSE otherwise.
is_dgc <- function(x) inherits(x, "dgCMatrix")


is_csc <- function(x) inherits(x, "matrix.csc")


is_dgr <- function(x) inherits(x, "dgRMatrix")


is_csr <- function(x) inherits(x, "matrix.csr")


# Test if Simple Triplet Matrix
#
# Test if the object is a \code{"simple_triplet_matrix"}.
#
# @param x an object to be tested.
is_stm <- function(x) inherits(x, "simple_triplet_matrix")


# Test if Dense Matrix
#
# Test if the object is a dense matrix.
#
# @param x an object to be tested.
is_dmat <- function(x) inherits(x, "matrix") && length(dim(x)) == 2L



# Convert Matrix into Compressed Sparse Column Format
#
# Convert a matrix into compressed sparse column (CSC) format.
#
# @param x an object to be converted into the csc format.
#
# @references
# \url{https://www.gnu.org/software/gsl/doc/html/spmatrix.html}
as_csc <- function(x) {
    if (is_dgc(x)) {
        dim <- x@Dim
        list(col_ptr = x@p, row_id = x@i, value = x@x, nrow = dim[1L], ncol = dim[2L])
    } else if (is_stm(x)) {
        ind <- order(x$j, x$i)
        list(col_ptr = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
             row_id = x$i[ind] - 1L,
             value = x$v[ind],
             nrow = x[["nrow"]], ncol = x[["ncol"]])
    } else if (is_dmat(x)) {
        dim <- attributes(x)[["dim"]]
        ind <- which(x != 0, arr.ind = TRUE)
        list(col_ptr = c(0L, cumsum(tabulate(ind[, 2L], NCOL(x)))),
             row_id = ind[, 1] - 1L,
             value = x[ind],
             nrow = dim[1L], ncol = dim[2L])
    } else {
        stop(sprintf("unkown class %s", deparse(class(x))))
    }
}

