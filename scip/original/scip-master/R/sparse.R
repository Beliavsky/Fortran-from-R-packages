#' Convert a matrix to CSC (Compressed Sparse Column) format
#'
#' @param x A matrix, dgCMatrix, or simple_triplet_matrix
#' @return A list with components \code{i} (row indices, 0-based),
#'   \code{p} (column pointers, 0-based), \code{x} (values),
#'   \code{nrow}, \code{ncol}
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @keywords internal
make_csc_matrix <- function(x) {
    if (inherits(x, "dgCMatrix")) {
        list(i = x@i, p = x@p, x = x@x, nrow = x@Dim[1], ncol = x@Dim[2])
    } else if (inherits(x, "simple_triplet_matrix")) {
        m <- Matrix::sparseMatrix(i = x$i, j = x$j, x = x$v,
                                  dims = c(x$nrow, x$ncol),
                                  repr = "C")
        list(i = m@i, p = m@p, x = m@x, nrow = m@Dim[1], ncol = m@Dim[2])
    } else if (is.matrix(x)) {
        m <- methods::as(x, "dgCMatrix")
        list(i = m@i, p = m@p, x = m@x, nrow = m@Dim[1], ncol = m@Dim[2])
    } else if (inherits(x, "sparseMatrix")) {
        m <- methods::as(x, "dgCMatrix")
        list(i = m@i, p = m@p, x = m@x, nrow = m@Dim[1], ncol = m@Dim[2])
    } else {
        stop("Unsupported matrix type: ", class(x)[1])
    }
}
