#' Normalize constraint matrix
#'
#' @description it is not uncommon for quadprog to fail when there are large differences in 2-norm
#' between the columns of the constraint matrix (Amat).  It is possible to alleviate this issue in
#' some cases by normalizing the constraints (and their boundaries, defined by bvec).
#'
#' @param Amat constraint matrix as defined by solve.QP
#' @param bvec constraints as defined by solve.QP
#'
#' @seealso quadprog::solve.QP
#' @seealso quadprog::solve.QP.compact
#' @return a list with two elements: Amat and bvec that contain the normalized constraints.
normalizeConstraints <- function(Amat, bvec){

    norm2 <- sqrt(colSums(Amat ^ 2))
    if(any(norm2 == 0)){
        stop("At least one column of Amat has a zero 2-norm")
    }

    Anorm <- t.default(t.default(Amat) / norm2)
    
    list(
        Amat = Anorm,
        bvec = bvec / norm2
    )    
}
