##
##  a s s i g n m e n t . R  Linear Assignment Problem
##


assignment <- function(cmat, dir = "min") {
    stopifnot(is.numeric(cmat), is.matrix(cmat))
    if (! dir %in% c("min", "max"))
        stop("Argument 'dir' must be 'min' or 'max'.")
    n <- nrow(cmat); m <- ncol(cmat)
    if (n != m || n <= 1) 
        stop("Argument 'cmat' must be a quadratic, numeric matrix.")
    
    smat <- lp.assign(cmat, direction = dir)$solution
    inds <- which(round(smat) == 1, arr.ind = TRUE)
    
    ordr <- order(inds[, 1])
    perm <- inds[, 2][ordr]
    mini <- sum(cmat * smat)
    
    return(list(perm = perm, min = mini, err = 0))
}


# assignment <- function(cmat) {
#     stopifnot(is.numeric(cmat), is.matrix(cmat))
#     n <- nrow(cmat); m <- ncol(cmat)
#     if (n != m || n <= 1)
#         stop("Argument 'cmat' must be a quadratic, numeric matrix.")
#     if (any(floor(cmat) != ceiling(cmat))) {
#         warning("Matrix 'cmat' not integer; will take floor of it.")
#         cmat <- floor(cmat)
#     }
# 
#     a <- cbind(t(cmat), rep(0, n))
#     b <- rep(0, n)
#     t <- 0; ierr <- 0
#     # dummies
#     iwk <- vector("integer", 7*n+2)
# 
#     S <- .Fortran("assgn", as.integer(n), as.integer(a),
#                            b = as.integer(b), t = as.integer(t),
#                            as.integer(iwk), err = as.integer(ierr),
#                            PACKAGE = "adagio")
# 
#     if (S$err != 0)
#         warning("Integer overflow happened; result may not be correct.")
#     return(list(perm = S$b, min = S$t, err = S$err))
# }
