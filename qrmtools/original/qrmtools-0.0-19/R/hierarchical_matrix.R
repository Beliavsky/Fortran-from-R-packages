### Constructing hierarchical matrices #########################################

##' @title Constructing Hierarchical Matrices
##' @param x list of length 2 or 3 containing the homogeneous entry of the
##'        current block, the components belonging to the current block and,
##'        possibly, another (nested) list of the same type.
##'        See also 'nacList' as described on ?onacopulaL of the R package copula.
##' @param diagonal diagonal entries (default 1)
##' @return Hierarchical matrix
##' @author Marius Hofert
##' @note Based on ideas of copula::nacPairthetas
hierarchical_matrix <- function(x, diagonal = rep(1, d))
{
    stopifnot(is.list(x))
    ## Auxiliary function for determining d (size of return matrix M)
    len <- function(x) {
        if(is.list(x)) {
            l <- length(x)
            stopifnot(l == 2 || l == 3) # sanity check
            if(length(x) == 3) {
                slist <- if(is.list(x[[3]][[1]])) x[[3]] else x[3]
                length(x[[2]]) + sum(sapply(slist, len)) # current components + of all sub-lists
            } else sum(sapply(x[[2]], len)) # iterate over all sub-lists (in case there are multiple)
        } else {
            length(x)
        }
    }
    d <- len(x) # determine d
    ## Allocate the return matrix M
    M <- matrix(NA_real_, nrow = d, ncol = d)
    ## Fill M
    setM <- function(x) { # x = (nested) list
        ii <- x[[2]]
        if(length(ii) > 0) M[ii, ii] <<- x[[1]]
        if(length(x) > 2) { # ... there are sub-lists
            slist <- if(is.list(x[[3]][[1]])) x[[3]] else x[3] # to prevent unlist()ing errors in case there is only one sub-list
            for(lst in slist) { # recursion
                jj <- setM(lst) # determine the indices corresponding to the sub-list
                M[ii,jj] <<- M[jj,ii] <<- x[[1]] # set M (corresponding row and column)
                ii <- c(ii, jj) # append indices (vector of those indices for which we already determined M)
            }
        }
        ii # return (index of all (child) components for which we already determined M)
    }
    setM(x)
    ## Set diagonal of M and return
    diag(M) <- diagonal
    M
}
