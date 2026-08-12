##
##  o c c u r s . R  Finding subsequences
##


count <- function(x, sorted = TRUE) {
    # two times faster than `table`, returns integers
    x = c(x); n = length(x)
    if (n == 0) {
        stop("Argument 'x' must not be empty.")
    } else if (n == 1) {
        return(list(v = x[1], e = c(1)))
    }
    if (sorted && is.unsorted(x)) x = sort(x)
    
    v = c(); e = c()
    x0 <- x[1]; e0 = 1
    for (i in 2:n) {
        if (x[i] == x0) {
            e0 = e0 + 1
        } else {
            v = c(v, x0); e = c(e, e0)
            x0 = x[i]; e0 = 1
        }
    }
    v = c(v, x0); e = c(e, e0)
    return(list(v = v, e = e))
}


occurs <- function(subseq, series){
    # Find all indices i such that series[i, ..., i+m-1] == subseq
    stopifnot(is.numeric(subseq), is.numeric(series))
    m <- length(subseq)
    n <- length(series)
    if (m > n)
        return(as.integer(c()))

    inds <- seq.int(length.out = n-m+1)
    for (i in seq.int(length.out = m))
        inds <- inds[ subseq[i] == series[inds + i - 1] ]
    inds
}
