##
##  b i n p a c k i n g . r  Approximate Bin Packing
##


## =====================================================================
## first-fit and best-fit heuristic for bin packing

# Approximate solution to the bin packing problem
# Method: first_fit, input item size, and bin size
# Returns: no.of bins, assigned bin, bin allocation
# First-Fit, Best-Fit, decreasing = TRUE
bpp_approx <- function(S, cap,
                       method = c("firstfit", "bestfit", "worstfit")) {
    stopifnot(is.numeric(S), is.numeric(cap))
    method <- match.arg(method)
    if (length(cap) != 1 || any(S > cap))
        stop("Argument 'cap' must be a scalar and all S <= cap.")
    if (is.unsorted(rev(S))) {
        warning("Argument 'S' of item weights is not decreasingly sorted.",
                call. = FALSE, immediate. = TRUE)
    }

    n <- length(S)
    b <- rep(cap, n)
    p <- rep(0, n)
    
    if (method == "firstfit") {
        for (i in 1:n) {
            j  <- which(b >= S[i])[1]
            p[i] <- j
            b[j] <- b[j] - S[i]
        }
    } else if (method == "bestfit") {
        for (i in 1:n) {
            inds <- which(b >= S[i])
            j <- inds[which.min(b[inds])]
            p[i] <- j
            b[j] <- b[j] - S[i]
        }
    } else if (method == "worstfit") {
        for (i in 1:n) {
            inds <- which(b >= S[i] & b < cap)
            if (length(inds) == 0) inds <- which(b >= S[i])
            j <- inds[which.max(b[inds])]
            p[i] <- j
            b[j] <- b[j] - S[i]
        }
    } else
        stop("Unknown method; may be implemented in the course of time.")
    
    m <- max(p)
    f <- sum(S) / m / cap
    return(list(nbins = m, xbins = p, sbins = cap - b[1:m], filled = f))
}
