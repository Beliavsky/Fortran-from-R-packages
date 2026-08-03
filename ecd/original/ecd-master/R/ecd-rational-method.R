#' Utility to convert a numeric to a rational
#' 
#' Convert a numeric x to rational p/q, which is then used for polynomial construction.
#' It can be used for displaying the time as fraction of a year too.
#'
#' @param x numeric
#' @param cycles numeric, maximum number of steps, default is 10.
#' @param pref.denominator numeric, a list of preferred integer denominators to conform to, default is \code{numeric(0)}.
#' @param max.denominator numeric, maximum denominator when the loop of trial should stop, default is 500.
#' @param as.character logical, if specified, convert to character of p/q, default is \code{FALSE}.
#'
#' @return vector of two integers, representing numerator and denominator. 
#'         If as.character is true, then return character instead of the rational pair.
#'         If x is a vector and as.character is false, return a matrix of length(x) by 2.
#'
#' @keywords solve
#'
#' @export 
#'
#' @examples
#' pq1 <- ecd.rational(2.5)
#' pq2 <- ecd.rational(1/250)
### <======================================================================>
ecd.rational <- function(x, pref.denominator=numeric(0), cycles=10, max.denominator=500, as.character=FALSE)
{
    if (length(x) > 1) {
        f <- function(x) ecd.rational(x, cycles=cycles, 
                                      pref.denominator=pref.denominator,
                                      max.denominator=max.denominator, 
                                      as.character=as.character)
        if (!as.character) return(t(sapply(x,f)))
        else return(sapply(x,f))
    }
    
    if (length(x) != 1 | !is.finite(x)) {
        stop("Input must be length-one finite numeric!")
    }
    
    a0 <- rep(0, length(x))
    b0 <- rep(1, length(x))
    A <- matrix(b0)
    B <- matrix(floor(x))
    r <- as.vector(x) - drop(B)
    
    len <- 0
    while(any(r > 1/max.denominator) && (len <- len + 1) <= cycles) {
        a <- a0
        b <- b0
        a <- 1
        r <- 1/r
        b <- floor(r)
        r <- r - b
        A <- cbind(A, a)
        B <- cbind(B, b)
    }
    
    pq1 <- cbind(b0, a0)
    pq <- cbind(B[, 1], b0)
    len <- 1
    while((len <- len + 1) <= ncol(B)) {
        pq0 <- pq1
        pq1 <- pq
        pq <- B[, len] * pq1 + A[, len] * pq0
    }
    
    R <- unname(as.vector(pq))
    
    for (pd in pref.denominator) {
        if (floor(pd/R[2]) == pd/R[2]) {
            R <- R*(pd/R[2])            
            break
        }
    }
    
    if (!as.character) return(R) 
    else {
        if (R[1] == 0) return("0")
        if (R[1] == R[2]) return("1")
        return(sprintf("%d/%d", R[1], R[2]))
    }
}
