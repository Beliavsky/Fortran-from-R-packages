#' Incomplete gamma function and asymptotic expansion
#'
#' \code{ecld.gamma} is the wrapper for incomplete gamma function
#' \eqn{\Gamma(s,x)}. It is mainly to wrap around \code{pgamma}.
#' And \code{ecld.gamma_hgeo} is the asymptotic expansion of \eqn{\Gamma(s,x)}
#' using hypergeometric series, \eqn{e^{-x} x^{s-1} {}_2 F_0 (1,1-s;;-1/x)}.
#' It is mainly used in for star OGF \eqn{L^{*}(k;\lambda)}.
#' \code{ecld.gamma_2F0} is simply \eqn{{}_2 F_0 (1,1-s;;-1/x)}, which is used
#' in the star OGF expansion.
#'
#' @param s numeric vector, for the order of incomplete gamma function
#' @param x numeric or MPFR vector
#' @param order numeric, the order of the power series
#' @param na.stop logical, stop if NaN is generated. The default is \code{TRUE}.
#'
#' @return numeric
#'
#' @keywords gamma
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.gamma
#' @export ecld.gamma_hgeo
#' @export ecld.gamma_2F0
#'
#' @importFrom stats pgamma
#'
### <======================================================================>
"ecld.gamma" <- function(s, x=0, na.stop=TRUE) {
    # TODO pgamma can't handle MPFR
    # It also has a limitation on how large s can be!!!
    pg_lim = 100
    s = ecd.mp2f(s) 
    
    N = length(s)
    M = length(x)
    if (N>1 & M != N & M != 1) {
        stop("ERROR: ecld.gamma reqires length(s)=length(x) when s is a vector")
    }

    # s can't have NA
    if (any(is.na(s))) stop(paste("ecld.gamma: s must be valid numeric, instead I get", s))
    
    # for very large s, recurrence relation
    if (any(s > pg_lim)) {
        # pgamma has a limit of s=100, so we have to use recurrence relation
        # this is somewhat slow though
        if (N==1) {
            G <- ecld.gamma(pg_lim, x, na.stop=na.stop)
            for (s1 in pg_lim:(s-1)) {
                G <- x^s1*exp(-x) + s1*G
            }
            return(G)
        } else if (N>1) { # s is vector
            if (M==1) x <- rep(x,N) # expand x if it is a number
            G <- x*0
            for (i in 1:N) {
                G[i] <- ecld.gamma(s[i], x[i])
            }
            return(G)
            
        } else {
            stop("Length of s is unknown")
        }
    }
    
    # pgamma, up to pg_lim
    stopifnot(all(s <= pg_lim))
    y <- pgamma(ecd.mp2f(x), ecd.mp2f(s), lower.tail = FALSE)*gamma(s)
    if (na.stop & anyNA(y)) {
        stop(paste("ERROR: ecld.gamma yields NaN for s=", s))
    }
    x*0+y # x*0 will preserve MPFR class, from x
}
### <---------------------------------------------------------------------->
#' @rdname ecld.gamma
"ecld.gamma_hgeo" <- function(s, x, order) {
    p <- exp(-x) * x^(s-1)
    q <- ecld.gamma_2F0(s, x, order)
    return(p*q)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.gamma
"ecld.gamma_2F0" <- function(s, x, order) {
    if (order == 0) return(1)
    stopifnot(order > 0)
    
    q <- 1
    xi <- 1
    si <- 1
    for (i in 1:order) {
        si <- si * (s-i)
        xi <- xi * x
        q <- q + si/xi
    }
    return(q)
}
### <---------------------------------------------------------------------->

