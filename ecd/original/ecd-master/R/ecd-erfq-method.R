#' Quartic scaled error function
#'
#' The scaled error function in quartic pricing model that encaptulates
#' both scaled \code{erfi} and \code{erfc} functions
#' into a single representation. This is used to provide an elegant expression for
#' the MGF and local option prices, \eqn{L_{c,p}}.
#' When \code{sgn=-1}, it is \eqn{\sqrt{\pi}e^{-x^2} erfi(x)}, which twice of Dawson function.
#' When \code{sgn=1}, it is \eqn{\sqrt{\pi}e^{x^2} erfc(x)}.
#' \code{ecd.erfq_sum} is the summation implementation with truncation rule set forth
#' in the quartic pricing model. It achieves high precision when x > 4.5.
#'
#'
#' @param x numeric
#' @param sgn an integer of 1 or -1
#'
#' @return The \code{mpfr} object
#'
#' @keywords utility
#'
#' @export ecd.erfq
#' @export ecd.erfq_sum
#'
#' @examples
#' x <- ecd.erfq(c(5,10,15), 1)
#' y <- ecd.erfq(c(5,10,15), -1)
### <======================================================================>
ecd.erfq <- function(x, sgn) {
    if (sgn != 1 & sgn != -1) {
        stop("sgn must be 1 or -1")
    }
    
    x <- ecd.mp1*x # force MPFR type for precision
    c <- sqrt(ecd.mppi())
    y <- if (sgn == -1) ecd.dawson(x)*2 else c*ecd.erfcx(x)
    return(y)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.erfq
ecd.erfq_sum <- function(x, sgn) {
    if (sgn != 1 & sgn != -1) {
        stop("sgn must be 1 or -1")
    }

    c <- sqrt(ecd.mppi())
    y <- ecd.mp1*0*x
    for (i in 1:length(x)) {
        mx <- floor(ecd.mp2f(x[i]^2))
        if (mx > 25) mx <- 25
        m <- seq(0, mx, by=1)
        y[i] <- 1/c*sum(gamma(m+1/2)*(-sgn)^m/(x[i]^2)^m)
    }
    return(y/x)
}
### <---------------------------------------------------------------------->
