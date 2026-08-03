#' Utility to diff a vector of numeric or mpfr to get first derivative
#'
#' This utility uses diff to get first derivative dy/dx.
#' but it handles mpfr vector properly
#'
#' @param y a vector of numeric or mpfr
#' @param x a vector of numeric or mpfr
#' @param pad integer, to manage padding so that the output vector has
#'            the same length as the input. 0 for no padding,
#'            1 to repeat the first element, -1 to repeat the last element.
#'
#' @return the derivative vector
#'
#' @keywords utility
#'
#' @export
#'
#' @examples
#' d <- ecd.diff(c(10,20,30), c(1,2,3), pad=1)
#'
### <======================================================================>
ecd.diff <- function(y, x, pad=0)
{
    ret <- diff(y)/diff(x)
    
    stopifnot (pad %in% c(-1,0,1))
    if (pad==1) ret <- c(ret[1], ret)
    if (pad==-1) ret <- c(ret, ret[length(ret)-1])
    
    return(ret)
}
### <---------------------------------------------------------------------->
