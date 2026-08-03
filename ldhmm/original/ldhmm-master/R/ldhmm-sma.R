#' Simple moving average of a time series
#' 
#' This utility calculates simple moving average, with option to backfill for NA.
#'
#' @param x numeric, the time series.
#' @param order a positive integer to specify order of moving average. 
#' @param na.backfill logical, specify whether to backfill for NA. Default is \code{TRUE}.
#' 
#' @return numeric, simple moving average, same length as x.
#'
#' @keywords data
#'
#' @author Stephen H. Lihn
#'
#' @export 
#' 
#' @examples
#' x <- 1:100
#' a <- ldhmm.sma(x, 10)
### <======================================================================>
ldhmm.sma <- function(x, order, na.backfill=TRUE)
{
    if (length(order) != 1) stop("order must be an integer")
    if (order == 0) return (x)
    if (order < 0) stop("order must be positive")
    if (! is(x, "numeric")) x <- as.numeric(x)
    
    ma <- x*NA
    ma[1] <- x[1]
    for (i in 2:length(x)) {
        if (na.backfill) {
            if (is.na(x[i])) x[i] <- x[i-1] # backfill NA
        }
        j <- i-order+1
        if (j <= 0) j <- 1
        ma[i] <- mean(x[j:i], na.rm=TRUE)
    }
    return(ma)
}
