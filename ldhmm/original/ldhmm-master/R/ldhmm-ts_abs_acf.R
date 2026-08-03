#' Computing ACF of the absolute value of a time series
#' 
#' This utility computes the ACF of the absolute value of a time series as a proxy 
#' of the auto-correlation of the volatility. It allows to drop the largest N outliers
#' so that they would not skew the ACF calculation.
#'
#' @param x numeric, the observations.
#' @param drop a positive integer, specifying number of outliers to be dropped.
#' @param lag.max a positive integer, specifying number of lags to be computed.
#'
#' @return a vector of ACF
#'
#' @keywords acf
#'
#' @author Stephen H. Lihn
#'
#' @importFrom stats acf
#' 
#' @export
#' 
### <======================================================================>
ldhmm.ts_abs_acf <- function (x, drop=0, lag.max=100) {
    x <- as.numeric(abs(x))
    y <- rev(x[order(x)])
    `%notin%` <- function (x, table) match(x, table, nomatch = 0L) == 0L
    x1 <- x[x %notin% head(y,drop)] # acf can be skewed by largest movements
    a <- stats::acf(abs(x1), lag.max=lag.max, plot=FALSE)
    a$acf[-1]
}
### <---------------------------------------------------------------------->
