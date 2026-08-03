#' Computing the volatility forecast for next one period
#' 
#' This utility computes the volatility forecast based on the given future observations for next one period.
#'
#' @param object an ldhmm object
#' @param x numeric, the observations.
#' @param xf numeric, the future observations to be forecasted.
#' @param ma.order a positive integer or zero, specifying order of moving average. Default is zero.
#' @param days.pa a positive integer specifying trading days per year, default is 252.
#'
#' @return matrix of future observations and volatilities, size of 2 times length of \code{xf}.
#'
#' @keywords forecast
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
#' @importFrom utils tail
#' 
### <======================================================================>
ldhmm.forecast_volatility <- function(object, x, xf, ma.order=0, days.pa=252)
{
    n <- length(x)
    nxf <- length(xf)
    df <- matrix(0, nrow=2, ncol=nxf, byrow=TRUE)
    for (i in 1:nxf) {
        hd <- ldhmm.decoding(object, c(as.numeric(x),xf[i]))
        vd <- ldhmm.decode_stats_history(hd, ma.order=ma.order)[,2]
        fv <- utils::tail(vd,1)*sqrt(days.pa)*100
        df[1,i] <- xf[i]
        df[2,i] <- fv
    }
    df
}
### <---------------------------------------------------------------------->

