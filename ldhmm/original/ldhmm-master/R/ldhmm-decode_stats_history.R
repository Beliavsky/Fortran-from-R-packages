#' Estimating historical statistics (mean, volatility and kurtosis)
#' 
#' This utility estimates historical statistics (mean, volatility and kurtosis) according to the state probabilities.
#' The ldhmm object must have been decoded by running through \code{ldhmm.decoding} function.
#' Note that kurtosis is naively implemented as the linear sum from each state weighted by state probabilities.
#' It is subject to change to more rigorous formula in future releases.
#' 
#' @param object a decoded ldhmm object
#' @param ma.order a positive integer or zero, specifying order of moving average. Default is zero.
#' @param annualize logical, to annaulize the sd and mean to V (xsqrt(days.pa)x100) and R (xdays.pa). 
#'                  Default is \code{FALSE}.
#' @param days.pa a positive integer, specifying number of days per year, default is 252.
#' 
#' @return an matrix of statistics history, size of observations times size of 3
#'
#' @keywords forecast
#'
#' @author Stephen H. Lihn
#'
#' @export 
#' 
### <======================================================================>
ldhmm.decode_stats_history <- function(object, ma.order=0, 
                                       annualize=FALSE, 
                                       days.pa=252) 
{
    if (length(ma.order) != 1) stop("ma.order must an integer")
    if (ma.order < 0) stop("ma.order must be positive")
    n <- length(object@observations)
    ld_stats <- ldhmm.ld_stats(object)
    
    calc_st <- function(j, si) {
        prob <- object@states.prob[,j]
        st <- ld_stats[,si]
        if (si==2) {
            mn <- calc_st(j,1)
            v2 <- (ld_stats[,1]-mn)^2+ld_stats[,2]^2
            var <- as.numeric(prob %*% v2)
            return (sqrt(var))
        } 
        else as.numeric(prob %*% st)
    }
    stats_with_ma <- function(si) {
        v <- sapply(1:n, function(j) calc_st(j,si))
        ldhmm.sma(v, order=ma.order)
    }
    rs <- data.frame(mean = stats_with_ma(1), 
                     sd = stats_with_ma(2), 
                     kurtosis = stats_with_ma(3))
    if (annualize) {
        rs$mean <- rs$mean * days.pa
        rs$sd <- rs$sd * sqrt(days.pa) * 100
        colnames(rs) <- c("R", "V", "kurtosis")
    }
    return (as.matrix(rs))
}
### <---------------------------------------------------------------------->

