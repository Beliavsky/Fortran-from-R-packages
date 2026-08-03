#' Computes the theoretical statistics per state
#' 
#' This utility computes the statistics (mean, sd, kurtosis) based on the lambda distribution.
#' This is used to compare to the statistics from observations for each state.
#' \code{alloc} is a short-hand for Merton's optimal allocation, \code{mean/sd^2}.
#'
#' @param object an ldhmm object
#' @param annualize logical, to annaulize the sd and mean to V (xsqrt(days.pa)x100) and R (xdays.pa). 
#'                  Default is \code{FALSE}.
#' @param days.pa a positive integer, specifying number of days per year, default is 252.
#'
#' @return a matrix of statistics for each state, size of states times 4
#'
#' @keywords pdf
#'
#' @author Stephen H. Lihn
#'
#' @export 
#' 
### <======================================================================>
ldhmm.ld_stats <- function(object, 
                           annualize=FALSE, 
                           days.pa=252) {
    m <- object@m
    s <- matrix(NA, m, 4, byrow=TRUE)
    lds <- ldhmm.state_ld(object)
    for (j in 1:m) {
        s[j,1] <- ecld.mean(lds[[j]])
        s[j,2] <- ecld.sd(lds[[j]])
        s[j,3] <- ecld.kurtosis(lds[[j]])
        s[j,4] <- s[j,1]/s[j,2]^2
    }
    colnames(s) <- c("mean", "sd", "kurtosis", "alloc")
    
    if (annualize) {
        s[,1] <- s[,1] * days.pa
        s[,2] <- s[,2] * sqrt(days.pa) * 100
        s[j,4] <- s[j,1]/(s[j,2]/100)^2
        colnames(s) <- c("R", "V", "kurtosis", "alloc")
    }
    
    return(s)
}
### <---------------------------------------------------------------------->
