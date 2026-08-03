#' Simulating auto-correlation (ACF)
#' 
#' This utility simulates the auto-correlation. The first few lag of ACF should match 
#' the ACF from the market data fairly well. This is a major validation of a successful HMM.
#' Be aware this is a CPU intensive calculation. It uses the multi-core functionality.
#'
#' @param object an ldhmm object that can supply m, param.nbr and stationary.
#' @param n a positive integer specifying number of observations to simulate.
#' @param lag.max a positive integer, specifying number of lags to be computed.
#' @param debug logical, specifying to print progress message or not. Default is \code{FALSE}. 
#'
#' @return a vector of ACF
#'
#' @keywords acf
#'
#' @author Stephen H. Lihn
#' 
#' @importFrom stats cor
#'
#' @export 
#'
### <======================================================================>
ldhmm.simulate_abs_acf <- function(object, n=10000, lag.max=5, debug=FALSE) {
    ac <- c()
    h2 <- h1 <- ldhmm.simulate_state_transition(object, init=n)
    if (debug) print(paste(Sys.time(), "Finished init"))
    for (k in 1:lag.max) {
        h2 <- ldhmm.simulate_state_transition(h2)
        ac[k] <- stats::cor(abs(h1@observations), abs(h2@observations))
        if (debug) print(paste(Sys.time(), "Finished lag", k))
    }
    return(ac)
}
### <---------------------------------------------------------------------->

