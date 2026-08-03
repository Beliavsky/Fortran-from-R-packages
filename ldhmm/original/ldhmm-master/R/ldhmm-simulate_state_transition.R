#' Simulating state transition
#' 
#' This utility allows to simulate the states and obervations over time. 
#' Be aware this is a CPU intensive calculation. It uses the multi-core functionality.
#'
#' @param object an ldhmm object that can supply m, param.nbr and stationary.
#' @param init a positive integer specifying number of observations to simulate initially.
#'             The default is NULL, indicating that the simulation should use the (local) states and observations
#'             from within the object, and simulate the next set of random states and observations according to gamma.
#'             When init is an integer, the utility will generate random states and observations according to delta.   
#'
#' @return an ldhmm object containing the simulated states and observations.
#'         The observations are stored in the \code{observations} slot.
#'         The states are stored in the \code{states.local} slot.
#'
#' @keywords simulation
#'
#' @author Stephen H. Lihn
#'
#' @importFrom stats runif
#' @importFrom stats uniroot
#'
#' @export 
#'
### <======================================================================>
ldhmm.simulate_state_transition <- function(object, init=NULL) {
    m <- object@m 
    
    lds <- ldhmm.state_ld(object)
    sim_one_X <- function(i) {
        c <- stats::runif(1)
        x <- stats::uniroot(function(x) ecld.cdf(lds[[i]],x)-c, lower=-1, upper=1)$root
        return(x)
    }
    
    app <- function(x,f) simplify2array(parallel::mclapply(x,f))
    
    # generate initial states
    if (! is.null(init)) {
        cum_delta <- cumsum(object@delta)
        init_state <- function(x) length(which(x > cum_delta))+1
        h <- object
        h@states.local <- sapply(runif(init), init_state)
        h@observations <- app(h@states.local, sim_one_X)
        return (h)
    }
    
    # simulate subseqent states
    cum_prob <- matrix(NA, nrow=m, ncol=m, byrow=TRUE)
    for (i in 1:m) cum_prob[i,] <- cumsum(object@gamma[i,])
    
    next_state <- function(i) {
        x <- runif(1)
        length(which(x > cum_prob[i,]))+1
    }
    h <- object
    h@states.local <- sapply(object@states.local, next_state)
    h@observations <- app(h@states.local, sim_one_X)
    return (h)
}
### <---------------------------------------------------------------------->

