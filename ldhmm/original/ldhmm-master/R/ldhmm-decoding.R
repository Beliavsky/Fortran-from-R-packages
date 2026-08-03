#' Computing the minus log-likelihood (MLLK)
#' 
#' This utility computes the state probabilities, uses local and global decoding to calculate the states.
#' The results are saved to the returned \code{ldhmm} object.
#'
#' @param object       an ldhmm object
#' @param x            numeric, the observations.
#' @param do.global    logical, if \code{TRUE} (default), perform Viterbi decoding.
#' @param do.stats     logical, if \code{TRUE} (default), calculate stats.
#'
#' @return an ldhmm object containing results of decoding
#'
#' @keywords mllk
#'
#' @author Stephen H. Lihn
#'
#' @export ldhmm.decoding
#' 
### <======================================================================>
ldhmm.decoding <- function(object, x, do.global=TRUE, do.stats=TRUE)
{
    m <- object@m
    x <- as.numeric(x)
    n <- length(x)
    state_probs <- matrix(NA, nrow=m, ncol=n)

    la <- ldhmm.log_forward(object, x)
    lb <- ldhmm.log_backward(object, x)
    c  <- max(la[,n])
    llk <- c + log(sum(exp(la[,n]-c)))
    for (i in 1:n) state_probs[,i] <- exp(la[,i]+lb[,i]-llk)

    # local decoding for most likely states
    get_state <- function(i) {
        p <- c(na.exclude(state_probs[,i]))
        if (length(p) == 0) return(1) # fall back to the first state (check???)
        min(which(p==max(p))) # min to ensure scalar
    }
    local_states <- sapply(1:NCOL(state_probs), get_state)
    
    object@observations <- x
    object@states.prob <- state_probs
    object@states.local <- local_states
    if (do.global) object@states.global <- ldhmm.viterbi(object, x)
    
    # stats
    if (do.stats) {
        object@states.local.stats <- ldhmm.calc_stats_from_obs(object, use.local=TRUE)
        if (do.global) object@states.global.stats <- ldhmm.calc_stats_from_obs(object, use.local=FALSE)
    }
    return(object)
}
### <---------------------------------------------------------------------->
