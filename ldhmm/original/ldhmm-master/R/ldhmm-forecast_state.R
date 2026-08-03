#' Computing the state forecast
#' 
#' This utility computes the state forecast, given the sequence of observations in the past.
#'
#' @param object an ldhmm object
#' @param x numeric, the observations.
#' @param h integer, time steps to forecast.
#'
#' @return matrix of probabilities per state (even if h=1), number of states times size of h
#'
#' @keywords forecast
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
### <======================================================================>
ldhmm.forecast_state <- function(object, x, h=1)
{
    n <- length(x)
    la <- ldhmm.log_forward(object, x)
    c  <- max(la[,n])
    llk <- c + log(sum(exp(la[,n] - c)))
    state_preds <- matrix(NA, nrow=object@m, ncol=h)
    phi <- exp(la[,n] - llk)
    for (i in 1:h) {
        phi <- phi %*% object@gamma
        state_preds[,i] <- phi
    }
    return(state_preds)
}
### <---------------------------------------------------------------------->

