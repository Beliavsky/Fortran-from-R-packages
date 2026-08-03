#' Transforming natural parameters to a linear working parameter array
#' 
#' This utility linearizes the natural parameters and transforms the contrained parameters
#' to unconstrained parameters. (Zucchini, 3.3.1)
#'
#' @param object an ldhmm object
#' @param mu.scale numeric, if provided, e.g. \code{mean(abs(x))}, 
#'                 it is used to scale up mu so that the scale is more friendly to the optimizer.
#'                 Default is 1. 
#'
#' @return numeric, linear working parameter array
#'
#' @keywords parameter
#'
#' @author Stephen H. Lihn
#'
#' @export ldhmm.n2w
#' 
#' @examples
#' param0 <- matrix(c(0.003, 0.02, 1, -0.006, 0.03, 1.3), 2, 3, byrow=TRUE)
#' gamma0 <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow=TRUE)
#' d <- ldhmm(m=2, param=param0, gamma=gamma0)
#' v <- ldhmm.n2w(d)
#'
### <======================================================================>
ldhmm.n2w <- function(object, mu.scale=1)
{
    m  <- object@m
    param.nbr <- object@param.nbr
    gamma <- object@gamma
    delta <- object@delta
    stationary <- object@stationary
    vparam <- as.vector(t(object@param))
    
    idx <- seq(0, m*param.nbr-1)
    
    # mu is unbounded, so no need to take log
    # sigma and lambda are bounded to positive numbers, so take log
    tparam <- ifelse(idx %% param.nbr == 0, vparam/mu.scale, log(abs(vparam)))
    if (m==1) return(tparam)
    
    # handle gamma, tgamma is length of m*(m-1), diag elemnts are removed
    lgamma <- log(gamma/diag(gamma))
    tgamma <- as.vector(lgamma[!diag(m)])
    
    # handle delta, tdetla is either NULL or length of m-1
    if (object@stationary) { tdelta  <- NULL }
    else { tdelta <- log(delta[-1]/delta[1]) }
    
    v <- c(tparam, tgamma, tdelta)
    return (v)
}
### <---------------------------------------------------------------------->
