#' Constructing the ecld objects per state
#' 
#' This utility constructs the ecld objects per state and return them in a list of easy query. 
#'
#' @param object an ldhmm object
#' @param state numeric, the states. 
#'
#' @return a list of ecld objects
#'
#' @keywords pdf
#'
#' @author Stephen H. Lihn
#'
#' @export 
#' 
### <======================================================================>
ldhmm.state_ld <- function(object, state=NULL) {
    if (is.null(state)) state <- 1:object@m
    
    param.nbr <- object@param.nbr
    param <- object@param

    get_ld <- function(st) {
        mu <- param[st,1]
        sigma <- param[st,2]
        
        if (mu >= 1e5) mu <- 1e5
        if (mu <= -1e5) mu <- -1e5
        if (sigma >= 1e5) sigma <- 1e5
        if (sigma <= 1e-8) sigma <- 1e-8
        
        lambda <- if (param.nbr==3) param[st,3] else 1
        if (lambda <= 1e-2) lambda <- 1e-2
        if (lambda >= 10) lambda <- 10
        
        ecld(lambda=lambda, sigma=sigma, mu=mu)
    }
    sapply(state, get_ld)
}
### <---------------------------------------------------------------------->
