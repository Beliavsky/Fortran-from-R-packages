#' Computing the PDF per state given the observations
#' 
#' Computing the PDF per state given the observations. Only one of state or x 
#' can be a vector per call. 
#'
#' @param object an ldhmm object
#' @param state numeric, the states. 
#' @param x numeric, the observations.
#'
#' @return a vector or matrix of PDF. The dimension of matrix is state times x
#'
#' @keywords pdf
#'
#' @author Stephen H. Lihn
#'
#' @export 
#'
#' @importFrom parallel mclapply
#' 
### <======================================================================>
ldhmm.state_pdf <- function(object, state, x) {
    m <- object@m
    n <- length(x)
    ms <- length(state)
    
    # obtain all ecld objects at once
    lds <- ldhmm.state_ld(object, state)
    
    # j is the index of the state vector (not the state itself)
    get_pdf <- function(j, x) ecld.pdf(lds[[j]], x)
    
    # vector result
    if (n==1 || ms==1) return(sapply(1:ms, function(j) get_pdf(j,x)))
    
    # matrix result
    pdf <- matrix(NA, ms, n)
    ps <- parallel::mclapply(1:ms, function(j) get_pdf(j,x))
    for (j in 1:ms) {
        p <- ps[[j]]
        for (k in 1:n) pdf[j,k] <- p[k]
    }
    return(pdf)
}
### <---------------------------------------------------------------------->
