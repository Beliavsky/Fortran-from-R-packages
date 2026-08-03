#' Computing the conditional probabilities
#' 
#' This utility computes the conditional probabilities that observation at time t 
#' equals xc, given all observations other than that at time t being the same.
#'
#' @param object an ldhmm object
#' @param x numeric, the observations.
#' @param xc numeric, the conditional observations.
#'
#' @return matrix of probabilities, size of xc times size of x.
#'
#' @keywords pdf
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
### <======================================================================>
ldhmm.conditional_prob <- function(object, x, xc)
{
    n         <- length(x)
    m         <- object@m
    nxc       <- length(xc)
    
    pdf <- ldhmm.state_pdf(object, 1:m, xc)
    if (is(pdf, "numeric")) pdf <- as.matrix(pdf, nrow=m, ncol=nxc)

    la        <- ldhmm.log_forward(object, x)
    lb        <- ldhmm.log_backward(object, x)
    la        <- cbind(log(object@delta), la)
    lafact    <- apply(la, 2, max)
    lbfact    <- apply(lb, 2, max)
    
    dxc <- matrix(NA, nrow=nxc, ncol=n)
    for (i in 1:n) {
        a <- exp(la[,i]-lafact[i])
        b <- exp(lb[,i]-lbfact[i])
        foo <- a %*% object@gamma * b
        foo <- foo/sum(foo)
        dxc[,i] <- foo %*% pdf
    }
    return(dxc)
}
### <---------------------------------------------------------------------->

