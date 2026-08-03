#' Computing the global decoding by the Viterbi algorithm
#' 
#' This utility computes the global decoding by the Viterbi algorithm.
#'
#' @param object an ldhmm object
#' @param x numeric, the observations.
#'
#' @return a vector of states
#'
#' @keywords viterbi
#'
#' @author Stephen H. Lihn
#'
#' @export ldhmm.viterbi
#' 
### <======================================================================>
ldhmm.viterbi <- function(object, x)
{
    m <- object@m
    n <- length(x)
    xi <- matrix(0, n, m)
    pdf <- ldhmm.state_pdf(object, 1:m, x)
    if (is(pdf, "numeric")) pdf <- as.matrix(pdf, nrow=m, ncol=n)

    phi <- object@delta * pdf[,1]
    xi[1,] <- phi/sum(phi)
    for (i in 2:n) {
        phi <- apply(xi[i-1,] * object@gamma, 2, max) * pdf[,i]
        xi[i,] <- phi/sum(phi)
    }
    iv <- numeric(n)
    iv_n <- which.max(xi[n,]) # this could return nothing in some cases
    if (length(iv_n) > 0) iv[n] <- iv_n
    for (i in (n-1):1) {
        if (!is.na(iv[i+1])) {
            iv_i <- which.max(object@gamma[,iv[i+1]] * xi[i,])
            if (length(iv_i) > 0) iv[i] <- iv_i # this could return nothing in some cases
        }
    }
    return(iv)
}
### <---------------------------------------------------------------------->

