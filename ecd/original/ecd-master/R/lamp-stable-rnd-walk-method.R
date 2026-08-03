#' Calculate the stable random walk
#' 
#' Calculate the stable random walk. There are 4 types of random walk you can specify:
#' 11. Laplace(0,1). No skewess.
#' 1. Experimental Laplace random walk via Gauss-Laplace transmutation.
#' 22. Normal distribution N(0, sqrt(n))*epsilon. No skewess.
#' 2. Binomial random walk, b*epsilon. This can produce skewness.
#'
#' @param object an object of lamp class
#' @param n numeric, number of items in Levy sums
#' @param b numeric, cumulative sum of signs in Levy sums
#' 
#' @return numeric, the value of the random walk
#'
#' @keywords simulation
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#' 
### <======================================================================>
"lamp.stable_rnd_walk" <- function(object, n, b) {

    lambda <- object@lambda
    alpha <- object@alpha
    rnd.walk <- object@rnd.walk
    beta <- object@beta
    
    if (rnd.walk %in% c(1,2)) {
        # this only works for symmetric case
        if (!(beta==0 | beta==1)) stop("rnd.walk mode does not work with skewness scenario")
    }
    
    if (rnd.walk==1) return(rlaplace0(1,1))

    epsilon <- object@T.inf^(-1/lambda)
    if (rnd.walk==2) return(epsilon*rnorm(1,0,sqrt(n)))
    
    if (rnd.walk==11) {
        scale <- sqrt(rexp(1,1)*2) # it needs to generate sd=sqrt(2) and K=6
        return(b/sqrt(n)*scale)
    }
    
    if (rnd.walk==22) return(b*epsilon) # random walk of normal type
    
    stop(paste("ERROR: invalid rnd.walk scenario:", rnd.walk))
}
### <---------------------------------------------------------------------->
