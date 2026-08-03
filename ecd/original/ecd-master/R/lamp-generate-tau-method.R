#' Generate tau from stable distribution
#' 
#' Generate tau, a random sequence representing the stable random walk process.
#'
#' @param object an object of lamp class
#'
#' @return an object of lamp class with tau populated, tau_i is set to 1.
#'
#' @keywords simulation
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#' lp <- lamp(4, rnd.n=10)
#' lp1 <- lamp.generate_tau(lp)
#' lp1@tau
#' 
### <======================================================================>
"lamp.generate_tau" <- function(object) {
    
    lambda <- object@lambda
    alpha <- object@alpha
    n <- object@rnd.n
    
    if (object@beta==1 & lambda<=2) stop("lambda=2 beta=1 is not supported")

    sd.factor <- if (object@sd.method==0) lamp.sd_factor(object) else 1
    
    g <- if (alpha < 1) cos(pi*alpha/2)^(1/alpha) else 1
    if (alpha==1) g <- 1/8 # TODO guestimate, why does it work?
    
    if (object@beta==1) {
        rs <- rstable(n, alpha, 1, gamma=g/sd.factor, pm=object@pm)
        bi <- sign(runif(n)-0.5)
        object@tau <- rs*bi # embed the sign into one-sided Levy
        object@tm <- Sys.time()
        return(object)
    }
    
    # S(alpha, beta; pm=1)    
    object@tau <- rstable(n, alpha, object@beta, gamma=g/sd.factor, pm=object@pm) 
    object@tau_i <- 1
    object@tm <- Sys.time()
    return(object)
}
### <---------------------------------------------------------------------->
