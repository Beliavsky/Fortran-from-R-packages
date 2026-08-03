#' Calculate sd adjustment factor
#' 
#' Calculate sd adjustment factor. 
#' For L2 random walk, it is the power of 1/(1+alpha/2).
#' For L1 random walk, it is the power of 1.
#' This factor can be used to adjust either the scale parameter of the stable distribution
#' or T.inf that cuts off the Levy sums.
#'
#' @param object an object of lamp class
#'
#' @return numeric, the sd factor
#'
#' @keywords simulation
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#' 
### <======================================================================>
"lamp.sd_factor" <- function(object) {
    
    if (is.na(object@sd)) return(1) # do nothing
    
    lambda <- object@lambda
    alpha <- object@alpha
    rnd.walk <- object@rnd.walk
    sd <- object@sd
    if (sd < 0) stop("sd must be positive in lamp")
    
    sd.factor <- 1 # sd adjustment factor
    sd0 <- ecld.sd(ecld(lambda))
    if (rnd.walk %in% c(1,11))  sd.factor <- (object@sd/sd0)
    else if (rnd.walk %in% c(2,22)) sd.factor <- (object@sd/sd0)^(1/(1+alpha/2))
    else stop("sd.factor: unsupported rnd.walk")

    return(sd.factor)
}
### <---------------------------------------------------------------------->
