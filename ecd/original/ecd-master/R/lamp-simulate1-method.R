#' Simulate one sequence of lambda process from stable distribution
#' 
#' Simulate lambda process from one random sequence representing the stable random walk process.
#'
#' @param object an object of lamp class
#' @param drop  numeric, number of tau to discard at the end. Default is 10.
#' @param keep.tau  numeric, 0 to clean up, 1 to return unused tau, 2 to return all tau. Default is 1.
#'
#' @return an object of lamp class with Z, B, N populated
#'
#' @keywords simulation
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#' 
#' @examples
#' lp <- lamp(4, T.inf=8640, rnd.n=100000)
#' lp1 <- lamp.simulate1(lp)
#' 
### <======================================================================>
"lamp.simulate1" <- function(object, drop=10, keep.tau=1) {

    lambda <- object@lambda
    alpha <- object@alpha
    beta <- object@beta
    T.inf <- object@T.inf 
    
    sd.factor <- if (object@sd.method==1) lamp.sd_factor(object) else 1
    T.inf2 <- T.inf * sd.factor # this determines the cut
    
    # -------------------------------------------
    p <- lamp.generate_tau(object)

    Z <- B <- N <- numeric(0)
    Z_i <- tau_i <- 1
    Tn <- n <- b <- 0

    one <- if (object@use.mpfr) ecd.mp1 else 1
    len <- object@rnd.n
    while (tau_i < len - drop) {
        t <- p@tau[tau_i] *one # we only want to index it once
        # skip large tau
        if (abs(t) > T.inf2) {
            # print(paste("extreme value encountered:", tau_i, floor(tau[tau_i]/T.inf)))
            tau_i <- tau_i+1
            next
        }
        Tn <- Tn + abs(t)
        n <- n + 1
        b <- b + sign(t)
        tau_i <- tau_i+1
        while (Tn >= T.inf2) {
            n <- n-1 # unwind
            b <- b - sign(t) # unwind
            cnt <- n^(lambda/2)/T.inf # TODO why is this not T.inf2?
            if (n > 0 & cnt >= object@N.lower & cnt <= object@N.upper) {
                N[Z_i] <- cnt
                B[Z_i] <- lamp.stable_rnd_walk(object, n, b)
                Z[Z_i] <- B[Z_i]*cnt 
                Z_i <- Z_i + 1
            }
            Tn <- Tn - T.inf2
            b <- sign(t)
            n <- if (Tn > 0) 1 else 0
            p@tau_i <- tau_i
        }
    }
    
    p@Z_i <- length(Z)
    p@Z <- Z
    p@B <- B
    p@N <- N
    p@tm <- Sys.time()

    p@tau_i <- tau_i
    if (keep.tau==0) p@tau <- numeric(0) # clean up
    if (keep.tau==1) p@tau <- p@tau[tau_i:len] # return unused tau
    
    return(p)
}
### <---------------------------------------------------------------------->
