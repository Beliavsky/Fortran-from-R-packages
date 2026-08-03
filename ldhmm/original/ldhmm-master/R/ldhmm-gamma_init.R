#' Initializing tansition probability paramter
#' 
#' This utility has multiple purposes. It can generate a simple transition probability matrix, 
#' using p1 and p2, if prob is left as NULL. The generated gamma is raw and not normalized.
#' If prob is provided as a vector, the utility converts it into a matrix as gamma.
#' Furthermore, if prob is provided as a vector or matrix, the utility 
#' applies \code{min.gamma}, and normalize the sum of t.p.m. rows to 1.
#' This is mainly an internal function used by MLE, not be concerned by external users.
#'
#' @param m numeric, number of states
#' @param p1 numeric, the first-neighbor transition probability, default is 0.04.
#' @param p2 numeric, the second-neighbor transition probability, default is 0.01.
#' @param prob numeric or matrix, a fully specified transition probability by user, default is \code{NULL}.
#'             If this is specified, p1, p2 would be ignored.
#' @param min.gamma numeric, a minimum transition probability added to gamma to avoid singularity, default is 0.
#'                  This is only used when \code{prob} is not \code{NULL}. 
#'
#' @return a matrix as gamma
#'
#' @keywords constructor
#'
#' @author Stephen H. Lihn
#'
#' @export 
#' 
#' @examples
#'   gamma0 <- ldhmm.gamma_init(m=3)
#'   prob=c(0.9, 0.1, 0.1, 
#'          0.1, 0.9, 0.0,
#'          0.1, 0.1, 0.8)
#'   gamma1 <- ldhmm.gamma_init(m=3, prob=prob)
#'   gamma2 <- ldhmm.gamma_init(m=2, prob=gamma1, min.gamma=1e-6)
#'
### <======================================================================>
ldhmm.gamma_init <- function(m, p1=0.04, p2=0.01, prob=NULL, min.gamma=0) {
    gamma <- matrix(NA, nrow=m, ncol=m, byrow=TRUE)
    if (! is.null(prob)) {
        if (is(prob, "numeric")) {
            gamma <- matrix(prob, nrow=m, ncol=m, byrow=TRUE)
        } 
        else if (inherits(prob, "matrix")) gamma <- prob
        else stop("prob must be a numeric or matrix")
        
        if (! inherits(gamma, "matrix")) stop("gamma must be a matrix")
        
        # normalize gamma so that transition probability sums up to one for each state
        # replace all zeros with min.gamma so that there is no singularity
        gamma <- ifelse(abs(gamma)==0, gamma + min.gamma, gamma) 
        one <- sapply(1:m, function(i) sum(gamma[i,]))
        for (i in 1:m) gamma[i,] <- gamma[i,]/one[i]
        return(gamma)
    }
    
    # prob = NULL, generate a simple t.p.m based on p1 and p2
    if (m==2) return(matrix(c(1-p1, p1, p1, 1-p1), m, m, byrow=TRUE))
    
    # sophisticated case
    t <- c(p2, p1, 0, p1, p2)
    if (m > 3) t <- c(0*(1:(m-3)),t)
    t <- c(t, 0*(1:m))

    for (i in 1:m) {
        j <- m-i+1
        t1 <- t[j:(j+m-1)]
        gamma[i,] <- t1
        gamma[i,i] <- 1-sum(t1)
    }
    return (gamma)
}
### <---------------------------------------------------------------------->

