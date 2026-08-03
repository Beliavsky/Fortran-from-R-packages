#' Constructor of ldhmm class
#' 
#' Construct an ldhmm class by providing the required parameters.
#'
#' @param m numeric, number of states
#' @param param matrix, the ecld parameters of states.
#' @param gamma numeric or matrix, the transition probability matrix, must be conformed to m by m.
#'              if provided as vector, it will be converted to a matrix with \code{byrow=TRUE}.
#' @param delta numeric, the initial distribution for each state, default is \code{NULL}.
#' @param stationary logical, specify whether the initial distribution is stationary or not,
#'                   default is \code{TRUE}.
#' @param mle.optimizer character, specify alternative optimizer, default is \code{nlm}.
#'
#' @return An object of ldhmm class
#'
#' @keywords constructor
#'
#' @author Stephen H. Lihn
#'
#' @export ldhmm
#' 
#' @examples
#' param0 <- matrix(c(0.003, 0.02, 1, -0.006, 0.03, 1.3), 2, 3, byrow=TRUE)
#' gamma0 <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow=TRUE)
#' d <- ldhmm(m=2, param=param0, gamma=gamma0)
#'
### <======================================================================>
"ldhmm" <- function(m, param, gamma, delta=NULL, stationary=TRUE, mle.optimizer="nlm")
{
    call <- match.call()
    param.nbr <- NCOL(param)
    
    if (length(m) != 1 | m %% 1 != 0) stop("m must be an integer")
    if (param.nbr %% 1 != 0) stop("param.nbr must be integer")
    if (param.nbr != 2 & param.nbr != 3) stop("param.nbr must be 2 or 3")
    if (NROW(param) != m) stop("param must have m rows")
    if (stationary == FALSE) {
        if (length(delta) != m) stop("delta must be length-m vector when stationary is not true")
    }

    if (is(gamma, "numeric")) {
        if(length(gamma) != m*m) stop(paste("gamma is a vector, but length is not m by m:", length(gamma)))
        gamma <- ldhmm.gamma_init(m, prob=gamma)
    }
    if (NROW(gamma) != m | NCOL(gamma) != m) stop("gamma must be mxm matrix")
    
    # give param column names
    if (param.nbr==2) colnames(param) <- c("mu", "sigma")
    if (param.nbr==3) colnames(param) <- c("mu", "sigma", "lambda")
    
    d <- new("ldhmm",
             call = call,
             m = unname(m),
             param.nbr = param.nbr,
             param = param,
             gamma = gamma,
             delta = delta,
             stationary = stationary,
             mle.optimizer = mle.optimizer
            )
       
    invisible(d)
}
### <---------------------------------------------------------------------->
