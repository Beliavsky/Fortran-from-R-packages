#' Constructor of lamp class
#' 
#' Construct an lamp class by providing the required parameters.
#' The default is the unit quartic lambda process.
#'
#' @param lambda numeric, the lambda parameter. Must be positive. Default is NaN.
#' @param alpha numeric, optional, if you don't like to use lambda. Default is NaN.
#'              Either lambda or alpha must be specified with a positive number.
#' @param beta  numeric, the skewness parameter. Default: 0.
#' @param rnd.walk numeric, random walk method, 1: Laplace, 2: Binomial/normal. Default is 1.
#' @param T.inf numeric, the infinite bound to cut off Levy sums. Default is 86400000.
#' @param rnd.n  numeric, the length of one rnd call. Default is 1000000.
#' @param sd  numeric, standard deviation adjustment. No adjustment if NaN. Default is NaN.
#' @param sd.method  numeric, methodology of sd adjustment. 0 means in scale parameter, 1 means in Levy sums. Default is 0.
#' @param N.lower  numeric, the lower bound of N to truncate the boundary effect. Default is 0.
#' @param N.upper  numeric, the upper bound of N to limit the outliers. Default is 1000.
#' @param file  character, file path to save the object and simulation result. Default is character(0).
#'
#' @return an object of lamp class
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export lamp
#'
#' @examples
#' lp <- lamp(4, T.inf=86400*1000000)
#'
### <======================================================================>
"lamp" <- function(lambda = NaN, T.inf = 86400*1000, rnd.n = 1000000,
                   alpha = NaN, beta = 0, rnd.walk = 1,
                   sd = NaN, sd.method = 0,
                   N.lower = 0, N.upper = 1000,
                   file = character(0))
{
    call <- match.call()

    if ((!is.na(alpha)) & is.na(lambda)) lambda <- 2/alpha
    if(!is.numeric(lambda)) stop("Parameter 'lambda' must be a numberic!\n")
    if(!is.numeric(beta)) stop("Parameter 'beta' must be a numberic!\n")

    # -------------
    if(!(lambda > 0)) stop("Parameter 'lambda' must be positive!\n")
    if(!(beta >= -1 & beta <= 1)) stop("Parameter 'beta' must be between -1 and 1!\n")
    if(!(rnd.walk %in% c(1,11,2,22))) stop("Parameter 'rnd.walk' is invalid")
    if(!(N.lower >= 0)) stop("Parameter 'N.lower' must be 0 or positive")
    if(!(N.upper >= 0)) stop("Parameter 'N.upper' must be positive")
    if(!(N.upper >= N.lower)) stop("Parameter 'N.upper' must be larger than 'N.lower'")
    
    # MPFR is specified indirectly through T.inf
    use.mpfr <- ifelse(class(T.inf)=="mpfr", TRUE, FALSE)
    if (use.mpfr) T.inf <- ecd.mp1 * T.inf

    # -------------
    lp <- new("lamp", call = call,
              lambda = unname(lambda),
              alpha = unname(2/lambda),
              beta  = unname(beta),
              pm = 1, # no other pm allowed
              rnd.walk = unname(rnd.walk),
              sd = unname(sd),
              sd.method = unname(sd.method),
              T.inf = unname(T.inf),
              rnd.n = unname(rnd.n),
              N.lower = unname(N.lower),
              N.upper = unname(N.upper),
              use.mpfr = unname(use.mpfr),
              file = unname(file)
          )
    
    # -------------
    invisible(lp)
}
### <---------------------------------------------------------------------->
