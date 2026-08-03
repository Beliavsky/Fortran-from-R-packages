#' Constructor of ecld class
#' 
#' Construct an \code{\link{ecld-class}} by providing the required parameters.
#' The default is the standard symmetric cusp distribution (lambda=3).
#'
#' @param lambda numeric, the lambda parameter. Must be positive. Default: 3.
#' @param sigma numeric, the scale parameter. Must be positive. Default: 1.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param verbose logical, display timing information, for debugging purpose, default is \code{FALSE}.
#'
#' @return an object of ecld class
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld
#'
#' @examples
#' ld <- ecld()
#' ld <- ecld(2, 0.01)

### <======================================================================>
"ecld" <- function(lambda = 3, sigma = 1, mu = 0, verbose=FALSE)
{
    call <- match.call()
    # -------------
    if(sigma <= 0){
        stop("Parameter 'sigma' must be positive!\n")
    }
    if(lambda <= 0){
        stop("Parameter 'lambda' must be positive!\n")
    }

    ld <- new("ecld", call = call,
               lambda = unname(lambda),
               sigma = unname(sigma),
               mu    = unname(mu)
          )
    # -------------
    if (verbose) print(paste(Sys.time(), "ecld constructor: done"))
   
    invisible(ld)
}

