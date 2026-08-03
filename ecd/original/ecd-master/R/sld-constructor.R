#' Constructor of sld class
#' 
#' Construct an \code{\link{sld-class}} by providing the required parameters.
#' The \code{qsld} constructor is also provided by hiding and defaulting lambda to 4.
#'
#' @param t numeric, the time parameter. Must be positive. Default: 1.
#' @param nu0 numeric, the floor volatility parameter. Must be positive. Default: 0.
#' @param theta numeric, the volatility scale parameter. Must be positive. Default: 1.
#' @param convo numeric, the convolution parameter. Must be positive. Default: 1.
#' @param beta.a  numeric, the skewness parameter. Default: 0.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param lambda numeric, the lambda parameter. Must be positive. Default: 4.
#'
#' @return an object of sld class
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export sld
#' @export qsld
#'
#' @examples
#' d <- sld()
#' d <- qsld()
#'
### <======================================================================>
"sld" <- function(t = 1, nu0 = 0, theta = 1, convo = 1, beta.a = 0, mu = 0,
                  lambda = 4)
{
    call <- match.call()

    # -------------
    if(t <= 0) stop("Parameter 't' must be positive!\n")
    if(convo <= 0) stop("Parameter 'convo' must be positive!\n")
    if(lambda <= 0) stop("Parameter 'lambda' must be positive!\n")

    if(nu0 < 0) stop("Parameter 'nu0' must be zero or positive!\n")
    if(theta < 0) stop("Parameter 'theta' must be zero or positive!\n")
    if(nu0+theta <= 0) stop("Parameters nu0 plus theta must be positive!\n")
    
    sum <- t+nu0+theta+convo+beta.a+mu
    if(!(abs(sum)>=0 & abs(sum) != Inf)) {
        stop("Parameters must be finite, known real numbers!")
    }

    sld <- new("sld", call = call,
               t = unname(t),
               nu0 = unname(nu0),
               theta = unname(theta),
               convo = unname(convo),
               beta.a  = unname(beta.a),
               mu    = unname(mu),
               lambda = unname(lambda)
          )
    
    invisible(sld)
}
### <---------------------------------------------------------------------->
#' @rdname sld
"qsld" <- function(t = 1, nu0 = 0, theta = 1, convo = 1, beta.a = 0, mu = 0)
{
    sld(t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4) 
}
### <---------------------------------------------------------------------->
