#' Compute statistics analytically for an sld object
#' 
#' Compute statistics for mean, var, skewness, kurtosis for SLD. 
#' These functions are just wrappers on \code{\link{ksl}}.
#' If you need to calculate the statistics in quantity, you should use \code{\link{ksl}} or \code{\link{kqsl}} directly.
#'
#' @param object an object of sld class
#'
#' @return numeric
#'
#' @keywords statistics
#'
#' @author Stephen H-T. Lihn
#'
#' @export sld.sd
#' @export sld.var
#' @export sld.mean
#' @export sld.skewness
#' @export sld.kurt
#' @export sld.kurtosis
#'
#' @examples
#' d <- qsld(nu0=10.4, theta=1.6, convo=2)
#' sld.sd(d)
#' sld.var(d)
#' sld.mean(d)
#' sld.skewness(d)
#' sld.kurt(d)
#'
### <======================================================================>
"sld.sd" <- function(object)
{
    sqrt(sld.var(object))
}
### <---------------------------------------------------------------------->
#' @rdname sld.sd
"sld.var" <- function(object)
{
    .sld.k(object)[2]
}
### <---------------------------------------------------------------------->
#' @rdname sld.sd
"sld.mean" <- function(object)
{
    .sld.k(object)[1]
}
### <---------------------------------------------------------------------->
#' @rdname sld.sd
"sld.skewness" <- function(object)
{
    .sld.k(object)[5]
}
### <---------------------------------------------------------------------->
#' @rdname sld.sd
"sld.kurtosis" <- function(object)
{
    .sld.k(object)[6]
}
### <---------------------------------------------------------------------->
#' @rdname sld.sd
"sld.kurt" <- function(object) sld.kurtosis(object)
### <---------------------------------------------------------------------->
".sld.k" <- function(object)
{
    d <- object
    ksl(t=d@t, nu0=d@nu0, theta=d@theta, convo=d@convo, beta.a=d@beta.a, mu=d@mu, lambda=d@lambda)
}
### <---------------------------------------------------------------------->



