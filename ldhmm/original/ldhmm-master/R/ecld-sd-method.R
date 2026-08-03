#' Compute statistics analytically for an ecld object
#' 
#' Compute statistics for mean, var, skewness, kurtosis.
#'
#' @param object an object of ecld class
#'
#' @return numeric
#'
#' @keywords statistics
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.sd
#' @export ecld.var
#' @export ecld.mean
#' @export ecld.skewness
#' @export ecld.kurt
#' @export ecld.kurtosis
#'
#' @examples
#' ld <- ecld(3)
#' ecld.sd(ld)
#' ecld.var(ld)
#' ecld.mean(ld)
#' ecld.skewness(ld)
#' ecld.kurt(ld)
#'
### <======================================================================>
"ecld.sd" <- function(object)
{
    sqrt(ecld.var(object))
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.var" <- function(object)
{
    s <- object@sigma
    l <- object@lambda
    
    # symmetric
    x <- gamma(l*3/2)
    y <- gamma(l/2)
    return(s^2*x/y)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.mean" <- function(object)
{
    return(object@mu)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.skewness" <- function(object)
{
    return(0)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.kurtosis" <- function(object)
{
    s <- object@sigma
    l <- object@lambda
    
    x <- gamma(l/2) * gamma(l*5/2)
    y <- gamma(l*3/2)^2
    return(x/y)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.kurt" <- function(object) ecld.kurtosis(object)
### <---------------------------------------------------------------------->



