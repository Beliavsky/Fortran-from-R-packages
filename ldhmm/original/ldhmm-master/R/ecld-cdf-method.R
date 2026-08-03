#' CDF and CCDF of ecld
#'
#' The analytic solutions for CDF and CCDF of ecld, if available.
#' \code{ecld.cdf_gamma} is a sub-module with the CDF expressed as
#' incomplete gamma function.
#' SGED is supported only in \code{ecld.cdf} and \code{ecld.ccdf}.
#'
#' @param object an object of ecld class
#' @param x a numeric vector of \code{x}
#'
#' @return The CDF or CCDF vector
#'
#' @keywords cdf
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.cdf
#' @export ecld.ccdf
#'
#' @importFrom gnorm pgnorm
#'
#' @examples
#' ld <- ecld(sigma=0.01)
#' x <- seq(-0.1, 0.1, by=0.01)
#' ecld.cdf(ld,x)
### <======================================================================>
"ecld.cdf" <- function(object, x)
{
    pgnorm(x, mu=object@mu, alpha=object@sigma, beta=2.0/object@lambda)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.cdf
"ecld.ccdf" <- function(object, x)
{
    1.0 - ecld.cdf(object, x)
}
