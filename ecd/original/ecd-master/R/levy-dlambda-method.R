#' Standard Lambda distribution
#' 
#' Standard Lambda distribution PDF that can take complex argument.
#'
#' @param x numeric, complex, mpfr, mpfc
#' @param lambda numeric. Default is 4, the quartic distribution.
#' 
#' @return PDF in the same type as x
#'
#' @keywords PDF
#'
#' @author Stephen H. Lihn
#'
#' @export 
#' 
#' @examples
#' x = seq(1,10)
#' y = levy.dlambda(x)
#'
### <======================================================================>
"levy.dlambda" <- function(x, lambda=4)
{
    C <- lambda*gamma(lambda/2)
    exp(-(x^2)^(1/lambda))/C
}
### <---------------------------------------------------------------------->
