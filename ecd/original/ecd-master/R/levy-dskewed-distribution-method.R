#' Skewed Levy distribution in Levy statistics
#' 
#' Skewed Levy distribution PDF. In our context, "skewed" means "completed asymmetric alpha-stable",
#' or called "one-sided alpha-stable". And we use lambda = 2/alpha.
#'
#' @param x numeric, complex, mpfr, mpfc
#' @param lambda numeric. Default is 4, the Levy distribution as is generally called.
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
#' y = levy.dskewed(x)
#'
### <======================================================================>
"levy.dskewed" <- function(x, lambda=4)
{
    if (lambda==4) {
        levy <- function(x) (2*sqrt(pi)*x^(3/2))^(-1)*exp(-1/4/x)
        return(levy(x))
    }
    stop("lambda is not supported")
    # TODO to implement either the series solution, or integral solution
}
### <---------------------------------------------------------------------->
