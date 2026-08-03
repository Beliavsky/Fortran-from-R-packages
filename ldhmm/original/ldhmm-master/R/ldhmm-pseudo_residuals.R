#' Computing pseudo-residuals
#' 
#' This utility computes  pseudo-residuals. (Zucchini, 6.2)
#'
#' @param object an ldhmm object
#' @param x numeric, the observations.
#' @param xc.length a positive integer specifying the length of \code{xc} 
#'                  when calculating conditional probabilities, default is 1000.
#'
#' @return a vector of normal quantiles
#'
#' @keywords residuals
#'
#' @author Stephen H. Lihn
#'
#' @importFrom stats qnorm
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#'   sr <- ldhmm.pseudo_residuals(object, x)
#'   hist(sr)
#'   acf(sr)
#'   qqnorm(sr, cex=0.5)
#'   L <- seq(-3,3,length.out=100)
#'   lines(L,L,col="red",lwd=2, lty=2)
#' }
### <======================================================================>
ldhmm.pseudo_residuals <- function(object, x, xc.length=1000)
{
    n <- length(x)
    xc <- seq(min(x)*1.25, max(x)*1.25, length.out=xc.length)
    dxc <- diff(xc)[1]
    cdists <- ldhmm.conditional_prob(object, x, xc) 
    cumdists <- rbind(rep(0, n), apply(cdists, 2, cumsum)) * dxc
    ulo <- uhi <- rep(NA,n)
    for (i in 1:n) {
        j <- max(which(xc <= x[i]))
        ulo[i]  <- cumdists[j,i]
        uhi[i]  <- cumdists[j+1,i]
    }
    umi <- 0.5*(ulo+uhi)
    sr <- stats::qnorm(umi)
    return(sr)
}
### <---------------------------------------------------------------------->
