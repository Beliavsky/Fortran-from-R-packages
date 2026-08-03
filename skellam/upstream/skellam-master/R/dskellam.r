#' The Skellam Distribution
#'
#' Density, distribution function, quantile function, and random generation for the Skellam distribution.
#'
#' The Skellam distribution describes the difference between two independent Poisson random variables.
#' This documentation covers:
#'
#' \strong{Density:} \cr
#' \preformatted{
#' dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)
#' }
#'
#' \strong{Distribution Function:} \cr
#' \preformatted{
#' pskellam(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
#' }
#'
#' \strong{Quantile Function:} \cr
#' \preformatted{
#' qskellam(p, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
#' }
#'
#' \strong{Random Generation:} \cr
#' \preformatted{
#' rskellam(n, lambda1, lambda2 = lambda1)
#' }
#'
#' \strong{Saddlepoint Approximations:} \cr
#' \preformatted{
#' dskellam.sp(x, lambda1, lambda2 = lambda1, log = FALSE)
#' pskellam.sp(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
#' }
#'
#'
#' @details
#' If \eqn{Y_1} and \eqn{Y_2} are Poisson variables with means \eqn{\mu_1} and \eqn{\mu_2} and correlation \eqn{\rho},
#' then \eqn{X = Y_1 - Y_2} is Skellam with parameters:
#'
#' \deqn{\lambda_1 = \mu_1 - \rho \sqrt{\mu_1 \mu_2}}
#'
#' \deqn{\lambda_2 = \mu_2 - \rho \sqrt{\mu_1 \mu_2}}
#'
#' The density is given by:
#'
#' \deqn{I(2 \sqrt{\lambda_1 \lambda_2}, |x|) (\lambda_1/\lambda_2)^{x/2} \exp(-\lambda_1-\lambda_2)}
#'
#' where \eqn{I(y,\nu)} is the modified Bessel function of the first kind.
#'
#' @param x,q For functions \code{dskellam}, \code{dskellam.sp}, and \code{pskellam.sp}: a numeric vector of quantiles.
#' @param p For \code{qskellam}: a numeric vector of probabilities.
#' @param n For \code{rskellam}: a non-negative integer specifying the number of observations.
#' @param lambda1,lambda2 Numeric vectors of (non-negative) means; \code{lambda2} defaults to \code{lambda1} if not provided.
#' @param log,log.p Logical; if TRUE, returns the logarithm of the computed value.
#' @param lower.tail Logical; if TRUE (default), returns \eqn{P(X \le x)}; otherwise, returns \eqn{P(X > x)}.
#'
#' @return
#' \itemize{
#'   \item \code{dskellam} returns the (log) density.
#'   \item \code{pskellam} returns the (log) cumulative distribution function.
#'   \item \code{qskellam} returns the quantile function.
#'   \item \code{rskellam} generates random deviates.
#' }
#'
#' Invalid \code{lambda} values will return \code{NaN} with a warning.
#'
#' @note
#' The \pkg{VGAM} package also provides Skellam functions. This implementation offers a broader working range,
#' correct handling when one rate parameter is zero, enhanced argument checking, and improved accuracy for \eqn{x < 0}
#' (in R versions prior to 2.9). Use \code{skellam::dskellam} or \code{VGAM::dskellam} to specify which implementation to use.
#'
#' @references
#' \itemize{
#'   \item Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
#'         Cambridge University Press.
#'   \item Johnson, N. L. (1959) On an extension of the connection between Poisson and \eqn{\chi^2} distributions.
#'         \emph{Biometrika} 46, 352-362.
#'   \item Johnson, N. L., Kotz, S., & Kemp, A. W. (1993) \emph{Univariate Discrete Distributions},
#'         2nd ed., John Wiley and Sons.
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates.
#'         \emph{Journal of the Royal Statistical Society, Series A} 109(3), 296.
#'   \item Strackee, J. & van der Gon, J. J. D. (1962) The frequency distribution of the difference
#'         between two Poisson variates. \emph{Statistica Neerlandica} 16(1), 17-23.
#'   \item Wikipedia: \url{https://en.wikipedia.org/wiki/Skellam_distribution}
#' }
#'
#' @examples
#' # Compare with Poisson when one lambda = 0
#' dskellam(0:10, 5, 0)
#' dpois(0:10, 5)
#'
#' # Both lambdas non-zero
#' dskellam(c(-1,1), c(12,10), c(10,12))
#' pskellam(c(-1,0), c(12,10), c(10,12))
#'
#' # Quantile function
#' qskellam(c(0.05, 0.95), 3, 4)
#'
#' # Random generation
#' rskellam(10, 8.5, 10.25)
#'
#' @name skellam
#' @aliases dskellam pskellam qskellam rskellam dskellam.sp pskellam.sp
NULL

NULL


#' @export
#' @rdname skellam
dskellam <- function(x, lambda1, lambda2=lambda1, log=FALSE){
 # density (PMF) of Skellam distriubition (difference of Poissons)
    if (missing(x)|missing(lambda1)) stop("first 2 arguments are required")
    lambdas <- c(lambda1,lambda2)
    oops <- !(is.finite(lambdas)&(lambdas>=0))
    if(any(oops)) {
        warning("NaNs produced")
        lambdas[oops] <- NaN
        lambda1 <- lambdas[1:length(lambda1)]
        lambda2 <- lambdas[(length(lambda1)+1):length(lambdas)]
    }
   # make all args the same length (for subsetting)
    lens <- c(length(x),length(lambda1),length(lambda2))
    len <- max(lens)
    if(len>min(lens)) {
        if (all(len/lens!=len%/%lens)) warning("longer object length is not a multiple of shorter object length", domain=NA)
        x <- rep(x,length.out=len)
        lambda1 <- rep(lambda1,length.out=len)
        lambda2 <- rep(lambda2,length.out=len)
    }
   # warn of non-integer x values (since support of PMF is integers)
    nonint <- x!=trunc(x)
    if (any(nonint)) {
        xreal <- x[nonint]
        for (i in 1:length(xreal)) warning(paste("non-integer x =",xreal[i]))
    }
    x <- trunc(x)
    ret <- rep(NaN,length.out=len)
   # handle a zero lambda separately (avoids division by zero & accuracy issues for large values of lambda or x)
    ret[lambda1==0] <- stats::dpois(-x[lambda1==0],lambda2[lambda1==0],log=log)
    ret[lambda2==0] <- stats::dpois( x[lambda2==0],lambda1[lambda2==0],log=log)    # corrects bug in VGAM 0.7-9
   # non-zero lambdas
    nz <- is.finite(lambda1)&is.finite(lambda2)&(lambda1>0)&(lambda2>0)
    L1 <- lambda1[nz];
    L2 <- lambda2[nz];
    sqL12 <- sqrt(L1*L2)
    xx <- x[nz]
    if (log[1]) {       # use abs(x) in besselI for improved accuracy (prior to 2.9) and S-PLUS compatibility
       # log(besselI(y,nu)) == y+log(besselI(y,nu,TRUE))
        ret[nz] <- log(besselI(2*sqL12, abs(xx),TRUE))+2*sqL12-L1-L2+xx/2*log(L1/L2)
    } else {
       # besselI(y,nu); exp(y)*besselI(y,nu,TRUE)
#       ret[nz] <- besselI(2*sqL12, abs(xx),TRUE)*(L1/L2)^(xx/2)*exp(2*sqL12-L1-L2)
        ret[nz] <- besselI(2*sqL12, abs(xx),TRUE)*exp(2*sqL12-L1-L2+xx/2*log(L1/L2))
    }
    chk <- nz & (!is.finite(ret) | (!log)&(ret<1e-308)) # use saddlepoint approximation to detect possible over/underflow
    if (length(chk[chk])>0) {
        L1 <- lambda1[chk];
        L2 <- lambda2[chk];
            sqL12 <- sqrt(L1*L2)
        xx <- x[chk]
        s <- log(0.5*(xx+sqrt(xx^2+4*L1*L2))/L1)# the saddlepoint
        K <- L1*(exp(s)-1)+L2*(exp(-s)-1)       # CGF(s)
        K2 <- L1*exp(s)+L2*exp(-s)              # CGF''(s)
        spd <- exp(K-xx*s)/sqrt(2*pi*K2)         # saddlepoint density, was x in place of xx
        usp <- (spd>1e-308)&is.finite(spd)      # don't trust the existing result
        if (length(usp[usp])>0) {   # add another term to the saddlepoint approximation
            su <- s[usp]
            K2u <- K2[usp]
            c <- (1-((L1[usp]*exp(su)-L2[usp]*exp(-su))/K2u)^2*5/3)/K2u*0.125+1
            ret[chk][usp] <- exp(K[usp]-x[usp]*su)/sqrt(2*pi*K2u)*(1+c)*0.5
        }
    }
    ret
}
