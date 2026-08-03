#' Analytic solutions on the statistics of quartic stable lambda distribution
#'
#' Several analytic solutions on the statistics of quartic stable lambda distribution (QSLD) are 
#' implemented. These functions provide precise validation on the distribution.
#'
#' @param x numeric, vector of responses.
#' @param nu numeric, vector of nu in the pdf integrand, starting from 0 (not nu0).
#' @param t numeric, the time parameter, where the variance is t, default is 1.
#' @param nu0 numeric, the location parameter, default is 0.
#' @param theta numeric, the scale parameter, default is 1.
#' @param convo numeric, the convolution number, default is 1.
#' @param beta.a numeric, the skewness parameter, default is 0. This number is annualized by sqrt(t).
#' @param mu numeric, the location parameter, default is 0.
#'
#' @return numeric
#'
#' @keywords stable-Lambda
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats rexp
#' @importFrom stats rnorm
#'
#' @export qsl_variance_analytic
#' @export qsl_kurtosis_analytic
#' @export qsl_skewness_analytic
#' @export qsl_std_pdf0_analytic
#' @export qsl_pdf_integrand_analytic
#' 
#' 
#' @references
#'   For more detail, see Appendix C of 
#'   Stephen Lihn (2017). 
#'   \emph{A Theory of Asset Return and Volatility 
#'   under Stable Law and Stable Lambda Distribution}.
#'   SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
#'   
#' @examples
#'   # obtain the variance for SPX 1-day distribution
#'   var <- qsl_variance_analytic(t=1/250, nu0=6.92/100, theta=1.17/100, convo=2, beta=-1.31)

### <======================================================================>
qsl_kurtosis_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    if (beta.a != 0) stop("qsl_kurtosis_analytic does not support non-zero beta")
    m <- convo
    beta <- beta.a*sqrt(t)
    nu1 <- nu0+6*theta
    A <- nu1^4 + (144+96*m)*theta^2*nu1^2 + 768*(1+m)*theta^3*nu1 + (4032+3456*m)*theta^4
    B <- (nu1^2+24*theta^2)^2
    return(3+(3/m)*A/B)
}
### <---------------------------------------------------------------------->
#' @rdname qsl_kurtosis_analytic
qsl_skewness_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    m <- convo
    beta <- beta.a*sqrt(t)
    nu1 <- nu0+6*theta

    P <- 144*beta^2 + 432 + 144*m*(2+beta^2)
    Q <- 192*beta^2*m^2 + 576*(2+beta^2)*m + 384*beta^2 + 1152

    A <- 2*(3+beta^2)*nu1^3 + P*theta^2*nu1 + Q*theta^3
    B <- (2+beta^2)*nu1^2 + 24*theta^2*(2+beta^2*(1+m))
    # print(paste(A,B^(3/2))) # debug this complex beast!
    return(beta/sqrt(m)*A/B^(3/2))
}
### <---------------------------------------------------------------------->
#' @rdname qsl_kurtosis_analytic
qsl_variance_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    m <- convo
    beta <- beta.a*sqrt(t)
    nu1 <- nu0+6*theta

    var_annl <- (nu1^2+24*theta^2) + m*beta^2/(2+beta^2)*(24*theta^2)
    return(var_annl*t)
}
### <---------------------------------------------------------------------->
#' @rdname qsl_kurtosis_analytic
qsl_std_pdf0_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    if (beta.a != 0) stop("qsl_std_pdf0_analytic does not support non-zero beta")

    if (convo > 100) convo <- ecd.mp1*convo # to ensure gamma function works

    m <- convo
    C <- gamma(m-1/2)*sqrt(2*m)/gamma(m)/8/pi
    n <- nu0/theta
    pdf0 <- 2*sqrt(pi)-pi*sqrt(n)*exp(n/4)*ecd.erfc(sqrt(n/4))
    var_annl <- (n+6)^2+24
    ecd.mp2f(C*pdf0*sqrt(var_annl)) # render numeric always
}
### <---------------------------------------------------------------------->
#' @rdname qsl_kurtosis_analytic
qsl_pdf_integrand_analytic <- function(x, nu, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    if (convo > 100) convo <- ecd.mp1*convo # to ensure gamma function works

    m <- convo
    beta <- beta.a*sqrt(t)
    B0 <- sqrt(1+beta^2/4)
    sm <- sqrt(t/m/(2+beta^2))
    w <- B0*(x-mu)/sm/(nu+nu0)

    C <- 4*pi*gamma(m)*sm*theta^(3/2)
    A <- sqrt(nu)/(nu+nu0)
    K <- abs(w/2/B0^2)^(m-1/2) * besselK(abs(w),m-1/2)
    B <- exp(w*beta/2/B0 - nu/4/theta)
    ecd.mp2f(1/C*A*K*B) # render numeric always
}
### <---------------------------------------------------------------------->
