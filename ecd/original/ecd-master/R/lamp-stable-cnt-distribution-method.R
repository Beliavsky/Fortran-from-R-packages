#' Stable Count distribution
#'
#' Implements the stable count distribution
#' (based on stabledist package) for stable random walk simulation.
#' Quartic stable distribution is implemented through gamma distribution.
#'
#' @param n numeric, number of observations.
#' @param x numeric, vector of responses.
#' @param q numeric, vector of quantiles.
#' @param s numeric, vector of responses for characteristic function.
#' @param alpha numeric, the shape parameter, default is NULL. User must provide
#'              either alpha or lambda.
#' @param nu0 numeric, the location parameter, default is 0.
#' @param theta numeric, the scale parameter, default is 1.
#' @param lambda numeric, alternative shape parameter, default is NULL.
#'
#' @return numeric, standard convention is followed:
#'         d* returns the density,
#'         p* returns the distribution function,
#'         q* returns the quantile function, and
#'         r* generates random deviates.
#'         The following are our extensions:
#'         k* returns the first 4 cumulants, skewness, and kurtosis,
#'         cf* returns the characteristic function.
#'
#' @keywords Stable
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stabledist dstable
#' @importFrom stabledist pstable
#' @importFrom stats dgamma
#' @importFrom stats pgamma
#' @importFrom stats rgamma
#' @importFrom stats qgamma
#'
#' @export dstablecnt
#' @export pstablecnt
#' @export qstablecnt
#' @export rstablecnt
#' @export cfstablecnt
#' @export kstablecnt
#'
#' @section Details:
#'   The stable count distribution is the conjugate prior of the stable distribution.
#'   The density function is defined as \deqn{
#'     \mathit{N}_{\alpha}\left(\nu;\nu_{0},\theta\right) 
#'     \equiv\frac{\alpha}{\mathit{\Gamma}\left(\frac{1}{\alpha}\right)}\,
#'     \frac{1}{\nu-\nu_{0}}\,L_{\alpha}\left(\frac{\theta}{\nu-\nu_{0}}\right),
#'     \:\mathrm{where}\,\nu>\nu_{0}.
#'   }{
#'     N_\alpha(\nu; \nu_0, \theta) = \alpha/\Gamma(1/\alpha) * 
#'     1/(\nu-\nu_0) * L_\alpha(1/(\nu-\nu_0))
#'   }
#'   where \eqn{\nu>\nu_0}. \eqn{\alpha} is the stability index, 
#'   \eqn{\nu_0} is the location parameter, and \eqn{\theta} is the scale parameter.
#'   \cr
#'   At \eqn{\alpha=0.5} aka \eqn{\lambda=4}, it is called "quartic stable count distribution",
#'   which is a gamma distribution with shape of 3/2. It has the closed form of \deqn{
#'     \mathit{N}_{\frac{1}{2}}\left(\nu;\nu_{0},\theta\right)
#'     \equiv\frac{1}{4\sqrt{\pi}\,\theta^{3/2}}
#'     \left(\nu-\nu_{0}\right)^{\frac{1}{2}}
#'     e^{-\frac{\nu-\nu_{0}}{4\theta}}
#'   }{
#'     N_\alpha(\nu; \nu_0, \theta) = 1/(4 sqrt(\pi) \theta^1.5)
#'     (\nu-\nu_0)^0.5 exp(-(\nu-\nu_0)/(4\theta))
#'   }
#' 
#' @references
#'   For more detail, see Section 2.4 and Section 3.3 of 
#'   Stephen Lihn (2017). 
#'   \emph{A Theory of Asset Return and Volatility 
#'   under Stable Law and Stable Lambda Distribution}.
#'   SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
#'   This distribution is also documented formally 
#'   in Wikipedia: \url{https://en.wikipedia.org/wiki/Stable_count_distribution}.
#'   
#' @examples
#'   # generate the pdf of the VIX distribution
#'   x <- c(0, 100, by=0.1)
#'   pdf <- dstablecnt(x, nu0=10.4, theta=1.6, lambda=4)
#'   
### <======================================================================>
dstablecnt <- function(x, alpha=NULL, nu0=0, theta=1, lambda=NULL) {
    if (is.null(alpha) & !is.null(lambda)) alpha <- 2/lambda
    .stablecnt.check(alpha, nu0, theta)
    if (alpha==1/2) {
        x0 <- ifelse(x >= nu0, x-nu0, NaN)
        return(stats::dgamma(x0, shape=3/2, scale=4*theta))
    }

    x0 <- ifelse(x > nu0, x-nu0, NaN)
    c <- alpha / gamma(1/alpha)
    g <- cos(pi*alpha/2)^(1/alpha)
    c/x0*stabledist::dstable(theta/x0, alpha, beta=1, gamma=g, pm=1)
}
### <---------------------------------------------------------------------->
#' @rdname dstablecnt
pstablecnt <- function(x, alpha=NULL, nu0=0, theta=1, lambda=NULL) {
    if (is.null(alpha) & !is.null(lambda)) alpha <- 2/lambda
    .stablecnt.check(alpha, nu0, theta)
    if (alpha==1/2) {
        x0 <- ifelse(x >= nu0, x-nu0, NaN)
        return(stats::pgamma(x0, shape=3/2, scale=4*theta))
    }

    cdf1 <- function(x) {
        fn <- function(x) dstablecnt(x, alpha=alpha, nu0=nu0, theta=theta, lambda=lambda)
        integrate(fn, lower=0.001+nu0, upper=x)$value
    }
    sapply(x, cdf1)
}
### <---------------------------------------------------------------------->
#' @rdname dstablecnt
rstablecnt <- function(n, alpha=NULL, nu0=0, theta=1, lambda=NULL) {
    if (is.null(alpha) & !is.null(lambda)) alpha <- 2/lambda
    .stablecnt.check(alpha, nu0, theta)
    if (alpha==1/2) {
        return(nu0 + rgamma(n, shape=3/2, scale=4*theta))
    }
    u <- runif(n, 0.001, 1-0.001)
    qstablecnt(u, alpha=alpha, nu0=nu0, theta=theta, lambda=lambda)
}
### <---------------------------------------------------------------------->
#' @rdname dstablecnt
qstablecnt <- function(q, alpha=NULL, nu0=0, theta=1, lambda=NULL) {
    if (is.null(alpha) & !is.null(lambda)) alpha <- 2/lambda
    .stablecnt.check(alpha, nu0, theta)
    if (alpha==1/2) {
        return(nu0 + qgamma(q, shape=3/2, scale=4*theta))
    }
    qnt1 <- function(q) {
        if (q < 0.001) return(nu0)
        if (q > 1-0.001) return(Inf)
        fn <- function(x) pstablecnt(x, alpha=alpha, nu0=nu0, theta=theta, lambda=lambda) - q
        uniroot(fn, lower=0.001+nu0, upper=nu0+100*theta)$root
    }
    simplify2array(parallel::mclapply(q, qnt1))
}
### <---------------------------------------------------------------------->
#' @rdname dstablecnt
cfstablecnt <- function(s, alpha=NULL, nu0=0, theta=1, lambda=NULL) {
    if (is.null(alpha) & !is.null(lambda)) alpha <- 2/lambda
    .stablecnt.check(alpha, nu0, theta)
    if (alpha==1/2) {
        return(exp(1i*s*nu0) * (1-4i*s*theta)^(-3/2))
    }
    stop(paste("ERROR: cfstablecnt is not supported for alpha:", alpha))
}
### <---------------------------------------------------------------------->
#' @rdname dstablecnt
kstablecnt <- function(alpha=NULL, nu0=0, theta=1, lambda=NULL) {
    if (is.null(alpha) & !is.null(lambda)) alpha <- 2/lambda
    .stablecnt.check(alpha, nu0, theta)
    if (alpha==1/2) {
        kN <- function(i) {
            if (length(i)>1) return(sapply(i,kN))
            if (i==1) return(nu0+6*theta)
            if (i==2) return(24*theta^2)
            if (i==3) return(192*theta^3)
            if (i==4) return(2304*theta^4)
            stop
        }
        skw <- kN(3)/kN(2)^(3/2)
        kur <- kN(4)/kN(2)^2+3

        k <- c(kN(1:4), skw, kur)
        names(k) <- c("mean", "var", "k3", "k4", "skewness", "kurtosis")
        return(k)
    }
    stop(paste("ERROR: rstablecnt is not supported for alpha:", alpha))
}
### <---------------------------------------------------------------------->
.stablecnt.check <- function(alpha, nu0, theta) {
    if (length(alpha)!=1) stop("Either lambda or alpha must be provided as length-one numeric")
    if (length(theta)!=1) stop("theta must be length-one")
    if (alpha>=1 | alpha<=0) stop("alpha is only supported between 0 and 1")
    if (nu0<0) stop("nu0 must be positive")
    if (theta<0) stop("theta must be positive")
}
### <---------------------------------------------------------------------->
