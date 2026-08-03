#' Standardized Laplace process and distribution
#'
#' Implements the standardized Laplace process and distribution.
#' Be aware of the performance concerns:
#' (a) The cumulative density function is implemented by direct integration over the density.
#' (b) The quantile function is implemented by root finding on cumulative density function.
#'
#' @param n numeric, number of observations.
#' @param x numeric, vector of responses.
#' @param q numeric, vector of quantiles.
#' @param s numeric, vector of responses for characteristic function.
#' @param t numeric, the time parameter, of which the variance is t.
#' @param mu numeric, location parameter, default is 0.
#' @param beta numeric, skewness parameter according to skewed lambda distribution, default is 0.
#' @param convo numeric, the convolution number, default is 1, which is Laplace without convolution.
#'              There is a special provision in \code{rstdlap}, where it will simulate the Wiener process
#'              if \code{convo=Inf} and \code{beta=0}.
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
#' @keywords Laplace
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats runif
#' @importFrom stats rexp
#' @importFrom stats rnorm
#'
#' @export dstdlap
#' @export pstdlap
#' @export qstdlap
#' @export rstdlap
#' @export cfstdlap
#' @export kstdlap
#' @export dstdlap_poly
#'
#' @section Details:
#'   The Lihn-Laplace distribution is the stationary distribution of Lihn-Laplace process.
#'   The density function is defined as \deqn{
#'     f_{\mathit{L}}^{\left(m\right)}\left(x;t,\beta,\mu\right)
#'     \equiv\frac{1}{\sqrt{\pi}\Gamma\left(m\right)\sigma_{m}}\,
#'     \left|\frac{x-\mu}{2B_{0}\sigma_{m}}\right|^{m-\frac{1}{2}}
#'     K_{m-\frac{1}{2}}\left(\left|\frac{B_{0}(x-\mu)}{\sigma_{m}}\right|\right)
#'     e^{\frac{\beta (x-\mu)}{2\sigma_{m}}}
#'   }{
#'     f_L^(m)(x;t,\beta,\mu) =
#'     1/(sqrt(\pi)\Gamma(m)\sigma_m)
#'     |(x-\mu)/2B0 \sigma_m|^(m-1/2)
#'     K_(m-1/2)(|B0 (x-\mu)/\sigma_m|)
#'     exp(\beta (x-\mu)/2\sigma_m)
#'   }
#'   where \deqn{
#'     \sigma_{m}\equiv\sqrt{\frac{t}{m\left(2+\beta^{2}\right)}},
#'     \:
#'     B_{0}\equiv\sqrt{1+\frac{1}{4}\beta^{2}}.
#'   }{
#'     \sigma_m = sqrt(t/m/(2+\beta^2)),
#'     B0 = sqrt(1+\beta^2/4).
#'   }
#'   \eqn{K_n(x)} is the modified Bessel function of the second kind.
#'   \eqn{t} is the time or sampling period,
#'   \eqn{\beta} is the asymmetric parameter,
#'   \eqn{\mu} is the location parameter.
#'
#' @references
#'   For more detail, see Section 5.4 of
#'   Stephen Lihn (2017).
#'   \emph{A Theory of Asset Return and Volatility
#'   under Stable Law and Stable Lambda Distribution}.
#'   SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
#'
#' @examples
#'   # generate the pdf at time t=1 for the second convolution and beta = 0.1 for skewness
#'   x <- c(-10, 10, by=0.1)
#'   pdf <- dstdlap(x, t=1, convo=2, beta=0.1)
#'
### <======================================================================>
dstdlap <- function(x, t=1, convo=1, beta=0, mu=0) {
    if (convo<=0) stop(paste("ERROR: Lihn-Laplace pdf is not supported for convolution:", convo))

    x <- x-mu # take care of location
    m <- convo
    Sm <- sqrt(t/(2+beta^2)/m)
    B0 <- sqrt(1+beta^2/4)
    w <- abs(B0*x/Sm)

    C <- 1/gamma(m)/sqrt(pi)/Sm
    K <- (w/2)^(m-1/2) * besselK(w,m-1/2)
    C0 <- if (convo > 1/2) C*gamma(m-1/2)/2 else NaN
    CK <- ifelse(w==0, C0, C*K) # K(0) is Inf
    Eb <- (B0^-2)^(m-1/2) * exp(beta*x/2/Sm)
    return(CK*Eb)
}
### <---------------------------------------------------------------------->
#' @rdname dstdlap
pstdlap <- function(x, t=1, convo=1, beta=0, mu=0) {
    if (length(x)>1) {
        ld <- ecld(2)
        g <- function(x) pstdlap(x, t=t, convo=convo, beta=beta, mu=mu)
        v <- ecd.mcsapply(ld, x, g)
        return(v)
    }

    f <- function(x) dstdlap(x, t=t, convo=convo, beta=beta, mu=mu)
    if (x <= mu) {
        v <- integrate(f, lower=-Inf, upper=x)
        return(v$value)
    }
    v1 <- integrate(f, lower=-Inf, upper=mu)
    v2 <- integrate(f, lower=mu, upper=x)
    return(v1$value + v2$value)
}
### <---------------------------------------------------------------------->
#' @rdname dstdlap
qstdlap <- function(q, t=1, convo=1, beta=0, mu=0) {
    if (length(q)>1) {
        ld <- ecld(2)
        g <- function(q) qstdlap(q, t=t, convo=convo, beta=beta, mu=mu)
        v <- ecd.mcsapply(ld, q, g)
        return(v)
    }
    if (q > 1 | q < 0) stop("quantile out of range")
    if (q == 1) return(Inf)
    if (q == 0) return(-Inf)

    f <- function(x) pstdlap(x, t=t, convo=convo, beta=beta, mu=mu)-q

    # check out of bound condition
    xmax = 20*sqrt(t)
    if (f(mu-xmax) >= 0) return(-Inf)
    if (f(mu+xmax) <= 0) return(Inf)

    r <- uniroot(f, lower=mu-xmax, upper=mu+xmax)
    return(r$root)
}
### <---------------------------------------------------------------------->
#' @rdname dstdlap
rstdlap <- function(n, t=1, convo=1, beta=0, mu=0) {
    if(convo==Inf & beta==0) {
        # just Wiener process
        return(rnorm(n, mean=mu, sd=sqrt(t)))
    }
    if (convo==1) {
        if (beta==0) {
            L <- stats::rnorm(n)*sqrt(stats::rexp(n,1))*sqrt(t)
            return(L+mu)
        }
        # skewed Laplace implementation
        B0 <- sqrt(1+beta^2/4)
        B1 <- B0+beta/2 # B plus
        B2 <- B0-beta/2 # B minus
        W <- B2/B0/2
        sgn <- sign(stats::runif(n)-W)
        sigma <- ifelse(sgn > 0, B1, B2)
        L <- stats::rexp(n, 1/sigma)*sgn
        sd <- sqrt(2+beta^2)
        return(L/sd*sqrt(t)+mu)
    }
    # compose convolution recursively
    if (convo < 1) stop("fractional convolution is not supported")

    if (convo > 1 & convo < 2 & beta==0) {
        # we use bivariate LLP's effective convolution to produce the result
        a <- 0.3374
        meff <- function(rho) (1+a*rho^2+(1-a)*rho^4) - convo
        rho <- uniroot(meff, lower=0, upper=1)$root

        c <- sign(rho)*sqrt(1/2*(1-sqrt(1-rho^2)))
        SS <- matrix(c(sqrt(1-c^2), c, c, sqrt(1-c^2)), nrow=2, ncol=2, byrow=TRUE)
        L <- matrix(NaN, nrow=2, ncol=n, byrow=TRUE)
        L[1,] <- rstdlap(n, t=t, convo=1)
        L[2,] <- rstdlap(n, t=t, convo=1)
        for (i in 1:n) L[,i] <- SS %*% L[,i]

        return(L[1,])
    }

    if (convo!=ceiling(convo)) warning("real-number convolution is only experimental")

    L <- rep(0,n)
    frac <- convo+1-ceiling(convo)
    if (frac==0) frac <- 1
    convo1 <- ceiling(convo)
    for (i in 1:convo1) {
        f <- if (i==convo1) frac else 1
        L <- L + rstdlap(n, t*f, beta=beta)
    }
    return(L/sqrt(convo)+mu)
}
### <---------------------------------------------------------------------->
#' @rdname dstdlap
cfstdlap <- function(s, t=1, convo=1, beta=0, mu=0) {

    m <- convo
    Sm <- sqrt(t/(2+beta^2)/m)

    cf_mu <- exp(1i*mu*s)
    cf1 <- 1 - 1i*beta*s*Sm + s^2*Sm^2
    return(cf_mu*cf1^(-m))
}
### <---------------------------------------------------------------------->
#' @rdname dstdlap
kstdlap <- function(t=1, convo=1, beta=0, mu=0) {
    kL <- function(i) {
        if (length(i)>1) return(sapply(i,kL))
        B <- 2+beta^2
        if (i==1) return(mu+beta*sqrt(convo*t/B))
        if (i==2) return(t)
        if (i==3) return(2*t^(3/2)/sqrt(convo)*beta*(3+beta^2)/B^(3/2))
        if (i==4) return(6*t^2/convo*(2+4*beta^2+beta^4)/B^2)
        stop
    }
    skw <- kL(3)/kL(2)^(3/2)
    kur <- kL(4)/kL(2)^2 + 3

    k <- c(kL(1:4), skw, kur)
    names(k) <- c("mean", "var", "k3", "k4", "skewness", "kurtosis")
    return(k)
}
### <---------------------------------------------------------------------->
#' @rdname dstdlap
dstdlap_poly <- function(x, t=1, convo=1, beta=0, mu=0) {
    x <- x-mu # take care of location
    m <- convo
    Sm <- sqrt(t/(2+beta^2)/m)
    B0 <- sqrt(1+beta^2/4)
    B1 <- B0+beta/2
    B2 <- B0-beta/2
    w <- abs(B0*x/Sm)

    B <- ifelse(x>0, B1, B2)
    Ey <- exp(-abs(x/Sm/B))

    if (convo==1) return(Ey/Sm/(2*B0))
    if (convo==2) return(Ey/Sm*(w+1)/(4*B0^3))
    if (convo==3) {
        p <- w^2+3*w+3
        return(Ey/Sm*p/(16*B0^5))
    }
    if (convo==4) {
        p <- w^3+6*w^2+15*w+15
        return(Ey/Sm*p/(96*B0^7))
    }

    stop(paste("ERROR: Lihn-Laplace pdf is not supported for convolution:", convo))
}
### <---------------------------------------------------------------------->
