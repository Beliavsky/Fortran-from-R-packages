#' Stable lambda distribution
#'
#' Implements the stable lambda distribution (SLD) and the quartic stable lambda distribution (QSLD).
#' Be aware of the performance concerns:
#' (a) The cumulative density function is implemented by direct integration over the density.
#' (b) The quantile function is implemented by root finding on cumulative density function.
#'
#' @param n numeric, number of observations.
#' @param x numeric, vector of responses.
#' @param q numeric, vector of quantiles.
#' @param s numeric, vector of responses for characteristic function.
#' @param t numeric, the time parameter, where the variance is t, default is 1.
#' @param nu0 numeric, the location parameter, default is 0.
#' @param theta numeric, the scale parameter, default is 1.
#' @param convo numeric, the convolution number, default is 1.
#' @param beta.a numeric, the skewness parameter, default is 0. This number is annualized by sqrt(t).
#' @param mu numeric, the location parameter, default is 0.
#' @param lambda numeric, the shape parameter, default is 4.
#' @param method character, method of characteristic function (CF) calculation. Default is "a".
#'               Method a uses cfstdlap x dstablecnt.
#'               Method b uses dstdlap x cfstablecnt.
#'               Method c uses direct integration on PDF up to 50 stdev.
#'               They should yield the same result.
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
#' @keywords stable-Lambda
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats rexp
#' @importFrom stats rnorm
#'
#' @export rqsl
#' @export dqsl
#' @export pqsl
#' @export qqsl
#' @export cfqsl
#' @export kqsl
#' @export rsl
#' @export dsl
#' @export psl
#' @export qsl
#' @export cfsl
#' @export ksl
#'
#' @section Details:
#'   The stable lambda distribution is the stationary distribution for financial asset returns.
#'   It is a product of the stable count distribution and the Lihn-Laplace process.
#'   The density function is defined as \deqn{
#'     P_{\Lambda}^{\left(m\right)}\left(x;t,\nu_{0},\theta,\beta_{a},\mu\right)
#'     \equiv\int_{\nu_{0}}^{\infty} \frac{1}{\nu}\,
#'     f_{\mathit{L}}^{\left(m\right)}\left(\frac{x-\mu}{\nu};t,\beta_{a}\sqrt{t}\right)\,
#'     \mathit{N}_\alpha\left(\nu;\nu_{0},\theta\right)d\nu
#'   }{
#'     P_\Lambda^(m) (x;t,\nu_0,\theta,\beta_a,\mu)
#'     = integrate_\nu_0^\infty 1/\nu
#'     f_L^(m)((x-\mu)/\nu; t,\beta_a sqrt(t))
#'     N_\alpha(\nu; \nu_0,\theta) d\nu
#'   }
#'   where
#'   \eqn{f_L^{(m)}(.)} is the Lihn-Laplace distribution and
#'   \eqn{N_\alpha(.)} is the quartic stable count distribution.
#'   \eqn{t} is the time or sampling period,
#'   \eqn{\alpha} is the stability index which is \eqn{2/\lambda},
#'   \eqn{\nu_0} is the floor volatility parameter,
#'   \eqn{\theta} is the volatility scale parameter,
#'   \eqn{\beta_a} is the annualized asymmetric parameter,
#'   \eqn{\mu} is the location parameter.
#'   \cr
#'   The quartic stable lambda distribution (QSLD) is a specialized distribution with \eqn{\lambda=4} aka \eqn{\alpha=0.5}.
#'   The PDF integrand has closed form, and all the moments have closed forms.
#'   Many financial asset returns can be fitted by QSLD precisely up to 4 standard deviations.
#'
#' @seealso
#'   \code{\link{dstablecnt}} for \eqn{N_\alpha(.)},
#'   and \code{\link{dstdlap}} for \eqn{f_L^{(m)}(.)}.
#'
#' @references
#'   For more detail, see Section 8.2 of
#'   Stephen Lihn (2017).
#'   \emph{A Theory of Asset Return and Volatility
#'   under Stable Law and Stable Lambda Distribution}.
#'   SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
#'
#' @examples
#'   # generate the quartic pdf for SPX 1-day distribution
#'   x <- c(-0.1, 0.1, by=0.001)
#'   pdf <- dqsl(x, t=1/250, nu0=6.92/100, theta=1.17/100, convo=2, beta=-1.31)
#'
### <======================================================================>
dsl <- function(x, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    integrand_xp <- function(x, nu) {
        L <- dstdlap((x-mu)/(nu+nu0), t=t, convo=convo, beta=beta.a*sqrt(t))
        N <- dstablecnt(nu, alpha=2/lambda, nu0=0, theta=theta)
        v <- N*L/(nu+nu0)
        return(ifelse(is.na(v), 0, v)) # bypass the singular points
    }

    pdf <- function(x) {
        f <- function(nu) integrand_xp(x, nu)
        S <- sqrt(t*((nu0+6*theta)^2 + 24*theta^2)) # scale

        Va <- integrate(f, lower=0.000001*S, upper=0.001*S)$value # very small and problematic
        V0 <- integrate(f, lower=0.001*S, upper=0.1*S)$value
        V1 <- integrate(f, lower=0.1*S, upper=1*S)$value
        V2 <- integrate(f, lower=1*S,   upper=10*S)$value
        V3 <- integrate(f, lower=10*S,  upper=100*S)$value
        V4 <- integrate(f, lower=100*S, upper=1000*S)$value
        return(Va+V0+V1+V2+V3+V4)
    }

    ld <- ecld(4)
    V <- ecd.mcsapply(ld, x, pdf)
    return(V)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
dqsl <- function(x, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
rqsl <- function(n, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    rsl(n, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
rsl <- function(n, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    L <- rstdlap(n, t=t, convo=convo, beta=beta.a*sqrt(t))
    N <- rstablecnt(n, alpha=2/lambda, nu0=nu0, theta=theta)
    return(L*N+mu)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
pqsl <- function(x, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    psl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
psl <- function(x, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    if (length(x)>1) {
        ld <- ecld(4)
        g <- function(x) psl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=lambda)
        v <- ecd.mcsapply(ld, x, g)
        return(v)
    }

    S <- sqrt(t*((nu0+6*theta)^2 + 24*theta^2)) # scale
    xmax = 50*S
    f <- function(x) dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=lambda)
    if (x <= mu) {
        if (x < -xmax) return(0)
        v <- integrate(f, lower=-xmax, upper=x)
        return(v$value)
    }
    v1 <- integrate(f, lower=-xmax, upper=mu)
    v2 <- integrate(f, lower=mu, upper=x)
    return(v1$value + v2$value)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
qqsl <- function(q, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    qsl(q, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
qsl <- function(q, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    if (length(q)>1) {
        ld <- ecld(2)
        g <- function(q) qsl(q, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=lambda)
        v <- ecd.mcsapply(ld, q, g)
        return(v)
    }
    if (q > 1 | q < 0) stop("quantile out of range")
    if (q == 1) return(Inf)
    if (q == 0) return(-Inf)

    f <- function(x) psl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=lambda)-q

    # check out of bound condition
    S <- sqrt(t*((nu0+6*theta)^2 + 24*theta^2)) # scale
    xmax = 50*S
    if (f(mu-xmax) >= 0) return(-Inf)
    if (f(mu+xmax) <= 0) return(Inf)

    r <- uniroot(f, lower=mu-xmax, upper=mu+xmax)
    return(r$root)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
kqsl <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    ksl(t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
ksl <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    kL <- kstdlap(t=t, convo=convo, beta=beta.a*sqrt(t))
    kN <- kstablecnt(alpha=2/lambda, nu0=nu0, theta=theta)
    mC <- k2mnt(kN)*k2mnt(kL)
    kC <- mnt2k(mC)
    skw <- kC[3]/kC[2]^(3/2)
    kur <- kC[4]/kC[2]^2 + 3

    kC[1] <- kC[1]+mu
    k <- c(kC[1:4], skw, kur)
    names(k) <- c("mean", "var", "k3", "k4", "skewness", "kurtosis")
    return(k)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
cfqsl <- function(s, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, method="a") {
    cfsl(s, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4, method=method)
}
### <---------------------------------------------------------------------->
#' @rdname dsl
cfsl <- function(s, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4, method="a") {
    # TODO the upper bounds of the integrals are all based on quartic distribution
    # User needs to be careful when extending to other lambda's
    # Future development!!!
    if (method=="a") {
        N <- 100000
        nu <- nu0 + seq(0, 100, length.out=N)*theta
        M <- exp(1i*mu*s)
        L <- cfstdlap(s*nu, t=t, convo=convo, beta=beta.a*sqrt(t))
        N <- dstablecnt(nu, alpha=2/lambda, nu0=nu0, theta=theta)
        dnu <- nu[2]-nu[1]
        return(M*sum(L*N)*dnu)
    }
    if (method=="b") {
        N <- 100000
        x <- seq(-50, 50, length.out=N+1)*sqrt(t)
        dx <- x[2]-x[1]
        M <- exp(1i*mu*s)
        L <- dstdlap(x, t=t, convo=convo, beta=beta.a*sqrt(t))
        N <- cfstablecnt(s*x, alpha=2/lambda, nu0=nu0, theta=theta)
        return(M*sum(L*N)*dx)
    }
    if (method=="c") {
        # TODO this is a hack for other lambda, but our primary use case is lambda=4
        sd <- qsl_variance_analytic(t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a)^0.5
        M <- exp(1i*mu*s)
        p_re <- function(x) cos(s*x)*dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a)
        p_im <- function(x) sin(s*x)*dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a)
        intg <- function(f,a,b) integrate(f, lower=sd*a, upper=sd*b)$value
        cf_re <- c(intg(p_re, 0, 5), intg(p_re, 5, 10), intg(p_re, 10, 50),
                   intg(p_re, -5, 0), intg(p_re, -10, -5), intg(p_re, -50, -10))

        cf_im <- c(intg(p_im, 0, 5), intg(p_im, 5, 10), intg(p_im, 10, 50),
                   intg(p_im, -5, 0), intg(p_im, -10, -5), intg(p_im, -50, -10))

        # print(paste(sd, cf_re, cf_im))
        return(M*(sum(cf_re) + 1i*sum(cf_im)))
    }

    stop(paste("cfsl: method is not supported", method))
}
### <---------------------------------------------------------------------->
