#' The moments and MGF of ecld
#' 
#' Compute the moments and MGF of ecld for \code{mu=0} (centered),
#' via analytical result whenever is available.
#' SGED is supported.
#'
#' @param object an object of ecd class
#' @param order numeric, order of the moment to be computed
#' @param t numeric, for MGF
#' @param ignore.mu logical, disgard mu; otherwise, stop if mu is not zero.
#'
#' @return numeric
#'
#' @keywords moment
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.moment
#' @export ecld.mgf
#' @export ecld.mgf_by_sum
#' @export ecld.mgf_quartic
#'
#' @examples
#' ld <- ecld(lambda=3, sigma=0.01*ecd.mp1)
#' ecld.moment(ld, 2)
#' ecld.mgf(ld)
### <======================================================================>
"ecld.moment" <- function(object, order, ignore.mu=TRUE)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for consistent type
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta * one
    n <- order
    
    if (! ignore.mu) {
        if (object@mu!=0) {
            stop("This function only handles moments with mu=0.")
        }
    }

    # SGED
    if (object@is.sged) {
        sn = s*(1-b)
        sp = s*(1+b)
        m1 <- sp^(n+1) + (-one)^n*sn^(n+1)
        m2 <- gamma(lambda*(n+1)/2)/s/2/gamma(lambda/2)
        return(m1*m2)
    }

	# symmetric
    if (object@beta==0) {
        if (floor(n/2)!=n/2) return(0) # odd moments are zero
        x <- gamma(lambda*(n+1)/2)
        y <- gamma(lambda/2)
        return(s^n*x/y)
    }

	# asymmetric
    if (lambda==2) {
    	B0 <- ecld.laplace_B(b, 0)
    	Bp <- ecld.laplace_B(b, 1)
    	Bn <- ecld.laplace_B(b, -1)
        x <- gamma(n+1)/(2*B0)
        y <- Bp^(n+1) + (-one)^n*Bn^(n+1)
        return(s^n*x*y)
    }
    
    if (lambda > 2) {
        # MPFR is channelled through sigma=1
        # since we are using unit distribution, either way should be fine
        ld0 <- ecld(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one)
        d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
        e_y2 <- function(x) {
            x^n * (exp(ecld.solve(ld0,x)) + (-one)^n*exp(ecld.solve(ld0,-x)))
        }
        M <- ecd.integrate(d0, e_y2, 0, Inf)
        if (M$message != "OK") {
            stop("Failed to integrate moment from unit distribution")
        }
        C <- ecld.const(ld0)
        return(ecd.mpnum(object, s^n/C*M$value))
    }
    
    stop("Unknown analytic formula for moment")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.moment
"ecld.mgf" <- function(object, t=1)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for consistent type

    lambda <- object@lambda
    s <- object@sigma * one
    b <- object@beta
    mu <- object@mu * one
    
    # SGED
    if (object@is.sged) {
        if (lambda==2) {
            y <- 1 - 2*b*s*t - (1-b^2)*s^2*t^2
            if (y<0) {
                stop("SGED lambda=2: MGF can not converge!")
            }
            return(exp(t*mu)/y)
        }
        
        # summation
        nmax <- floor(ecd.mp2f(ecld.mgf_trunc(object)))
        if (is.na(nmax)) {
            stop("NA found in mgf_trunc! Is sigma too large?")
        }
        
        if (nmax > 5000) {
            nmax <- 5000 # cap the length of summation
        }
        
        n <- seq(1, nmax) # use all terms, no skip
        terms <- ecld.mgf_term(object, n)
        if (any(is.na(terms))) {
            stop("NA found in mgf_term. Consider using MPFR to improve precision!")
        }
        return(1+sum(terms))
    }

    # normal
    if (lambda==1) {
        if (b != 0) {
            stop("lambda=1: beta must be zero")
        }
        return(exp(t*mu + s^2*t^2/4))
    }
    if (lambda==2) {
        y <- 1 - b*s*t - s^2*t^2
        if (y<0) {
            stop("lambda=2: MGF can not converge!")
        }
        return(exp(t*mu)/y)
    }
    #if (lambda==4 & b==0 & s*t > 0.005) {
    #    return(ecld.mgf_quartic(object,t=t))
    #}
    
    # use truncated summation for symmetric
    if (lambda>2 & b==0) {
        return(ecld.mgf_by_sum(object,t=t))
    }
    
    stop("Unknown analytic formula for MGF")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.moment
"ecld.mgf_by_sum" <- function(object, t=1)
{
    nmax <- ecd.mp2f(ecld.mgf_trunc(object, t=t))
    if ((!is.numeric(nmax)) | is.na(nmax)) {
        stop(paste("nmax is not valid numeric: ", nmax))
    }
    if (nmax > 5000) {
        nmax <- 5000 # cap the length of summation
    }
    #if (object@lambda==4 & object@sigma*t < 0.01) {
    #    # quartic dist only needs a few terms for small sigma
    #    z = 1/(4*object@sigma*t)
    #    N1 = z^2-1 # truncation rule
    #    N2 = 10/log10(z^2) # 10-digit precision
    #    nmax = max( c(min(c(N1,N2)), 12)) # min of precision, but at least 12 sums
    #
    
    n <- seq(2, floor(nmax), by=2) # skip odd terms
    terms <- ecld.mgf_term(object, order=n, t=t)
    if (any(is.na(terms))) {
        stop("NA found in mgf_term. Consider using MPFR to improve precision!")
    }
    return(1+sum(terms))
}
### <---------------------------------------------------------------------->
#' @rdname ecld.moment
"ecld.mgf_quartic" <- function(object, t=1)
{
    ecld.validate(object, sged.allowed=FALSE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for consistent type
    
    lambda <- object@lambda
    s <- object@sigma * one
    b <- object@beta
    mu <- object@mu * one

    st <- s*t
    z = 1/sqrt(4*s*t)
    dz = z^2*exp(-z^2)
    Mz <- z^3*(ecd.erfq(z,-1) -ecd.erfq(z,1))
    return(exp(t*mu) * (Mz+dz))
}
### <---------------------------------------------------------------------->
