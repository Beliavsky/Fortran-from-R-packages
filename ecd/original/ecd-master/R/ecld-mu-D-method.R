#' mu_D of ecld
#'
#' The analytic solutions for risk-neutral drift. If analytic form doesn't
#' exist, it uses integral of unit distribution. This is different from
#' \code{ecld.mgf} where series summation is used.
#'
#' @param object an object of ecld class
#' @param validate logical, if true (default), stop when the result is NaN or infinite.
#'
#' @return numeric
#'
#' @keywords option-pricing
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.mu_D
#' @export ecld.mu_D_quartic
#' @export ecld.mu_D_by_sum
#' @export ecld.mu_D_integrate
#'
#' @examples
#' ld <- ecld(sigma=0.01*ecd.mp1)
#' ecld.mu_D(ld)
### <======================================================================>
"ecld.mu_D" <- function(object, validate=TRUE)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
	
    # validate
    pass <- function(x) {
        if (validate) {
            if (is.na(x)) {
                stop(paste("mu_D is NaN:",
                     "lambda=", ecd.mp2f(lambda), "beta=", ecd.mp2f(b), "sigma=", ecd.mp2f(s)))
            }
            if (abs(x)==Inf) {
                stop(paste("mu_D is infinite:",
                     "lambda=", ecd.mp2f(lambda), "beta=", ecd.mp2f(b), "sigma=", ecd.mp2f(s)))
            }
        }
        return(x)
    }
    
    # SGED
    if (object@is.sged) {
        # MPFR is channelled through sigma=1
	    # since we are using unit distribution, either way should be fine
	    d0 <- ecd(lambda=ecd.mp2f(lambda), sigma=one, bare.bone=TRUE)
        
        sn <- s*(1-b)
        sp <- s*(1+b)
        
        # sp
        e_y1 <- function(xi) exp(-xi^(2/lambda))*(exp(sp*xi)-1)
        xt <- ecld.y_slope_trunc(object)/sp
        if (xt > .ecd.mpfr.N.sigma) xt <- .ecd.mpfr.N.sigma
	    M1 <- ecd.integrate(d0, e_y1, 0, xt)
        
        # sn
        e_y2 <- function(xi) exp(-xi^(2/lambda))*(exp(-sn*xi)-1)
        M2 <- ecd.integrate(d0, e_y2, 0, Inf)
        
        if (M1$message != "OK") {
	        stop("Failed to integrate right tail from unit distribution")
	    }
        if (M2$message != "OK") {
	        stop("Failed to integrate left tail from unit distribution")
	    }
	    C <- ecld.const(object)
        
	    return(pass(ecd.mpnum(object, -log(1 + (sp*M1$value + sn*M2$value)/C))))
    }

    # normal
    if (lambda==1) {
        if (b != 0) {
            stop("lambda=1: beta must be zero")
        }
        return(pass(-s^2/4))
    }
    
    if (lambda==2) {
        return(pass(log(1-b*s-s^2)))
    }
    
    if (lambda==4 & b==0) {
        object@mu <- 0
        return(ecld.mu_D_by_sum(object))
    }

    # must integrate to get M(1)
	if (TRUE) {
        return(ecld.mu_D_integrate(object, validate=validate))
	}
	
    stop("Unknown analytic formula for mu_D")

}
### <======================================================================>
#' @rdname ecld.mu_D
"ecld.mu_D_quartic" <- function(object)
{
    s <- object@sigma
    z = 1/2/sqrt(s)
    dz = z^2*exp(-z^2)
    Mz <- z^3*(ecd.erfq(z,-1) -ecd.erfq(z,1))
    return(-log(Mz + dz))
}
### <======================================================================>
#' @rdname ecld.mu_D
"ecld.mu_D_by_sum" <- function(object)
{
    return(-log(ecld.mgf_by_sum(object))+object@mu)
}
### <======================================================================>
#' @rdname ecld.mu_D
"ecld.mu_D_integrate" <- function(object, validate=TRUE)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
    
    # validate
    pass <- function(x) {
        if (validate) {
            if (is.na(x)) {
                stop(paste("mu_D is NaN:",
                "lambda=", ecd.mp2f(lambda), "beta=", ecd.mp2f(b), "sigma=", ecd.mp2f(s)))
            }
            if (abs(x)==Inf) {
                stop(paste("mu_D is infinite:",
                "lambda=", ecd.mp2f(lambda), "beta=", ecd.mp2f(b), "sigma=", ecd.mp2f(s)))
            }
        }
        return(x)
    }
    
    # MPFR is channelled through sigma=1
    # since we are using unit distribution, either way should be fine
    ld0 <- ecld(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one)
    d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
    
    e_y1 <- function(xi) exp(ecld.solve(ld0,xi))*(exp(s*xi)-1)
    xt <- ecld.y_slope_trunc(object)/s
    if (xt > .ecd.mpfr.N.sigma) xt <- .ecd.mpfr.N.sigma
    M1 <- ecd.integrate(d0, e_y1, 0, xt)
    
    e_y2 <- function(xi) exp(ecld.solve(ld0,-xi))*(exp(-s*xi)-1)
    M2 <- ecd.integrate(d0, e_y2, 0, Inf)
    
    if (M1$message != "OK") {
        stop("Failed to integrate right tail from unit distribution")
    }
    if (M2$message != "OK") {
        stop("Failed to integrate left tail from unit distribution")
    }
    C <- ecld.const(object)/s
    
    return(pass(ecd.mpnum(object, -log(1+(M1$value+M2$value)/C))))

}
### <---------------------------------------------------------------------->
