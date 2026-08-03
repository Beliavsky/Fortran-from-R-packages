#' Constructor of ecld class
#' 
#' Construct an \code{\link{ecld-class}} by providing the required parameters.
#' The default is the standard symmetric cusp distribution.
#' The default also doesn't calculate any ecd extension.
#' \code{ecld.from} allows you to pass the parameters from an existing ecd object.
#' \code{ecld.validate} checks if an object is ecld class.
#' \code{ecld.quartic} is a convenient constructor designed for quartic distribution.
#' \code{ecld.from_sd} calculates sigma from a given sd and renders a vanila ecld object.
#'
#' @param lambda numeric, the lambda parameter. Must be positive. Default: 3.
#' @param sigma numeric, the scale parameter. Must be positive. Default: 1.
#' @param beta  numeric, the skewness parameter. Default: 0.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param with.ecd logical, also calculate the ecd object, default is \code{FALSE}.
#' @param with.mu_D logical, also calculate the ecd risk-neutral drift, default is \code{FALSE}.
#'                  If \code{TRUE}, this flag supercedes \code{with.ecd}.
#'                  Also \code{mu} must set to zero.
#' @param with.RN logical, also calculate the risk-neutral ecd object, default is \code{FALSE}.
#'                If \code{TRUE}, this flag supercedes \code{with.mu_D}.
#' @param object an object of ecld class
#' @param is.sged logical, if \code{TRUE}, interpret parameters as SGED.
#' @param verbose logical, display timing information, for debugging purpose, default is \code{FALSE}.
#' @param sged.allowed logical, used in \code{ecld.validate} to indicate if the function allows SGED.
#' @param sged.only logical, used in \code{ecld.validate} to indicate if the function is only for SGED.
#' @param mu_plus,mu_plus_ratio numeric, excess value in addition to \code{mu_D}.
#'                When ratio is provided, it is relative to the stdev.
#' @param epsilon The supplemental residual premium for lambda transformation.
#'                It is default to NaN in ecld constructor since its meaning is not defined.
#' @param rho The supplemental momentum shift for lambda transformation.
#'                It is default to NaN in ecld constructor since its meaning is not defined.
#' @param sd numeric, the scale parameter expressed in stdev instead of sigma. Internally,
#'                It is converted to sigma via \code{uniroot} on \code{ecld.sd}.
#'                Must be positive. Default: 1.
#'
#' @return an object of ecld class
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld
#' @export ecld.from
#' @export ecld.quartic
#' @export ecld.validate
#' @export ecld.from_sd
#'
#' @examples
#' ld <- ecld()
#' ld <- ecld(2, 0.01)
#' ld <- ecld.from_sd(3, 0.1)

### <======================================================================>
"ecld" <- function(lambda = 3, sigma = 1, beta = 0, mu = 0,
                   epsilon = NaN, rho = NaN,
                   with.ecd = FALSE, with.mu_D = FALSE,
                   with.RN = FALSE, is.sged = FALSE, verbose=FALSE)
{
    call <- match.call()

    if(!is.numericMpfr(lambda)){
      stop("Parameter 'lambda' must be numericMpfr!\n")
    }
    if(!is.numericMpfr(sigma)){
        stop("Parameter 'sigma' must be numericMpfr!\n")
    }
    if(!is.numericMpfr(beta)){
        stop("Parameter 'beta' must be numericMpfr!\n")
    }
    if(!is.numericMpfr(mu)){
        stop("Parameter 'mu' must be numericMpfr!\n")
    }

    # -------------
    if(sigma <= 0){
        stop("Parameter 'sigma' must be positive!\n")
    }
    if(lambda <= 0){
        stop("Parameter 'lambda' must be positive!\n")
    }

    sum <- lambda+sigma+beta+mu
    if(!(abs(sum)>=0 & abs(sum) != Inf)) {
        stop(paste("Parameters must be finite, known real numbers!",
             "lambda=", ecd.mp2f(lambda), "sigma=", ecd.mp2f(sigma),
             "beta=", ecd.mp2f(beta), "mu=", ecd.mp2f(mu)
            ))
    }
    use.mpfr <- ifelse(is(sum, "mpfr"), TRUE, FALSE)

    if (use.mpfr) {
        sigma <- ecd.mp1 * sigma
    }

    ld <- new("ecld", call = call,
               lambda = unname(lambda),
               sigma = unname(sigma),
               beta  = unname(beta),
               mu    = unname(mu),
               epsilon = unname(epsilon),
               rho = unname(rho),
               use.mpfr = unname(use.mpfr),
               is.sged = is.sged,
               ecd = new("ecd"),
               ecd_RN = new("ecd"),
               status = 1L
          )
    
    # -------------
    # SGED
    if (is.sged) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: SGED"))
        if (with.ecd | with.mu_D | with.RN) {
            stop("SGED: with.ecd, with.mu_D, with.RN not supported")
        }
        return(invisible(ld))
    }
    
    # -------------
    # ecd
    if (with.ecd | with.mu_D | with.RN) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: calc ecd"))
        ld@ecd <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(beta),
                      sigma=sigma, mu=mu, verbose=verbose)
        ld@status <- bitwOr(ld@status, 2L)
    }
    if (with.mu_D | with.RN) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: calc mu_D"))
        ld@mu_D <- -log(ecd.imgf(ld@ecd, verbose=verbose))
        ld@status <- bitwOr(ld@status, 4L)
    }
    if (with.RN) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: calc ecd_RN"))
        ld@ecd_RN <- ecd(lambda=ld@ecd@lambda, beta=ld@ecd@beta,
                         sigma=ld@ecd@sigma, mu=ld@mu_D, verbose=verbose)
        ld@status <- bitwOr(ld@status, 8L)
    }

    if (verbose) print(paste(Sys.time(), "ecld constructor: done"))
   
    invisible(ld)
}
### <---------------------------------------------------------------------->
#' @rdname ecld
"ecld.from" <- function(object,
                        with.ecd = FALSE, with.mu_D = FALSE,
                        with.RN = FALSE, verbose=FALSE)
{
    if (! is(object, "ecd")) {
        stop("Must come from an ecd object")
    }
    if (object@alpha != 0 | object@gamma != 0) {
        stop("Must come from an ecd object with lambda-only parametrization")
    }
    
    ecld(lambda = object@lambda, sigma = object@sigma, beta = object@beta,
         mu = object@mu, # neither epsilon nor rho is in ecd object
         with.ecd = with.ecd, with.mu_D = with.mu_D,
         with.RN = with.RN, verbose=verbose)
}
### <---------------------------------------------------------------------->
#' @rdname ecld
"ecld.validate" <- function(object, sged.allowed=FALSE, sged.only=FALSE)
{
    if (! is(object, "ecld")) {
        stop(paste("Object must be an ecld object, found:", class(object)))
    }
    if (sged.only) {
        if (!object@is.sged) { # ECLD
            stop("Object must be an SGED object")
        } else {
            return(TRUE)
        }
    }
    if (object@is.sged) {
        if (sged.allowed) return(TRUE)
        stop("SGED object is not allowed here")
    }
    return(TRUE)
}
### <---------------------------------------------------------------------->
#' @rdname ecld
"ecld.quartic" <- function(sigma, epsilon, rho, mu_plus_ratio=NaN, mu_plus=NaN)
{
    # mu_plus and mu_plus_ratio are mutually exclusive
    if ((! is.na(mu_plus)) & (! is.na(mu_plus_ratio))) {
        stop("Only one of mu_plus and mu_plus_ratio can have value")
    }

    sigma <- sigma * ecd.mp1 # force to MPFR

    ld0 <- ecld(lambda=4, sigma=sigma)
    if (! is.na(mu_plus_ratio)) {
        sd0 <- ecld.sd(ld0)
        mu_plus <- mu_plus_ratio * sd0
    }
    if (is.na(mu_plus)) stop("Can not determine mu_plus")
    
    mu_D <- ecld.mu_D_quartic(ld0)
    ld1 <- ecld(lambda=4, sigma=sigma, mu=mu_D+mu_plus)
    ld1@mu_D <- mu_D
    ld1@epsilon <- epsilon
    ld1@rho <- rho
    
    return(ld1)
}
### <======================================================================>
#' @rdname ecld
"ecld.from_sd" <- function(lambda = 3, sd = 1, beta = 0, mu = 0)
{
    one <- sd*0 + 1
    
    sigma = ecd.mp2f(sd * sqrt(gamma(lambda/2)/gamma(lambda*3/2)))
    if (beta != 0) {
        sigma0 = sigma
        match_sd <- function(sigma) {
            ld0 <- ecld(lambda=lambda, sigma=sigma, beta=beta)
            ecd.mp2f(ecld.sd(ld0))-sd
        }
        rs <- uniroot(match_sd, lower=sigma0*0.1, upper=sigma0*10)
        sigma <- rs$root
    }
    sigma <- sigma * one
    object <- ecld(lambda=lambda, sigma=sigma, beta=beta, mu=mu)
    object@mu_D <- tryCatch(
        ecld.mu_D(object),
        error = function(cond) { return(NaN) }
    )

    ecld.validate(object, sged.allowed=TRUE)
    return(object)
}
### <---------------------------------------------------------------------->

