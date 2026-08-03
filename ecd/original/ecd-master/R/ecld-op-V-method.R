#' The O, V, U operators in option pricing model
#'
#' The O operator takes a vector of implied volatility \eqn{\sigma_1(k)}
#' and transforms them to a vector of normalized option prices.
#' The V operator takes a vector of normalized option prices and transforms
#' them to a vector of implied volatility \eqn{\sigma_1(k)}.
#' If \code{ttm} is provided, \eqn{\sigma_1(k)} will be divided by square root of \code{2 ttm} and yield Black-Scholes implied volatility.
#' The U operator calculates the log-slope of the option prices.
#' The op_VL_quartic operator is the quartic composite of V x OGF, assuming epsilon and rho are deposited in the ecld object.
#' The \code{RN} parameter for OGF is not available here. It is always assumed to be \code{FALSE}.
#'
#' @param sigma1 numeric, a vector of implied volatility (without T)
#' @param L numeric, a vector of normalized local option prices
#' @param k numeric, a vector of log-strike
#' @param ttm numeric, time to expiration (maturity), measured by fraction of year.
#'        If specified, V operator will adjust \eqn{\sigma_1(k)} to Black-Scholes implied volatility.
#'        Default is NaN.
#' @param rho numeric, specify the shift in the global mu.
#' @param sd numeric, the stdev of the distribution.
#'           Instead, if an ecld or ecd object is provided,
#'           the stdev will be calculated from it.
#' @param n numeric, number of lags in \code{ecld.op_U_lag}.
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#' @param stop.on.na logical, to stop if fails to find solution.
#'                   Default is to use NaN and not stop.
#' @param use.mc logical, to use mclapply, or else just use for loop.
#'               Default is \code{TRUE}.
#'               For loop option is typically for debugging.
#' @param object an object of ecld class created from \code{ecld.quartic}.
#'               This object contains the full quartic lambda model spec in order
#'               to be used in \code{ecld.op_VL_quartic}
#'
#' @return a numeric vector
#'
#' @keywords ogf
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.op_O
#' @export ecld.op_V
#' @export ecld.op_U_lag
#' @export ecld.op_VL_quartic
#'
### <======================================================================>
"ecld.op_V" <- function(L, k, otype="c", ttm=NaN, rho=0,
                        stop.on.na=FALSE, use.mc=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    use.mpfr <- ifelse(class(L)=="mpfr" | class(k)=="mpfr", TRUE, FALSE)

    len <- length(L)
    if (len != length(k)) {
        stop("Length of L and k must match!")
    }
    if (len > 1 & use.mc==TRUE) {
        f <- function(i) ecld.op_V(L[i], k[i], ttm=ttm, rho=rho,
                                   otype=otype, stop.on.na=stop.on.na)
        s1 <- parallel::mclapply(1:len, f)
        s1 <- if (use.mpfr) ecd.mpfr(s1) else simplify2array(s1)
        return(s1)
    }
    # handle length-one numeric if use.mc
    
    s1 <- L*NaN
    for (i in 1:length(s1)) {
        df <- function(s) ecld.op_O(s, k[i], otype=otype, rho=rho) - L[i]
        if (is.na(L[i])) next
        if (L[i] <= 0) next
        
        retry <- 5
        lower = 0.01*ecd.mp1 # for implied volitility this is already very low
        while (retry > 0 & lower > 0.0000001*L[i] & df(lower) > 0) {
            lower <- lower/10
            retry <- retry-1
        }
        if (df(lower) > 0 & stop.on.na) {
            stop(paste("Failed to find starting lower for i=", i, "L=", L[i], "at k=", k[i]))
        }
        
        upper = L[i]*100**ecd.mp1
        if (upper < lower) upper <- lower*10
        while ( upper < 100 & df(upper) < 0) {
            upper <- upper*10
        }
        if (df(upper) < 0 & stop.on.na) {
            stop(paste("Failed to find starting upper for i=", i, "L=", L[i], "at k=", k[i]))
        }
        if (upper < lower & stop.on.na) {
            stop(paste("Failed to find starting lower/upper:", lower, upper, "for i=", i, "L=", L[i], "at k=", k[i]))
        }
        # print(c(i*ecd.mp1, L[i], k[i], lower, upper)) # debug

        # we use two iterations of uniroot to improve the precision to 10^-5 no matter how small sigma1 is
        if (df(lower) * df(upper) < 0 & upper > lower) {
            rt <- ecd.uniroot(df, lower=lower, upper=upper, use.mpfr=use.mpfr)
            s_appr <- rt$root
            df2 <- function(sn) ecld.op_O(sn*s_appr, k[i], otype=otype, rho=rho) - L[i]
            rt2 <- ecd.uniroot(df2, lower=0.9, upper=1.1, use.mpfr=use.mpfr)
            s1[i] <- s_appr * rt2$root
        } else {
            s1[i] <- NaN
        }
        if (is.na(s1[i]) & stop.on.na) {
            stop(paste("Failed to find root:", lower, upper, "for i=", i, "L=", L[i], "at k=", k[i]))
        }

    }
    
    if (is.na(ttm)) {
        return(s1)
    } else {
        return(s1/sqrt(2*ttm))
    }
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_V
"ecld.op_O" <- function(sigma1, k, otype="c", rho=0)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    use.mpfr <- ifelse(class(sigma1)=="mpfr" | class(k)=="mpfr", TRUE, FALSE)
    
    # use MPFR for erf
    sigma1 <- sigma1 * ecd.mp1
    M1k <- exp(rho)-exp(k)
    p <- -exp(rho)/2 * ecd.erf((k-rho)/sigma1 - sigma1/4)
    q <- exp(k)/2 * ecd.erf((k-rho)/sigma1 + sigma1/4)

    sgn <- if (otype=="c") 1 else -1
    L <- p + q + sgn*M1k/2
    
    if (use.mpfr) return(L) else return(ecd.mp2f(L))
    
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_V
"ecld.op_U_lag" <- function(L, k, sd, n=2)
{
    if (is(sd,"ecld")) sd <- ecld.sd(sd)
    else if (is(sd,"ecd")) sd <- ecd.sd(sd)
    
    n2 <- floor(n/2)
    logL <- log(L)
    dlogL <- ecd.lag(logL, -n2) - ecd.lag(logL, n-n2)
    dk <- ecd.lag(k,-n2) - ecd.lag(k, n-n2)
    slope <- dlogL/dk*sd
    ecd.mp2f(slope)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_V
"ecld.op_VL_quartic" <- function(object, k, otype="c", ttm=NaN,
                                 stop.on.na=FALSE, use.mc=TRUE)
{
    stopifnot(object@lambda==4 & object@beta==0)

    L <- NULL
    epsilon <- if(is.na(object@epsilon)) 0 else object@epsilon
    rho <- if(is.na(object@rho)) 0 else object@rho
    
    L <- ecld.ogf_quartic(object, k-rho, otype=otype, RN=FALSE)
    L <- L + epsilon
    IV <- ecld.op_V(L, k-rho, otype=otype, ttm=ttm, stop.on.na=stop.on.na, use.mc=use.mc)
    return(IV)
}
### <---------------------------------------------------------------------->

