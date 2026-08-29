#### R functions for the -- noncentral chisquare density -- dchisq(.,., ncp > 0)

### ---> ./chisq-nonc-ex.R
###      ~~~~~~~~~~~~~~~~~

## source("~/R/MM/NUMERICS/dpq-functions/dpq-h.R")# the macros

## R's builtin  dchisq(x, df, ncp) is now {since P.Dalgaards improvement}
## -----------                     definitely better than this:
dnoncentchisq <- function(x, df, ncp, kmax = floor(ncp/2 + 5 * (ncp/2)^0.5))
{
    ## The noncentral-chisquare density at x
    ## for df degrees of freedom, noncentrality parameter ncp.
    ## (Only) x may be a vector.
    ##
    if(length(ncp) > 1 || ncp < 0)
        stop("noncentrality parameter ncp must be scalar >= 0!")
    if(length(df) > 1 || df <= 0)
        stop("df must be scalar > 0!")
    kv <- 0:kmax
    poiv <- if(ncp == 0) 1 else dpois(kv, ncp/2)
    n <- length(x)
    fv <- matrix(0, n, 1 + kmax)
    for(k in kv) ## FIXME: compute with one dchisq() call, outer(kv, x, dchisq(.)) ?!!
        fv[, 1L + k] <- dchisq(x, 2 * k + df)
    c(fv %*% poiv)
}

dchisqAsym <- function(x, df, ncp, log = FALSE)
{
  ## Purpose: Asymptotic approximation of NON-central by central Chi^2
  ## ----------------------------------------------------------------------
  ## Arguments: as for dchisq()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  1 Apr 2008, 10:34
    nl <- df + ncp
    n2l <- nl + ncp
    ic <- nl/n2l # = 1/c = 1 / (1+b) -- 'b' of Abramowitz & Stegun, p.942 [26.4.27]
    ff <- dchisq(x*ic, df = nl*ic, log = log)
    if(log) log(ic) + ff else ic * ff
}

dnchisqBessel <- function(x, df, ncp, log = FALSE)
{
    ## Purpose: Implement Fisher(1928) = Johnson,Kotz & Bala. (29.4) [p.436]
    ##          == Exact Formula for non-central dchisq(), using besselI()
    ## ----------------------------------------------------------------------
    ## Arguments: as for dchisq()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr 2008, 21:09

    if (anyNA(x) || anyNA(df) || anyNA(ncp))
	return(x + df + ncp)
    if (any(ncp <= 0) || any(df < 0))
	stop("must have ncp > 0 (here), df >= 0")
    if (!all(is.finite(df)) || !all(is.finite(ncp)))
	stop("some (ncp, df) values are not finite")

    ## Recycle:
    n <- max(length(x), length(df), length(ncp))
    if(n > 1) {
        if(length(x) < n) x <- rep(x, length = n)
        if(length(df) < n) df <- rep(df, length = n)
        if(length(ncp) < n) ncp <- rep(ncp, length = n)
    }

    nu. <- (df - 2)/2
    y <- sqrt(ncp*x)
    ## Iy := besselI(y, ., expon.scaled=TRUE) = exp(-y) * I_.(y)
    Iy <- besselI(y, nu = nu., expon.scaled = TRUE)
    ## NOTE: Bessel(y) for large y --> ../bessel-large-x.R
    ## NOTE 2: besselI(..) also needs a  'log = FALSE' argument <<< !
    if(log) {
        (y -(ncp+x)/2) + (nu./2) * log(x/ncp) + log(Iy/2)
    }
    else { ## not log
        exp(y -(ncp+x)/2) * (x/ncp)^ (nu./2) * Iy/2
    }
}

.checkFun4 <- function(F, chF) {
    if(!is.function(F))
        stop(chF, "is not a function")
    nff <- names(formals(F))
    if(!(("..." %in% nff) || length(nff) >= 4L))
        stop(chF, "() must allow 4 arguments (*, df, ncp, log)")
    ## else NULL
}

pl2curves <- function(fun1, fun2, df, ncp, log = FALSE,
                      from = 0, to = 2*ncp, p.log = "", n = 2001,
                      leg = TRUE,
                      col2 = 2, lwd2 = 2, lty2 = 3, ...)
{
    ## Purpose: Comparison plot of two curves to compare them,
    ##    e.g.  dchisq() and its Bessel-approximation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Aug 2019
    cF1 <- deparse(substitute(fun1)); .checkFun4(fun1, cF1)
    cF2 <- deparse(substitute(fun2)); .checkFun4(fun2, cF2)
    x <- NULL # -Wall(codetools)
    curve(fun1(x, df=df, ncp=ncp, log=log), from=from, to=to, n=n,
          ylab = paste0(paste0(c(cF1,cF2),"(x, *)", collapse = " vs. ")),
          main = deparse(sys.call()), log=p.log, ...)
    curve(fun2(x, df=df, ncp=ncp, log=log), n=n, ##=
          col=col2, lty=lty2, lwd=lwd2, add=TRUE)
    if(leg)
        legend("topright", paste0(c(cF1,cF2),"()"),
               col = c(col2,1L), lwd=c(lwd2,1L), lty=c(lty2,1L), inset = .02)
}

if(FALSE) ## previously -- now just use
    ## pl2curves(dnchisqBessel, dchisq, df, ncp, log,  from, to, ....)
p.dnchiB <- function(df, ncp, log=FALSE, from=0, to = 2*ncp, p.log="", ...)
{
    ## Purpose: Comparison plot of dchisq() and its Bessel-approximation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr 2008, 21:40
    x <- NULL # -Wall(codetools)
    curve(dnchisqBessel(x, df=df, ncp=ncp, log=log), from=from, to=to,
          n=2001, ylab = "dnchisq[Bessel](x,*)",
          main = deparse(sys.call()), log=p.log, ...)
    curve(dchisq       (x, df=df, ncp=ncp, log=log), col=2,   n=4001,
          lty=3, lwd=2, add=TRUE)
    legend("topright", c("dchisq()", "dnchisqBessel()"),
           col = 2:1, lwd=2:1, lty=c(3,1), inset = .02)
}



###-- The R version of R's C implementation ~/R/D/r-devel/R/src/nmath/dnchisq.c

## TODO: add 'verbose' or 'info = FALSE'; if true, return info about i_upper, i_max, i_lower, etc
## TODO(2): use dpois_raw() instead of dpois(), at least optionally, as R's C code does so too
dnchisqR <- function(x, df, ncp, log = FALSE,
                     eps = 5e-15, termSml = 1e-10, ncpLarge = 1000)
{
    ## allow x to be "general" (potentially "mpfr" here)
    stopifnot(length(df) == 1, length(ncp) == 1, length(log) == 1,
              is.numeric(df), is.numeric(ncp), is.logical(log))

    if(is.na(df) || is.na(ncp) || length(x) == 0)
	return(x + df + ncp)
    if (ncp < 0 || df <= 0 || !is.finite(df) || !is.finite(ncp))
	return(x + NaN)
    ncp2 <- 0.5 * ncp

    ## For better debugging etc, we *name* the function here:
    dnchisqR1 <- function(x) {
	if(is.na(x)) return(x)
	if(x < 0) return(.D_0(log))
	if(x == 0 && df < 2.)
	    return(Inf)
	##     if(ncp == 0)
	## 	return(dchisq(x, df, log=log))

	##/* find max element of sum */
	imax <- ceiling((-(2+df) + sqrt((2-df) * (2-df) + 4 * ncp * x)) / 4)
	if (imax < 0)
	    imax <- 0
	if(is.finite(imax)) {
	    dfmid  <- df + 2 * imax
	    ## mid = dpois_raw(imax, ncp2, FALSE) * dchisq(x, dfmid, FALSE)
	    mid <- dpois(imax, ncp2, log=FALSE) * dchisq(x, dfmid, log=FALSE)
	} else mid <- 0

	if(mid == 0) {
	    ##/* underflow to 0 -- maybe numerically correct; maybe can be more accurate,
	    ## particularly when  give_log = TRUE */
	    ##/* Use  central-chisq approximation formula when appropriate;
	    ##    * ((FIXME: the optimal cutoff also depends on (x,df);  use always here? )) */
	    if(log || ncp > ncpLarge) {
		nl <- df + ncp; ic <- nl/(nl + ncp) ##/* = "1/(1+b)" Abramowitz & St.*/
		return(dchisq(x*ic, nl*ic, log=log))
	    } else return(.D_0(log))
	}
	sum <- mid

	##/* errorbound := term * q / (1-q)  now subsumed in while() / if() below: */

	##/* upper tail */
	term <- mid; df <- dfmid; i <- imax
	repeat {
	    i <- i+1
	    q <-  x * ncp2 / i / df;
	    df <- df+2
	    term <- term*q
	    sum <- sum + term
	    if(!(q >= 1 || term * q > (1-q)*eps || term > termSml*sum)) break
	}
	##/* lower tail */
	term <- mid; df <- dfmid; i <- imax
	while(i) {
	    df <- df - 2
	    q <- i * df / x / ncp2
	    i <- i-1
	    term <- term * q
	    sum <- sum + term
	    if (q < 1 && term * q <= (1-q)*eps) break
	}
	## return
	.D_val(sum, log)
    }## end{dnchisqR1}

    ## apply on all x[] :
    vapply(x, FUN.VALUE = x[[1L]], dnchisqR1)
}
