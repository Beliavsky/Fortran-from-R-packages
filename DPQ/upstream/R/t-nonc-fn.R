#### Functions for the Non-central t-distribution
####                   ===============================================
#### Notation etc from Johnson, Kotz and Balakrishnan (1995)
####                   Chapter 31, '5 Distribution Function', p.514 ff
####                   ===============================================
### Ref.:   (( \cite{JohNKB95} [/u/sfs/bib/Assbib.bib] ))
###  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995).
###  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley Series Prob...

#### Used to be part of ./t-nonc-approx.R, see also ./pt-ex.R
####                      ---------------             -------

## for .DT_val():
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")

### TODO:  E.pnt = E [ . ] = delta * sqrt(nu/2)*gamma(.5*(nu-1))/gamma(.5*nu)
###        -----  is very similar to  b_chi (and can use same asymptotic)

## was  b.nu <-
b_chi <- function(nu, one.minus = FALSE, c1 = 341, c2 = 1000)
{
    ## Purpose: Compute  E[ chi_nu ] / sqrt(nu)    --- useful for non-central t
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Mächler, Date:  6 Feb 1999, 11:21;  Jan.2, 2015
    ## REF: Johnson, Kotz,... , p.520, after (31.26a)
    stopifnot(c1 > 0, c2 >= c1)
    b1 <- function(nu) {
        r <- sqrt(2/nu)*gamma(.5*(nu+1))/gamma(.5*nu)
        if(one.minus) 1-r else r
    }
    b2 <- function(nu) {
        logr <- log(2/nu)/2 + lgamma(.5*(nu+1)) - lgamma(.5*nu)
        ## NOTA BENE: should be much better than 'b1()'  iff one.minus and r << 1
        if(one.minus) -expm1(logr) else exp(logr)
    }
    ## Using boolean vars instead of inefficient  ifelse() :
    r <- nu # will be the result
    r[nu==0] <- if(one.minus) 1 else 0
    ## 3 principal regions
    B1 <- 0  < nu & nu <= c1
    B2 <- c1 < nu & nu <= c2
    BL <- c2 < nu  # "L" : Large
    r[B1] <- b1(r[B1])
    r[B2] <- b2(r[B2])
    r[BL] <- b_chiAsymp(nu[BL], one.minus=one.minus)
    ##       ----------
    r
}

## was  b.nu.asymp <-
b_chiAsymp <- function(nu, order = 2, one.minus = FALSE)
{
  ## Purpose: Compute  E[ chi_nu ] / sqrt(nu)  --- for "LARGE" nu
  ## ----------------------------------------------------------------------
  ## Arguments: nu >=0  (degrees of freedom)
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 11 Mar 1999 (order = 2);  Aug 21 2018 (order in 1..5)

  ## Abramowitz & Stegun, 6.1.47 (p.257) for a= 1/2, b=0 : __ for  z --> oo ___
  ## gamma(z + 1/2) / gamma(z) ~ sqrt(z)*(1 - 1/(8z) + 1/(128 z^2) + O(1/z^3))
  ## i.e. for b_chi), z = nu/2;  b_chi) = sqrt(1/z)* gamma((nu+1)/2) / gamma(nu/2)
  ##  b_chi) = 1 - 1/(8z) + 1/(128 z^2) + O(1/z^3) ~ 1 - 1/8z * (1 - 1/(16z))
  ##
  ## For order >= 3, I've used Maple expansion etc ==> ~/maple/gamma-exp2.txt
  stopifnot(length(order) == 1L, order == as.integer(order), order >= 1)
  r <- 1/(4*nu)
  r <- r * switch(order, # polynomial order
                  1,			# 1
                  1 - r/2,              # 2
                  1 - r/2*(1+ r* 5),	# 3
                  1 - r/2*(1+ r*(5 - r/4* 21)),	# 4
                  1 - r/2*(1+ r*(5 - r/4*(21 + 399*r))),# 5
                  stop("Need 'order <= 5', but order=",order))
  if(one.minus) r else 1-r
}

## log(b_chi(nu)) --- Direct log( sqrt(2/nu)*gamma(.5*(nu+1))/gamma(.5*nu) )
lb_chi00 <- function(nu) {
    n2 <- nu/2
    log(gamma(n2 + 0.5)/ gamma(n2) / sqrt(n2))
}
lb_chi0 <- function(nu) { # notably when 'nu' is "mpfr"
    n2 <- nu/2
    lgamma(n2 + 0.5) - lgamma(n2) - log(n2)/2
}

##
lb_chiAsymp <- function(nu, order)
{
    ## Purpose: Asymptotic expansion (nu -> Inf) of  log(b_chi(nu))
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Aug 23 2018; Apr 2021

    stopifnot(length(order) == 1L, order == as.integer(order), order >= 1)
    ## You can derive the first term from
    ## Abramowitz & Stegun, 6.1.47 (p.257) for a= 1/2, b=0 , or probably
    ## using the Stirling formula for log(Gamma(.))
    if(order == 1)
        return(- 1/4/nu )
    r <- 1/(2*nu) # --> 0
    ## I've used Maple expansion etc ==> ~/maple/gamma-exp2.mw (or *.txt versions)
    ## from tex export of Maple [ ~/maple/gamma-exp2.tex ] :
    ## cH := r/2*(-1+(2/3+(-16/5+(272/7+(-7936/9+(353792/11+(-22368256/13+1903757312*r^2*(1/15))*r^2)*r^2)*r^2)*r^2)*r^2)*r^2)
    ##     = r/2*(-1+(2/3+(-16/5+(272/7+(-7936/9+(353792/11+(-22368256/13+1903757312/15*r^2)*r^2)*r^2)*r^2)*r^2)*r^2)*r^2)
    rr <- r*r # := r^2
    ##     = r/2*(-1+(2/3+(-16/5+(272/7+(-7936/9+(353792/11+(-22368256/13+1903757312/15*rr)*rr)*rr)*rr)*rr)*rr)*rr)
    O <- rr*0 # in correct (Rmpfr) precision
    -r/2 *
        switch(order, # polynomial order {written to use full prec w Rmpfr}
               1,			     # 1  (degree 1)
               1 - rr*2/3,                   # 2  (deg.   3)
               1 - rr*((O+2)/3 - rr*16/5),   # 3  (deg.   5)
               1 - rr*((O+2)/3 - rr*((O+16)/5 - rr* 272 /7))   # 4  (deg. 7)
             , 1 - rr*((O+2)/3 - rr*((O+16)/5 - rr*((O+272)/7 - rr* 7936/9))) # 5 (deg. 9)
             , 1 - rr*((O+2)/3 - rr*((O+16)/5 - rr*((O+272)/7 - rr*((O+7936)/9 - rr* 353792/11)))) # 6 (deg. 11)
             , 1 - rr*((O+2)/3 - rr*((O+16)/5 - rr*((O+272)/7 - rr*((O+7936)/9 - rr*((O+353792)/11
                 -rr* 22368256/13))))) # 7 (deg. 13)
             , 1 - rr*((O+2)/3 - rr*((O+16)/5 - rr*((O+272)/7 - rr*((O+7936)/9 - rr*((O+353792)/11
                 -rr*((O+22368256)/13+rr*1903757312/15)))))) # 8 (deg.15)
             , stop("Currently need 'order <= 8', but order=",order))
}


## was  'pt.appr1'
pntLrg <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    ## Approximate pt()  [non-central t] --- using the "extreme case" formula
    ##  from C's pnt.c  and  pntR() below:
    ##
    ## if (df > 4e5 || del*del > 2*M_LN2*(-(DBL_MIN_EXP))) {
    ## /*-- 2nd part: if del > sqrt(2*log(2)*1021) = 37.62.., then p=0 below
    ## FIXME: test should depend on `df', `tt' AND `del' ! */
    ## /* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */

    ## Vectorizing (in t) {{is this correct? -- FIXME??}}
    n <- max(length(t), length(df), length(ncp))
    if(length( t ) != n) t   <- rep_len(t,  n)
    if(length( df) != n) df  <- rep_len(df, n)
    if(length(ncp) != n) ncp <- rep_len(ncp,n)

    neg <- t < 0
    t  [neg] <- -   t[neg]
    ncp[neg] <- - ncp[neg]
    s <- 1/(4*df)
    pnorm(t*(1 - s), mean = ncp, sd = sqrt(1 + t*t*2*s),
          lower.tail = (lower.tail != neg), log.p=log.p)
}

## was called 'pt.appr'
pntJW39.0 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE)
{
  ## Purpose: Jennett & Welch (1939) approximation to non-central t
  ##          see Johnson, Kotz, Bala... Vol.2, 2nd ed.(1995)
  ##                                     p.520, after (31.26a)
  ## -- still works FAST for huge ncp
  ##    but has *wrong* asymptotic tail (for |t| -> oo)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?pt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date:  6 Feb 1999, 11 Mar 1999

  b <- b_chi(df)
  ##   =====
  ## FIXME:  (1 - b^2) below suffers from severe cancellation!
  pnorm((t*b - ncp)/sqrt(1+ t*t*(1 - b*b)),
        lower.tail = lower.tail, log.p = log.p)
}

pntJW39 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE)
{
  ## Purpose: Jennett & Welch (1939) approximation to non-central t
  ##          see Johnson, Kotz, Bala... Vol.2, 2nd ed.(1995)
  ##                                     p.520, after (31.26a)
  ## -- still works FAST for huge ncp
  ##    but has *wrong* asymptotic tail (for |t| -> oo)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?pt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date: Feb/Mar 1999; 1-b^2 improvement: Aug 2018
  ._1_b <- b_chi(df, one.minus=TRUE)# == 1 - b -- needed for good (1 - b^2)
  ##       =====
  b <- 1 - ._1_b
  ## (1 - b^2) == (1 - b)(1 + b) = ._1_b*(2 - ._1_b)
  pnorm((t*b - ncp)/sqrt(1+ t*t * ._1_b*(1 + b)),
        lower.tail = lower.tail, log.p = log.p)
}

## was  c.nu <-
c_dt <- function(nu) {
    ## Purpose: The constant in  dt(t, nu, log=TRUE) =
    ##         = \log f_{\nu}(t) = c_{\nu} - \frac{\nu+1}{2}\log(1 + x^2/\nu)
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Nov 2008, 14:26

    warning("FIXME: current c_dt() is poor -- base it on lb_chi(nu) !")

    ## Limit nu -> Inf: c_{\nu} \to - \log(2\pi)/2 = -0.9189385
    r <- nu
    nu. <- nu[I <- nu < 200]
    r[I] <- log(gamma((nu.+1)/2)/gamma(nu./2)/sqrt(pi*nu.))
    nu. <- nu[I <- !I & nu < 1e4]
    r[I] <- lgamma((nu.+1)/2) - lgamma(nu./2) - log(pi*nu.)/2
    I <- nu >= 1e4
    r[I] <- c_dtAsymp(nu[I])
    r
}

## was  c.nu.asymp <-
c_dtAsymp <- function(nu)
{
    ## Purpose: Asymptotic of c_dt -- via Abramowitz & Stegun, 6.1.47 (p.257)
    ## ----------------------------------------------------------------------
    ## Arguments: nu >=0  (degrees of freedom)
    ## ----------------------------------------------------------------------
    ##
    ## FIXME: This is trivially   -log(2*pi)/2  +  log(b_chi(nu))
    ##        and I've computed good asymptotics for
    ##    lb_chi(nu) := log(b_chi(nu))   above
    ##
    ##
    ## -log(2*pi)/2 -1/(4*nu) * (1 + 5/(96*nu^2))
    ##                               ^^^^^^^^^^^ not quite ok
    warning("this is poor -- use lb_chi(nu) !!")
    -log(2*pi)/2 - 1/(4*nu)
}

## was  c.pt
c_pt <- function(nu)
{
  ## Purpose: the asymptotic constant in log F_nu(-x) ~= const(nu) - nu * log(x)
  ##          where F_nu(x) == pt(x, nu)
### FIXME == Source / Reference for the above left tail statement?
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Nov 2008, 09:34
  warning("use better c_dt()") ## ==> use lb_chi(nu) !!
  c_dt(nu) + (nu-1)/2 * log(nu)
}


## MM:  My conclusion of experiments in ../tests/pnt-prec.R :
## -------              --------------------------
## Large t (>0 or < 0)  *MUST* get a new algorithm !
## -------              --------------------------
## >>>> ../man/pnt.Rd <<<<<<<<<
pntR1  <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                   use.pnorm =
                       (df > 4e5 ||
                        ncp^2 > 2*log(2) * 1021), # .Machine$double.min.exp = 1022
                   ## /*-- 2nd part: if del > 37.6403, then p=0 below
                   ## FIXME: test should depend on `df`, `t` AND `ncp`
                   itrmax = 1000, errmax = 1e-12, verbose = TRUE)
{
    ## Purpose: R version of the series used in pnt() in
    ##          ~/R/D/r-devel/R/src/nmath/pnt.c
    ## ----------------------------------------------------------------------
    ## Arguments: same as  pt()
    ## Author: Martin Mächler, Date:  3 Apr 1999, 23:43

    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1,
              df > 0, ncp >= 0) ## ncp == 0 --> pt()
    ## Make this workable also for "mpfr" objects --> use format() in cat()
    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(t, df, ncp)
    if(!isN) {
	if(verbose) cat("some 'non-numeric arguments .. fine\n")
	if(isMpfr) {
	    stopifnot(requireNamespace("Rmpfr"))
            getPrec <- Rmpfr::getPrec
	    prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
	    pi <- Rmpfr::Const("pi", prec = max(64, prec))
	} ##--
    }

    if (t >= 0) {
        negdel <- FALSE; tt <- t; del <- ncp
    } else {
        if(verbose) cat("t < 0  ==> swap delta := -ncp\n")
	negdel <- TRUE ; tt <- -t; del <- -ncp
    }

    if (use.pnorm) {
        ## Approx. from Abramowitz & Stegun 26.7.10 (p.949) -- FIXME (see above)
	s <- 1./(4.*df)
        pnt.. <- pnorm(tt*(1 - s), del, sqrt(1. + tt*tt*2.*s),
                       lower.tail = (lower.tail != negdel), log.p=log.p)
        cat("large 'df' or \"large\" 'ncp' ---> return()ing pnorm(*) =",
            format(pnt.., digits=16), "\n")
        return(pnt..)
    }

    ## initialize twin series */
    ## Guenther, J. (1978). Statist. Computn. Simuln. vol.6, 199. */

    x <- t * t
    rxb <- df/ (x + df) ## := (1 - x) {x : below} -- but more accurately
    x   <- x / (x + df) # in [0,1) -- x == 1 for very large t */
    ## (1 - x) = df / (x + df)

    if(verbose)
        cat("pnt(t=",format(t),", df=",format(df),", delta=",format(ncp),") ==> x= ",
            format(x),":")

    tnc <- 0 # tnc will be the result */
    if (x > 0) { ## <==>  t != 0 */
	lambda <- del * del
	p <- 0.5 * exp(-0.5 * lambda)
	if(verbose) cat("\t p=",format(p),"\n")

	if(p == 0.) { # underflow! */

## FIXME: should "analytically  factor-out the small 0.5 * exp(-0.5 * lambda)
## -----  and re-multiply at end --- where we need *logspace* addition of
## in      tnc <-  tnc*0.5 * exp(-0.5 * lambda)  +  pnorm(- del, 0, 1)
## in      tnc <-  exp(log(tnc*0.5 * exp(-0.5 * lambda))  + exp(pnorm(- del, 0, 1, log=TRUE))
## in      log.tnc <-  log(exp(AA) + exp(BB)) === copula:::lsum(rbind(AA, BB))
##            where  AA <- log(0.5) + log(tnc) - lambda/2
##            and    BB <- pnorm(- del, log.p=TRUE)


            ##========== really use an other algorithm for this case !!! */
            cat("\n__underflow__ -- p = 0 -- |ncp| = |delta|  too large\n\n")
            ## |delta| too large */
            return(0)
	}
        if(verbose >= 2)
            cat(  " it   1e5*(godd,  geven)         p          q         s ",
                ## 1.3 1..4..7.9 1..4..7.9|1..4..7.9  1..4..7.9  1..4..7.9_ */
                  "           pnt(*)    D(pnt)     errbd\n", sep=''
                ## 1..4..7..0..3..67 1..4..7.9 1..4..7.9*/
                )

	q <- sqrt(2/pi) * p * del
	a <- 0.5
	b <- 0.5 * df
	s <- 0.5 - p # but that may suffer from cancallation:
        ## s =  0.5 - p = 0.5*(1 - exp(-L)) =  -0.5*expm1(-L))
        if(s < 1e-7)
            s <- -0.5 * expm1(-0.5 * lambda)
	## was: rxb <- (1 - x) ^ b
        ## (1 - x) = df / (x + df)
        rxb <- rxb ^ b
	albeta <- .5*log(pi) + lgamma(b) - lgamma(0.5 + b)
	xodd <- if(isN) pbeta(x, a, b) else pbetaRv1(x, a, b)
        ##                                  -------- ./beta-gamma-etc/pbetaR.R
	godd <- 2 * rxb * exp(a * log(x) - albeta)
	xeven <- if(b*x <= .Machine$double.eps) b*x else 1 - rxb
        ## xeven = 1 - (1 - x)^b = b*x - choose(b,2) * x^2  + O((bx)^3)
	geven <- b * x * rxb
	tnc <- p * xodd + q * xeven
        errbd <- Inf

        ## repeat until convergence or iteration limit */
	for(it in 1:itrmax) {
	    a <- a+1
	    xodd  <- xodd - godd
	    xeven <- xeven - geven
	    godd  <- godd * x * (a + b - 1) / a
	    geven <- geven* x * (a + b - 0.5) / (a + 0.5)
	    p <- p * lambda / (2 * it)
	    q <- q * lambda / (2 * it + 1)
	    tnc <- tnc + p * xodd + q * xeven
	    s <- s - p
	    if(s < -1e-10){## happens e.g. for (t,df,delta)=(40,10,38.5), after 799 it.*/
                cat("IN loop (it=",it,"): precision NOT reached\n")
		## ML_ERROR(ME_PRECISION)
                if(verbose)
                    cat("s =",format(s)," < 0 !!! ---> non-convergence!!\n")
		break # goto finis
	    }
	    if(s <= 0 && it > 1) break # goto finis
	    errbd <- 2 * s * (xodd - godd)
            if(verbose >= 2)
                cat(sprintf(paste("%3d %#9.4g %#9.4g",
                                  "%#10.4g %#10.4g %#9.4g %#17.15g %#9.4g %#9.4g\n",
                                  sep="|"),
                            it, 1e5*godd, 1e5*geven,
                            p,q, s, tnc, p * xodd + q * xeven ,errbd))
	    if(abs(errbd) < errmax) break # convergence
	}
        ## non-convergence:*/
        if(abs(errbd) >= errmax && s > 0)
            cat("end loop: precision NOT reached\n") ## ML_ERROR(ME_PRECISION)
    } ## x > 0
    tnc0 <- tnc
    tnc <- tnc + pnorm(- del, 0, 1)

    ## was -- if (negdel) 1 - tnc else tnc
    lower.tail <- lower.tail != negdel ## xor

    if(verbose)
        cat(sprintf("%stnc{sum} = %.12g, tnc+pnorm(-del) = %.12g, lower.t = %d\n",
                    if(verbose == 1 && x > 0) sprintf("%d iter.: ", it) else "\\--> ",
                    tnc0, tnc, lower.tail))

    if(tnc > 1 - 1e-10 && lower.tail)
        cat("finis: precision NOT reached\n")

    .DT_val(min(tnc, 1.), lower.tail, log.p)
}# pntR1()

pntR  <- Vectorize(pntR1, c("t", "df", "ncp"))
##==

##' Simple vector version of  copula:::lsum() -- ~/R/Pkgs/copula/R/special-func.R
##' Properly compute log(x_1 + .. + x_n) for given log(x_1),..,log(x_n)
##' Here, x_i > 0  for all i
##'
##' @title Properly compute the logarithm of a sum
##' @param lx n-vector of values log(x_1),..,log(x_n)
##' @param l.off the offset to substract and re-add; ideally in the order of
##'        the maximum of each column
##' @return log(x_1 + .. + x_n) = log(sum(x)) = log(sum(exp(log(x))))
##'         = log(exp(log(x_max))*sum(exp(log(x)-log(x_max))))
##'         = log(x_max) + log(sum(exp(log(x)-log(x_max)))))
##'         = lx.max + log(sum(exp(lx-lx.max)))
##' @author Martin Maechler (originally joint with Marius Hofert)
##'
##' NB: If lx == -Inf for all should give -Inf, but gives NaN / if *one* is +Inf, give +Inf
lsum <- function(lx, l.off = max(lx)) {
    if (is.finite(l.off))
        l.off + log(sum(exp(lx - l.off)))
    else if(missing(l.off) || is.na(l.off) || l.off == max(lx))
        l.off
    else stop("'l.off  is infinite but not == max(.)")
}


##' Simple vector version of  copula:::llsum() -- ~/R/Pkgs/copula/R/special-func.R
##' Properly compute log(x_1 + .. + x_n) for given
##' log(|x_1|),.., log(|x_n|) and corresponding signs sign(x_1),.., sign(x_n)
##' Here, x_i is of arbitrary sign
##' @title compute logarithm of a sum with signed large coefficients
##' @param lxabs n-vector of values log(|x_1|),..,log(|x_n|)
##' @param signs corresponding signs sign(x_1), .., sign(x_n)
##' @param l.off the offset to substract and re-add; ideally in the order of max(.)
##' @param strict logical indicating if it should stop on some negative sums
##' @return log(x_1 + .. + x_n)
##'         log(sum(x)) = log(sum(sign(x)*|x|)) = log(sum(sign(x)*exp(log(|x|))))
##'         = log(exp(log(x0))*sum(signs*exp(log(|x|)-log(x0))))
##'         = log(x0) + log(sum(signs* exp(log(|x|)-log(x0))))
##'         = l.off   + log(sum(signs* exp(lxabs -  l.off  )))
##' @author Martin Maechler (originally joint with Marius Hofert)
lssum <- function (lxabs, signs, l.off = max(lxabs), strict = TRUE)
{
    sum. <- sum(signs * exp(lxabs - l.off))
    if (any(is.nan(sum.) || sum. <= 0))
        (if (strict) stop else warning)("lssum found non-positive sums")
    l.off + log(sum.)
}


### Simple inefficient but hopefully correct version of pntP94..()
### This is really a direct implementation of formula
### (31.50), p.532 of  Johnson, Kotz and Balakrishnan (1995)

pnt3150.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE, M = 1000,
                      verbose = TRUE) {
    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              is.numeric(M), M >= 1, M == round(M))
    if (t < 0)
        return(pnt3150.1(-t, df, ncp = -ncp, lower.tail = !lower.tail,
                         log.p=log.p, M=M))

    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(t, df, ncp)
    if(!isN) {
        if(isMpfr) {
            stopifnot(requireNamespace("Rmpfr"))
            ## if(!exists("pbetaRv1", mode="function"))
            ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
            pbeta <- Vectorize(pbetaRv1, "shape2") # so pbeta(x, p, <vector q>) works

            ## TODO ???? as with pntP94.1() below ???? <<<<<<<<<<<<<<<< ???
            ## ---------
            ## getPrec <- Rmpfr::getPrec
            ## prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
            ## pi <- Rmpfr::Const("pi", prec = max(64, prec))
            ## dbeta <- function(x, a,b, log=FALSE) {
            ##     lval <- (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b)
            ##     if(log) lval else exp(lval)
            ## }
        }
    }

    ## hence now, t >= 0 :
    x <- df / (df + t^2)
    lambda <- 1/2 * ncp^2
    j <- 0:M
    ## Cheap and stupid -- and not ok for moderately large ncp !
    ## terms <- (ncp/sqrt(2))^j / gamma(j/2 + 1) * pbeta(x, df/2, (j+1)/2)
    ## log(.):
    l.terms <- j*(log(ncp)- log(2)/2) - lgamma(j/2 + 1) + pbeta(x, df/2, (j+1)/2, log.p=TRUE)
    if(any(ina <- is.na(l.terms))) {
        ina <- which(ina)
        if(ina[length(ina)] == M+1 && all(diff(ina) == 1)) {
            ## NaN [underflow] only for j >= j_0
            cat("good: NaN's only for j >= ", ina[1],"\n")
            l.terms <- l.terms[seq_len(ina[1]-1)]
        }
    }
    if(verbose) {
        cat(sprintf(" log(terms)[1:%d] :\n", length(l.terms))); print(summary(l.terms))
        elt <- exp(l.terms)
        cat(sprintf("sum(exp(l.terms)) , sum(\"sort\"(exp(l.terms))) = (%.16g, %.16g)\n",
                    sum(elt), sum(exp(sort(l.terms)))))
        cat(sprintf("log(sum(exp(l.terms))) , lsum(l.terms), rel.Delta = (%.16g, %.16g, %.5e)\n",
                    log(sum(elt)), lsum(l.terms), 1 - log(sum(elt))/lsum(l.terms)))
        cat("exp(-delta^2/2) =", format(exp( - ncp^2/2)),"\n")
    }
    ## P(..) = 1 - exp(-lambda)/2 * Sum == 1 - exp(-LS)
    ## where  LS := lambda + log(2) - log(Sum)
    ##  exp(-lambda)/2 * Sum = exp(-lambda - ln(2) + log(Sum))

    ## For non-small ncp, and small t (=> very small P(.)),
    ##  e.g., pt(3, 5, 10) have huge cancellation here :
    LS <- lambda + log(2) - lsum(l.terms)

    if(log.p) {
        if(lower.tail) ## log(1 - exp(-LS)) = log1mexp(LS)
            log1mexp(LS)
        else ## upper tail: log(1 - (1 - exp(-LS))) = -LS
            -LS
    } else {
        if(lower.tail) -expm1(-LS) else exp(-LS)
    }
}
pnt3150  <- Vectorize(pnt3150.1, c("t", "df", "ncp"))

### New version of pntR1(), pntR() ----  Using the  Posten (1994) algorithm
pntP94.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                     itrmax = 1000, errmax = 1e-12, verbose = TRUE)
{
    stopifnot(length(t) == 1, length(df) == 1, length(ncp) == 1,
              df > 0, ncp >= 0) ## ncp == 0 --> pt()

    Cat <- function(...) if(verbose > 0) cat(...)
    ## Make this workable also for "mpfr" objects --> use format() in cat()
    isN <- is.numeric(t) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(t, df, ncp)
    if(!isN) {
        if(isMpfr) {
            stopifnot(requireNamespace("Rmpfr"))
            ## if(!exists("pbetaRv1", mode="function"))
            ##     source("~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/pbetaR.R")
            getPrec <- Rmpfr::getPrec
            prec <- max(getPrec(t), getPrec(df), getPrec(ncp))
            pi <- Rmpfr::Const("pi", prec = max(64, prec))
            pbeta <- pbetaRv1
            dbeta <- function(x, a,b, log=FALSE) {
                lval <- (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b)
                if(log) lval else exp(lval)
            }
        }
        ## if(FALSE) ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
	## .N <- if(isMpfr) Rmpfr::asNumeric else as.numeric
        Cat("Some 'non-numeric arguments .. fine\n")
    }

    neg.t <- (t < 0)
    if(neg.t) {
        if(verbose) cat("t < 0  ==> swap delta := -ncp\n")
	## tt <- -t
        del <- -ncp
    } else {
        ## tt <- t
        del <- ncp
    }

    x <- t * t
    ## x := df / (df + t^2)   <==>  (1-x) =: rxp = t^2 / (df + t^2)
    x   <- df/ (x + df) ## in (0, 1] :  x == 1 for very large t
    rxb <- x / (x + df) # == (1 - x)  in [0,1)
    Cat("pnt(t=",format(t),", df=",format(df),", delta=",format(ncp),") ==> x= ",
        format(x),":")

    lambda <- del * del / 2
    B  <- pbeta(x, df/2, 1/2)
    BB <- pbeta(x, df/2, 1  )
    x.1x <- x*rxb # ==  x (1 - x)
    S  <- 2*x.1x*dbeta(x, df/2, 1/2)
    SS <-   x.1x*dbeta(x, df/2, 1  )

    L <- lambda ## at the end multiply 'Sum' with exp(-L)
    ## "FIXME" e.g. for large lambda -- correct already: rescale (B,BB) and/or (S,SS)
    D <- 1
    E <- D*del*sqrt(2/pi)

    Sum <- D*B + E * BB # will be the result

    ## repeat until convergence or iteration limit */
    for(i in 1:itrmax) {
        B  <-  B + S
        BB <- BB + SS
        D <- lambda/ i        * D
        E <- lambda/(i + 1/2) * E
        term <- D*B + E * BB
        Sum <- Sum + term
        if(verbose >= 2)
            cat("D=",format(D),", E=",format(E),
                "\nD*B=",format(D*B),", E*BB=",format(E*BB),
                "\n -> term=",format(term),", sum=",format(Sum),
                ";\n rel.chng=", format(term/Sum),"\n")
        if(abs(term) <= errmax * abs(Sum))
            break ## convergence

        if(i < itrmax) {
            i2 <- i*2
            S  <- rxb * (df+i2-1)/(i2+1) * S
            SS <- rxb * (df+i2  )/(i2+2) * SS
        }
        else
            warning("Sum did not converge with ", itrmax, "terms")
    }

    if(neg.t) lower.tail <- !lower.tail

    ## FIXME 1: for large lambda [Sum maybe too large (have overflown) as well!]
    ## FIXME 2: if (log.p) do better anyway!

    ## iP = 1 - P(..) = exp(-L)*Sum/2 == exp(-LS), as
    ## exp(-L) * S/2 = exp(-L)*exp(log(S/2)) = exp(-L + log(S) - log(2)) =: exp(-LS)
    ## where  LS := L +log(2) - log(S)
    LS <- L + log(2) - log(Sum)
    if(verbose)
        cat( if(verbose == 1) sprintf("%d iter.: ", i) else "\\--> ",
            "Sum = ", format(Sum),
            "LS = -log(1-P(..))=", format(LS), "\n")

    if(log.p) {
        if(lower.tail) ## log(1 - exp(-L)* S/2) = log(1 - exp(-LS)) = log1mexp(LS)
            log1mexp(LS)
        else
            -LS
    } else {
        if(lower.tail) -expm1(-LS) else exp(-LS)
    }
}
pntP94 <- Vectorize(pntP94.1, c("t", "df", "ncp"))

### Should implement
pntChShP94.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                   itrmax = 1000, errmax = 1e-12, verbose = TRUE) {
    stop("not yet ...")
    ## ~/save/papers/Numerics/Chattamvelli+Sh_noncentr-t_1994.pdf

    ## Chattamvelli, R. and Shanmugam, R. (1994)
    ## An enhanced algorithm for noncentral t-distribution,
    ## \emph{Journal of Statistical Computation and Simulation} \bold{49}, 77--83.
}

pntChShP94 <- Vectorize(pntChShP94.1, c("t", "df", "ncp"))

### ==== Gil et al. (2023) "New asymptotic representations of the noncentral t-distribution"
###      ----------------  DOI 10.1111/sapm.12609 ==  https://doi.org/10.1111/sapm.12609
##' Gil A., Segura J., and Temme N.M. (2023) -- DOI 10.1111/sapm.12609
##' Theorem 6 -- asymptotic for zeta --> Inf;  zeta := \delta^2 y/2;  y = t^2/(df + t^2)
pntGST23_T6.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                          y1.tol = 1e-8, # <- somewhat arbitrary for now
                          Mterms = 20, # in the paper, they used Mterms = 11 for Table 1 (p.869)
                          alt = FALSE, # << TODO: find out "automatically" which (T/F) is better
                          verbose = TRUE)
{
    stopifnot(length(Mterms)==1L, Mterms >= 1, Mterms == (M <- as.integer(Mterms)),
              length(y1.tol)==1L, y1.tol >= 0)

    if(!lower.tail) { ## upper tail .. is *not* optimized for all cases
        if(!log.p)
            1 - pntGST23_T6.1(t, df, ncp, lower.tail=TRUE,
                              y1.tol=y1.tol, Mterms=Mterms, verbose=verbose)
        else { # log.p: log(1 - F(.))
           if(alt) ## use log scale
               ## alternatively, not equivalently -- see FIXME in log.p case below
               log1mexp(pntGST23_T6.1(t, df, ncp, lower.tail=TRUE, log.p=TRUE,
                                      y1.tol=y1.tol, Mterms=Mterms, verbose=verbose))
           else log1p(- pntGST23_T6.1(t, df, ncp, lower.tail=TRUE,
                                      y1.tol=y1.tol, Mterms=Mterms, verbose=verbose))
        }
    }
    t2 <- t^2
    t2df <- t2 + df
    y <- t2/t2df # in [0,1]; (mathematically would be in [0,1), but numerically gets == 1)
    l_y <- df/t2df # := 1-y  {but stable for df << t2  (when t2+df underflows to t2)}
    del2.2 <- ncp^2/2
    zeta <- del2.2 * y
    if(l_y <= y1.tol) {
        if(l_y == 0)
            stop("l_y = 1-y (stable form) = 0; the theorem 6 approximation cannot work")
        if(y == 1)
            warning("y = 1 (numerically); the approximation maybe very bad")
        else
            warning(sprintf("y = t^2/(t^2 + df) is too(?) close to 1: 1-y=%g", l_y))
    }
    if(verbose) cat(sprintf("pntGST23_T6(): ...."))
    ## a0 <- 1
    n2 <- df/2 # used "everywhere" below
    kk <- seq_len(M-1L) # 1,...
    ak <- choose(n2-1/2, kk) # a_k = choose(n2-1/2, k) = (n2+1/2-k)/k * a_{k-1} for k >= 1  (45), p.868
    ## The  c_k  and  f_k := (1 - n/2)_k {*rising* factorial [= *not* as I use Pochhammer in Rmpfr::pochMpfr]
    ## are defined recursively, during the summation of the infinite \sum_k=0^{\infty} sum which
    ## we approximate by the M terms sum k=0:M
    sign <- 1
    zet.k <- 1 # == zeta ^ k
    ## f0 <- 1 - n2
    nk <- -n2 # such that nk+1 = 1-n2 == f1
    fk <- 1 # == (1 - n/2)_0
    ck <- c0 <- 1/l_y
    Sum <- sign* fk * ck / zet.k # = summand for k=0
    if(verbose) cat(sprintf("k= 0: zeta=%13.9g, c0=%13.10g,   Sum= %-24.16g\n", zeta, c0, Sum))
    for(k in kk) { # 1, 2, .., M-1
        sign <- -sign
        fk <- fk*(nk <- nk + 1) # (a)_n = *rising* factorial  (7)
        ck <- (ak[k] + y * ck) / l_y  # (44)
        zet.k <- zet.k * zeta
        Sum <- Sum + sign* fk * ck / zet.k
        if(verbose) cat(sprintf("k=%2d: new summand= %18.12g --> new Sum= %-24.16g\n",
                                k, sign* fk * ck / zet.k, Sum))
    }
    if(verbose) cat("\n")

    if(log.p) {
        zeta - del2.2   +  log(y)/2  + n2*log(l_y) + (n2-1)*log(zeta) - lgamma(n2) + log(Sum)
    } else { ## lower.tail = TRUE,  log.p == FALSE
        if(alt) ## use log scale -- alternatively, not equivalently
            exp(zeta - del2.2  + log(y)/2  + n2*log(l_y) + (n2-1)*log(zeta) - lgamma(n2) + log(Sum))
        else
            exp(zeta - del2.2) * sqrt(y)   * l_y^n2      *  zeta^(n2-1)     / gamma(n2)  * Sum
    }
}
pntGST23_T6 <- Vectorize(pntGST23_T6.1, c("t", "df", "ncp"))


## "Simple accurate" pbeta(), i.e., I_x(p,q) from Gil et al. (2023) -- notably Nico Temme
##  *much* less sophisticated than R's TOMS 708 based pbeta(), *BUT* with potential to be used with Rmpfr
## from Nico Temme's "Table 1" Maple code ------------------------
## using R's beta() instead of their betaint() --> ../Misc/pnt-Gil_etal-2023/Giletal23_Ixpq.R
## to see some evidence, beta() is better.
##
## Vectorized version of I_x(p,q) == pbeta(x, p,q)
##                  vvv 1-x (with full accuracy)
Ixpq <- function(x, l_x = 1/8 - x + 7/8, p, q, tol = 3e-16, it.max = 100L, plotIt=FALSE) {
    stopifnot(is.numeric(tol), length(tol) == 1L, 0 < tol, tol < 1, it.max >= 4L)
    oneMinus <- function(x) 1/8 - x + 7/8 # slightly less cancellation than 1-x (?)
    if(missing(l_x))
        l_x <- oneMinus(x)
    else if(missing(x))
        x <- oneMinus(l_x)
    else
        stopifnot(length(x) == length(l_x),
                  "x and l_x do not add to 1" = abs(x + l_x - 1) <= 1e-15)
    r <- x
    r[N  <- (x <= 0)] <- 0
    r[G1 <- (x >= 1) & l_x <= 0] <- 1
    in01 <- !N & !G1 # inside (0,1)
    if(any(swap <- x[in01] > p/(p+q)))
        r[in01][swap] <- 1 - Ixpq(l_x[in01][swap], x[in01][swap], p = q, q = p,
                                  tol=tol, it.max=it.max, plotIt=plotIt)
    ## else: for indices which(in01)[!swap]
    if(length(ii <- which(in01)[!swap])) {
        x   <- x  [ii]
        l_x <- l_x[ii]
        pq <- p+q
        f <- 1
        fn1 <- 1; gn1 <- 1
        fn2 <- 0; gn2 <- 1
        p. <- p
        m <- 0L
        ## FIXME:  for *some* x, convergence is reached when others are still "far" away
        for(it in 1:it.max) {
            p. <- p. + 1 # p. == p+it
            if(it %% 2 == 0L) { # type(it,even))
                m <- m+1L
                dn <- m*x*(q-m)/((p.-1)*p.)
            } else {
                dn <- -x*(p+m)*(pq+m)/((p.-1)*p.)
            }
            g <- f # previous
            fn <- fn1 + dn*fn2
            gn <- gn1 + dn*gn2
            f <- fn/gn
            if (it > 3L) {
                err <- abs((f-g)/f)
                if(any(f0 <- f == 0)) # absolute error for  f = 0
                    err[f0] <- abs(g[f0])
                if(plotIt) {
                    if(it == 4) {
                        plot(x, err, type="l", log="y", ylim = range(tol, err[err > 0]))
                        legend("topleft", sprintf("Ixpq(x, p=%g, q=%g, tol=%g)", p,q, tol), bty="n")
                        abline(h = (1:2)*2^-52, lty=3, col="thistle")
                    } else
                        lines(x, err, col=adjustcolor(2, 1/4), lwd=2)
                }
                if(all(err <= tol)) ## converged
                    break
            }
            fn2 <- fn1; fn1 <- fn
            gn2 <- gn1; gn1 <- gn
        }
        if(it >= it.max && any(err > tol))
            warning(gettextf("continued fraction computation did not converge (tol=%g, it.max=%d)",
                             tol, it.max), domain=NA)
                                        # lprint(it)
        if(plotIt)
            legend("topleft", inset=1/20,
                   sprintf("--> #iter = %d; range(f) = [%g, %g]", it, min(f), max(f)), bty="n")
        r[ii] <- f * x^p * l_x^q / (p * beta(p,q))
    }
    r
} ## Ixpq()  {vectorized  (too simply, for "production" use)}


### ==== from Gil et al. (2023) ... eq. (1) .. but
##' Using the pbeta sums -- similar to Lenth(1989) but *not* suffering from underflow for large ncp:
##' ---- really based on the _Maple_ code Nico Temme sent to MM (May 2024) ----
pntGST23_1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                       j0max = 1e4, # for now
                       IxpqFUN  = Ixpq,
                       ## y1.tol = 1e-8, # <- somewhat arbitrary for now
                       alt = FALSE, # << TODO: find out "automatically" which (T/F) is better
                       verbose = TRUE,
                       ...  ## <-- further arguments passed to IxpqFUN
                       )
{
    stopifnot("'ncp' must be of length 1 ((for now))" = length(ncp) == 1L,
              length(j0max) == 1L,
              "IxpqFUN() must be an incomplete beta function (x, l_x, p, q, ....)" =
                  is.function(IxpqFUN) && length(formals(IxpqFUN)) >= 4 + ...length() &&
                  all.equal(0.6177193984, IxpqFUN(.4, .6, 4,7), tolerance = 1e-14))

    if(!lower.tail) { ## upper tail .. is *not* optimized for all cases
        if(!log.p)
            1 - pntGST23_1(t, df, ncp, lower.tail=TRUE,
                           ## y1.tol=y1.tol, Mterms=Mterms,
                           verbose=verbose)
        else { # log.p: log(1 - F(.))
           if(alt) ## use log scale
               ## alternatively, not equivalently -- see FIXME in log.p case below
               log1mexp(pntGST23_1(t, df, ncp, lower.tail=TRUE, log.p=TRUE,
                                   ## y1.tol=y1.tol, Mterms=Mterms,
                                   verbose=verbose))
           else log1p(- pntGST23_1(t, df, ncp, lower.tail=TRUE,
                                  ## y1.tol=y1.tol, Mterms=Mterms,
                                  verbose=verbose))
        }
    }


    ## was Fnxback() in  ../Misc/pnt-Gil_etal-2023/Giletal23_Ixpq.R
    ##     ~~~~~~~~~     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## MM: derive the "j0" (number of terms, as function of ncp=delta) from  GST 2023 (101):
    ## --  e^-z z^j  /  (e^-j j^j)        == eps            (101)  |  take log(.)
    ## <==>  -z + j log(z) +j - j log(j)  == log(eps)
    ## <==>  j0fn(j, z, eps) == 0  with j0fn(..) :=
    j0fn <- function(j, z, tol = 1e-16) j-z+j*(log(z) - log(j)) - log(tol)
    j0_solve1 <- function(z, tol, tolRoot)
        uniroot(function(j) j0fn(j, z=z, tol=tol),
                interval = c(z, max(30,2*z)), tol = tolRoot, extendInt = "down")

    ## p <- Pnxbackwrec(df, t, ncp)
    ## q <- Qnxbackwrec(df, t, ncp)
    ## p+q + pnorm(-ncp)

    ## pP <- Pnxbackwrec(df, t, ncp) -----------------------------------------
    t2 <- t^2
    t2df <- t2 + df
    y   <- t2/t2df # in [0,1]; (mathematically would be in [0,1), but numerically gets == 1)
    l_y <- df/t2df # := 1-y  {but stable for df << t2  (when t2+df underflows to t2)}
    z <- ncp^2/2 # 'del2.2'
    j0 <- ceiling(j0_solve1(z, tol = 1e-16, tolRoot = 1e-9)$root)
    stopifnot(2 <= j0, j0 <= j0max)

    p <- 1/2;  q <- df/2
    Ijj <- IxpqFUN(y, l_y, p+j0, q, ...) # Ij[j]   j=j0
    pj <- exp(-z+j0*log(z)-lgamma(j0+1))
    s <- pj*Ijj
    Ij1 <- IxpqFUN(y, l_y, p+j0-1, q, ...) # Ij[j-1]  j=j0
    pj <- pj*j0/z
    s <- s + pj*Ij1
    for(j in (j0-1L):1) { ## while(j > 0)
        ppj <- p+j; pjq.y <- (ppj+q-1)*y
        Ij_1 <- ((ppj+pjq.y)*Ij1 - ppj*Ijj) / pjq.y
        pj <- j*pj/z
        s <- s+ pj*Ij_1
        Ijj <- Ij1  # Ijj = Ij[<next j>  ] = Ij[j-1]
        Ij1 <- Ij_1 # Ij1 = Ij[<next j>-1]
    }
    pP <- s/2

    ## qQ <- Qnxbackwrec(df, t, ncp) -----------------------------------------
    p <- 1
    Ijj <- IxpqFUN(y, l_y, p+j0, q, ...) # Ij[j]   j=j0
    qj <- exp(-z + j0*log(z) - lgamma(j0+3/2))
    s <- qj*Ijj
    Ij1 <- IxpqFUN(y, l_y, p+j0-1, q, ...) # Ij[j-1]  j=j0
    qj <- qj*(j0+1/2)/z
    s <- s + qj*Ij1
    for(j in (j0-1L):1) { ## while(j > 0)
        ppj <- p+j; pjq.y <- (ppj+q-1)*y
        Ij_1 <- ((ppj+pjq.y)*Ij1 - ppj*Ijj) / pjq.y
        qj <- (j+1/2)*qj/z
        s <- s+ qj*Ij_1
        Ijj <- Ij1  # Ijj = Ij[<next j>  ] = Ij[j-1]
        Ij1 <- Ij_1 # Ij1 = Ij[<next j>-1]
    }
    qQ <- s * ncp/(2*sqrt(2))

    if(log.p) {
        log(pP+qQ + pnorm(-ncp))
    } else {
            pP+qQ + pnorm(-ncp)
    }
}
## pntGST23_1 <- Vectorize(pntGST23_1.1, "ncp")





##
### ==== Gil et al. (2023) -- DOI 10.1111/sapm.12609 (see above)
###      ----------------
##' Lemma 2, p.861 + Numerics, p.880 (102) --- Integrate Phi(.) = pnorm(.) ==> "...Inorm"
##'                            ...... ^^^^
pntInorm.1 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                       integr = c("exp", # trafo e^s = t (102)
                                  "direct"), #  (13)
                       logA = log.p || df >= 286.032165, # ~ empirical boundary
                       ..., # arguments passed to integrate()
                       alt = FALSE,
                       verbose = TRUE)
{
    stop("NOT YET implemented; see  ../Misc/pnt-Gil_etal-2023/ ")
### ../Misc/pnt-Gil_etal-2023/Giletal23_Ixpq.R
### ../Misc/pnt-Gil_etal-2023/Giletal23_FN.R
### ../Misc/pnt-Gil_etal-2023/Giletal23_e102.R  # << for int(-inf, inf; .. e^s ds (102) implementation

    if(!lower.tail) { ## upper tail .. is *not* optimized for all cases
        if(!log.p)
            1- pntInorm.1(t, df, ncp, lower.tail=TRUE,
                          ..., verbose=verbose)
        else { # log.p: log(1 - F(.))
           if(alt) ## use log scale
               ## alternatively, not equivalently -- see FIXME in log.p case below
               log1mexp(pntInorm.1(t, df, ncp, lower.tail=TRUE, log.p=TRUE,
                                   integr=integr, logA=logA, ..., verbose=verbose))
           else log1p(- pntInorm.1(t, df, ncp, lower.tail=TRUE,
                                   integr=integr, logA=logA, ..., verbose=verbose))
        }
    }

    integr <- match.arg(integr)
    n2 <- df/2 # used "everywhere" below
    ## A_n from (13) .. must use log-scale as soon as n=df is "large"
    if(logA)
        lA <- n2*log(n2) - lgamma(n2)
    else
        An <- n2^n2 / gamma(n2)

    switch(integr,
           "exp" = { ## formula (102)
               ## integrand <- function(....) ...... #_______________________________ FIXME _________________
               ## I <- integrate(integrand, -Inf, Inf, ...)
           },
           "direct" = { ## formula (13)
               ## integrand <- function(....) ......
               ## I <- integrate(integrand, 0, Inf, ...)
           },
           stop("invalid integration method 'integr'=",integr))

    if(verbose) { cat("Integration:\n"); str(I, digits.d = 10) }
    if(log.p) {
        lA + log(I$value)
    } else { ## lower.tail = TRUE,  log.p == FALSE
        if(logA) exp(lA + log(I$value))
        else An  * I$value
    }
}
pntInorm <- Vectorize(pntInorm.1, c("t", "df", "ncp"))

## Author:  Viktor Witkovsky
## Title:   A Note on Computing Extreme Tail Probabilities of the Noncentral $T$~Distribution with Large Noncentrality Parameter
## Submitted: June 21, 2013.
## Revised:   September 3, 2013

##== Details; experiments ==>      ../Misc/pnt-Gil_etal-2023/Witkovsky_2013/
##   ~~~~~~~  Paper (+ Boost URL): ../Misc/pnt-Gil_etal-2023/Witkovsky_2013/arXiv-1306.5294v2/
##-- ------------------ R fun+tst  ../Misc/pnt-Gil_etal-2023/Witkovsky_2013/nctcdfVW.R  ("translated from Matlab")
pntVW13 <- function(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                      keepS = FALSE, verbose = FALSE)
{
  # Check input lengths  -- MM FIXME -- allow recycling, notably for 2 scalars and 1 non-scalar
  if ((n <- length(t)) != (nnu <- length(df)) | nnu != (nncp <- length(ncp))) {
      n3 <- sort(c(n,nnu,nncp))
      if(!(identical(n3[1:2], c(1L, 1L)) || (n3[1] == 1 && n3[2] == n3[3])))
        stop("Input vectors t, df, and ncp must be of length 1 or the same length, but are",
             n, ", ", nnu, ", ", nncp)
      ## now ensure all three are recycled to full length(.) == n :
      if(n != n3[3])
          n <- length(t <- rep_len(t, n3[3]))
      if(n != nnu ) df <- rep_len(df,  n)
      if(n != nncp) ncp<- rep_len(ncp, n)
  } # now   n == length(t) == length(df) == length(ncp)

  if(verbose) {
     op <- options(str = strOptions(digits.d = 12))# strict.width = "wrap"
     on.exit(options(op))
     mydisp <- function(x) {
         if(n > 2)
             print(x, digits = 16)
         else
             cat(switch(n,
                        sprintf("%20.16g\n", x), # n == 1
                        sprintf("%20.16g %20.16g\n", x), # n == 2
                        ## otherwise.. should not happen!
                        ))
     }
  }

  chooseTail <- function(x, ncp, isLowerGamma, isLowerTail, neg) {
    large <- (x >= ncp)
    isLowerGamma[ large] <- TRUE
    isLowerGamma[!large] <- FALSE
    if (neg) {
        isLowerTail[ large] <- TRUE
        isLowerTail[!large] <- FALSE
    } else {
        isLowerTail[ large] <- FALSE
        isLowerTail[!large] <- TRUE
    }
    return(list(isLowerGamma = isLowerGamma, isLowerTail = isLowerTail))
  }

  integrate <- function(x, df, ncp, isLowerGamma) {
      ## Auxiliary Functions:  limits(), GKnodes(),
      ##' called once from integrate():
      limits <- function(x, df, ncp, NF) { # (x, df, ncp) : num. vector of length NF
        const <- -1.612085713764618
        logRelTolBnd <- -1.441746135564686e+02
        A <- B <- MOD <- numeric(NF)
        ## incpt <- c(1, 1, 1)
        for (i in 1:NF) {
          X <- x[i]
          ##              explicit recycling (no longer needed):
          NU  <- df [i] # [1L + (i-1L) %% length(df)]
          NCP <- ncp[i] # [1L + (i-1L) %% length(ncp)]
          X2 <- X*X
          MOD[i] <- (X * sqrt(4*NU * max(1, NU - 2) + X2 * max(4, NCP*NCP + 4*NU - 8)) -
                     NCP * (X2 + 2*NU)
                    ) / (2 * (X2 + NU))
          dZ <- min(0.5 * abs(MOD[i] + NCP), 0.01)
          dMOD <- c(-dZ, 0, dZ) + MOD[i]
          q. <- (dMOD + NCP)^2 / X2 # q. = q / NU ('q' previously)
          dM.2 <- dMOD^2
          fMOD <- const + 0.5 * (max(1, (NU - 2)) * log(q.) + NU - q.*NU - dM.2)
          ## abc <- lm(fMOD ~ dMOD + I(dMOD * dMOD))$coefficients
          ## ##     ^^     ^^^
          ## abc <- rev(abc) ## chatGPT  got this wrong, too !
          ## Faster than the above; via  "correct" (same matlab/octave code) order of X-columns:
          abc <- .lm.fit(x = cbind(dM.2, dMOD, 1, deparse.level=0L), y = fMOD)$coefficients
          if(verbose >= 2) { cat("inside limits(): abc[1:3] = \n"); print(abc) }

          logAbsTol <- fMOD[2] + logRelTolBnd
          ## Integration limits A and B | Solution to the quadratic EQN:
          ## a*Z^2 + b*Z + c = logAbsTol  where Z ~ dMOD
          D <- sqrt(abc[2] * abc[2] - 4 * abc[1] * (abc[3] - logAbsTol))
          C <- 2 * abc[1]
          if (!is.nan(D)) {
            A[i] <- max(-NCP, -abc[2] / C + D / C)
            B[i] <- max(A[i], -abc[2] / C - D / C)
          } else {
            A[i] <- max(-NCP, MOD[i])
            B[i] <- A[i]
          }
        }

        ## return
        cbind(A = A, B = B, MOD = MOD)
      } ## end  limits()

      integrand <- function(z, x, nu, ncp, isLowerGamma) {
          ##    Z is assumed to be a row vector,
          ##    X, NU, NCP are assumed to be column vectors of equal size
          ##    isLowerGamma is logical vector (indicator of GAMMAINC tail)
          ##    of equal size as X, NU, NCP.
          ##
          ##   In particular, FUNC evaluates the following function
          ##         F = gammainc(Q/2, NU/2, tail) * exp(-0.5 * Z^2) / sqrt(2*pi),
          ##   where
          ##	     Q = NU * (Z + NCP)^2 / X^2.

        const <- 0.398942280401432677939946 ## == phi(0), in R dnorm(0) = 1/sqrt(2*pi)
        FV <- array(0, dim = (d <- dim(z))) # chatGPT got this wrong !
        n <- length(x)
        stopifnot(identical(d[1], n), length(nu) == n, length(ncp) == n, length(isLowerGamma) == n)
        ## NB dim(z) = (n, 240) = (n, 16*15);
        if(verbose >= 3) {
            cat("inside integrand(), before sweep():  str(z+ncp):\n"); str(z+ncp)
            cat("str(x):\n"); str(x)
        }
        nuH <- nu / 2
        q <- sweep(z + ncp, 1L, x, FUN = "/")^2 * nuH
        ##                  --- chatGPT wrongly had '2' here ...
        ## Matlab/octave :  q = bsxfun(@times, bsxfun(@plus,z,ncp).^2, nu./(2*x.*x));
        ## df <- 0.5 * nu  -- use 'nuH = nu/2' below
        z2 <- -0.5 * z^2

        if (any(isL <- isLowerGamma)) {
            FV[isL,] <- pgamma(q[isL,], shape = nuH, lower.tail = TRUE) * exp(z2[isL,])
        }

        if (any(notL <- !isLowerGamma)) {
            FV[notL,] <- pgamma(q[notL,], shape = nuH, lower.tail = FALSE) * exp(z2[notL,])
        }
        ## return
        const * FV
      } ## integrand()

      ## GETSUBS (Sub-division of the integration intervals)
      ## -------- [SUBS,Z,M] = GetSubs(SUBS,XK,G,NK)

      ## Sub-division of the integration intervals for adaptive Gauss-Kronod quadrature
      GetSubs <- function(SUBS, XK, jG, simple=FALSE) {
          ## first call:  SUBS is  2 x 2 matrix (after fixing ChatGPT4o's bug)
        if(verbose && !simple) {
            cat("GetSubs(): input SUBS = \n"); print(SUBS)
        }
        M <- 0.5 * (SUBS[2,] - SUBS[1,]) # length 2 {Matlab: 1 x 2}
        C <- 0.5 * (SUBS[2,] + SUBS[1,]) #   "    2    "       "
        ## Matlab:  Z = XK * M + ones(NK,1)*C;  {MM: 15 x 2 + 15 x 2} -->  15 x 2 matrix
        ##  MM:        {15 x 1} * {1 x 2}
        ## ChatGPT "Quatsch":  Z <- sweep(XK, 2, M, FUN = "*") + C
        ## MM:
        NK <- length(XK)
        ## MM: correct is one of
        ## Z <- XK %o% M + matrix(C, NK, 2, byrow=TRUE)
        ## or
        ## Z <- XK %o% M + rep(C, each = NK)
        ## or
        Z <- outer(XK, M) + rep(C, each = NK)
        if(simple) {
            list(M = M, Z = Z)
        } else {   ## (on first call)   ;; if (nargout > 0) { <-- another GPT bug
          ZjG <- Z[jG, ]
          L <- rbind(SUBS[1,], ZjG)
          U <- rbind(ZjG, SUBS[2,])
          if(verbose) {
              ## cat("GetSubs(*, simple=FALSE): L = \n"); print(L)
              ## cat("           U = \n"); print(U)
              cat("GetSubs(*, simple=FALSE): --> t(SUBS) = \n"); print(cbind(L = c(L), U = c(U)))
          }
          ## return SUBS =
          rbind(c(L), c(U)) # row [1,] = L(ower)
                           ## row [2,] = U(pper)
        }
      } ## end GetSubs() ---


      ##  Begin { integrate() } -------------------------------------------------------

      ## The (G7,K15)-Gauss-Kronod (non-adaptive) quadrature nodes & weights :
      ## GKnodes <- function() {
        nodes <- c(0.2077849550078984676006894, 0.4058451513773971669066064,
                   0.5860872354676911302941448, 0.7415311855993944398638648,
                   0.8648644233597690727897128, 0.9491079123427585245261897,
                   0.9914553711208126392068547)
        wt <- c(0.2044329400752988924141620, 0.1903505780647854099132564,
                0.1690047266392679028265834, 0.1406532597155259187451896,
                0.1047900103222501838398763, 0.0630920926299785532907007,
                0.0229353220105292249637320)
        wt7 <- c(0.3818300505051189449503698, 0.2797053914892766679014678,
                 0.1294849661688696932706114)

        XK <- c(-rev(nodes), 0, nodes)
        WK <- c(rev(wt),  0.2094821410847278280129992, wt)
        WG <- c(rev(wt7), 0.4179591836734693877551020, wt7)
        jG <- seq(2, 15, by = 2) # = 2  4  6  8 10 12 14  -- was 'G' in Witkovsky's code

      rm(nodes, wt, wt7)
      ##   list(XK = XK, WK = WK, WG = WG, G = G)
      ## }

      ## Now the original begin of integrate() - - - - - - - - - - - - - - - - - - - - - -
      stopifnot((NF <- length(x)) >= 1)
      NSUB <- 16 # [N]umber of Gauss-K. [SUB]intervals
      NK <- 15   # [N]umber of Gauss-Konrod [K]nots
      NZ <- NSUB * NK # [N]umber of integrand evaluation points (for each x),  here 240
      z     <- matrix(0, nrow = NF, ncol = NZ)
      halfw <- matrix(0, nrow = NF, ncol = NSUB)
      CDF <- numeric(NF)

      lims <- limits(x, df, ncp, NF)
      ##      ------
      A <- lims[,"A"]
      B <- lims[,"B"]
      MOD <- lims[,"MOD"]
      if(verbose >= 2) { cat("lims <- limits() { = (A, B, MOD)}:\n"); print(lims, digits = 15) }

      ## First, evaluate CDF outside the integration limits [A,B], i.e.
      ## CDF = NORMCDF(A) (for 'upper' gamma tail)
      if (any(upG <- !isLowerGamma)) CDF[upG] <- pnorm(A[upG])
      ## other CDF[.] remain = 0

      ## Second, evaluate CDF inside the integration limits [A,B] by
      ## (G7,K15)-Gauss-Kronod (non-adaptive) quadrature over NSUB=16 subintervals
      ## GKs <- GKnodes()
      ## XK <- GKs$XK # XK, WK are [1:NK]
      ## WK <- GKs$WK
      ## WG <- GKs$WG # WG, jG are [1:7] ((1:(NK/2)))
      ## jG <- GKs$G
      if(verbose >= 2) {
          cat("GKnodes(): [XK WK] // [WG jG] :\n")
          ## print(as.data.frame(GKs[c("XK", "WK")]))
          ## print(as.data.frame(GKs[c("WG", "G" )]))
          print(data.frame(XK, WK))
          print(data.frame(WG, jG))
      }

      for (id in 1:NF) {
          Subs <- GetSubs(rbind(c(  A[id], MOD[id]),
                                c(MOD[id],   B[id])), XK, jG = jG) # 1st call (simple=FALSE)
          SubsL <- GetSubs(Subs, XK, jG = jG, simple=TRUE)         # 2nd call  simple=TRUE
          z    [id,] <- SubsL $ Z
          halfw[id,] <- SubsL $ M
      }
      if(verbose >= 2 && (NF <= 2 || verbose >= 3)) {
          cat("in integrate() after GetSubs():  z = \n"); mydisp(z) # str(z, vec.len = 20) # print(dim(z)); print(z)
          cat(" ...........  halfw[] = M = \n"); mydisp(halfw)
      }

    FV <- integrand(z, x, df, ncp, isLowerGamma)
          ##=======
      if(verbose >= 2) {
          cat(" .. FV <- integrand(): \n") # FV:  length(x) x 240  matrix
          str(FV, vec.len = 20)
          if(NF <= 3) mydisp(FV)
          else { cat("FV[1:3, 1:6]:\n"); print(FV, digits= 10) }
      }

    Q1 <- matrix(0, nrow = NF, ncol = NSUB)
    Q2 <- matrix(0, nrow = NF, ncol = NSUB)

    for (id in 1:NF) {
      F <- matrix(FV[id,], nrow = NK, ncol = NSUB)
      if(verbose) { cat("F ~ FV[id,]:\n"); str(F); cat("F[1:3, 1:6] :"); print(F[1:3, 1:6]) }
      Q1[id,] <- halfw[id,] * colSums(WK * F)
      Q2[id,] <- halfw[id,] * colSums(WG * F[jG, , drop=FALSE]) # <-- drop=FALSE !
      ## FIXME (faster)? Q2[id,] <- halfw[id,] *     sum(WG * F[jG,])
    }
    if(verbose) {
        cat(" .. after summation:  Q1, Q2 :\n")
        cat("Q1: "); mydisp(Q1)
        cat("Q2: "); if(NF==1) print(Q2) else str(Q2)
    }

    CDF <- CDF + rowSums(Q1)
    ErrBnd <- rowSums(abs(Q1 - Q2))

    list(CDF = CDF, ErrBnd = ErrBnd, A = A, B = B, Z = z, F = FV)
  } ## end integrate()

    ## Begin {body of} nctcdfVW(...) : ------------------------------------------------

    ## Initialize outputs
    cdfLower <- cdfUpper <- cdf <- numeric(n)

    isLowerTail <- isLowerGamma <- todo <- rep(TRUE, length(t))

    ## specialcases <- function(x, df, ncp, cdf, isLowerTail, todo) {
    nuneg <- df <= 0
    if (any(nuneg)) {
        cdf[nuneg] <- NA
        isLowerTail[nuneg] <- TRUE
        todo[nuneg] <- FALSE
    }

    if (any(idx <- (t == Inf & !nuneg))) {
        cdf[idx] <- 1
        isLowerTail[idx] <- FALSE
        todo[idx] <- FALSE
    }

    if (any(idx <- (t == -Inf & !nuneg))) {
        cdf[idx] <- 0
        isLowerTail[idx] <- TRUE
        todo[idx] <- FALSE
    }

    xzero <- (t == 0)
    ncppos <- (ncp >= 0)
    if (any(idx <- (xzero & !nuneg & ncppos))) { # t = 0  &  df > 0  &  ncp >= 0
        cdf[idx] <- pnorm(-ncp[idx])
        isLowerTail[idx] <- TRUE
        todo[idx] <- FALSE
    }

    if (any(idx <-(xzero & !nuneg & !ncppos))) { # t = 0  &  df > 0  &  ncp < 0
        cdf[idx] <- pnorm(ncp[idx])
        isLowerTail[idx] <- FALSE
        todo[idx] <- FALSE
    }
    ## end  { special cases }

    if (any(todo)) {
        neg <- (todo & (t < 0))
        pos <- (todo & (t >= 0))

        if (any(neg)) { ## modify {t, ncp, isLowerGamma, isLowerTail} [neg] :
            t  [neg] <- -t  [neg]
            ncp[neg] <- -ncp[neg]
            choose_result <- chooseTail(t[neg], ncp[neg],
                                        isLowerGamma[neg], isLowerTail[neg], TRUE)
            isLowerGamma[neg] <- choose_result$isLowerGamma
            isLowerTail [neg] <- choose_result$isLowerTail
        }

        if (any(pos)) { ## modify {t, ncp, isLowerGamma, isLowerTail} [pos] :
            choose_result <- chooseTail(t[pos], ncp[pos],
                                        isLowerGamma[pos], isLowerTail[pos], FALSE)
            isLowerGamma[pos] <- choose_result$isLowerGamma
            isLowerTail [pos] <- choose_result$isLowerTail
        }

    integrate_result <- integrate(t[todo], df[todo], ncp[todo], isLowerGamma[todo])
    ##                  ^^^^^^^^^
    cdf[todo] <- integrate_result$CDF
    ErrBnd <- integrate_result$ErrBnd
    A <- integrate_result$A
    B <- integrate_result$B
    Z <- integrate_result$Z
    F <- integrate_result$F
  } ## any(todo)

  if(keepS)
    details <- list(t = t, df = df, ncp = ncp,
                    cdf_computed = cdf,
                    tail_computed = isLowerTail,
                    gamma_lower = isLowerGamma,
                    error_upperBound = ErrBnd,
                    integration_limits = cbind(A, B),
                    integrand_abscissa = Z, # F[i,]  vs  Z[i,] -- to be plotted
                    integrand_values   = F) # (they *where* plotted in V.W.'s Matlab code)

  if(log.p) {
      cdfLower[ isLowerTail] <-     log(cdf[ isLowerTail])
      cdfLower[!isLowerTail] <- log1p(- cdf[!isLowerTail])
      cdfUpper[!isLowerTail] <-     log(cdf[!isLowerTail])
      cdfUpper[ isLowerTail] <- log1p(- cdf[ isLowerTail])
  } else {
      ## MM: I think the following is "inherently doubtful": '1 - cdf[..]' almost always loses some accuracy
      cdfLower[ isLowerTail] <-     cdf[ isLowerTail]
      cdfLower[!isLowerTail] <- 1 - cdf[!isLowerTail]
      cdfUpper[!isLowerTail] <-     cdf[!isLowerTail]
      cdfUpper[ isLowerTail] <- 1 - cdf[ isLowerTail]
  }
  ## MM: (read above) ==> only returning cdf  is "sub optimal"
  cdf <- if(lower.tail) cdfLower else cdfUpper
  ## return (possibly with details):
  if(keepS) {
      details$cdf_lowerTail <- cdfLower
      details$cdf_upperTail <- cdfUpper
      list(cdf = cdf, details = details)
  } else
      cdf
} ## end pntVW13



###--- dnt() --- for the non-central *density*
###                      =====================
### Wikipedia (or other online sources) have no formula for the density

### Johnson, Kotz and Balakrishnan (1995) [2nd ed.] have
###  (31.15) [p.516] and (31.15'), p.519 -- and they  contradict by a factor  (1 / j!)
###  in the sum ---> see below


if(FALSE) {## Just for historical reference :
           ## =============================
## was  dntR.1 <-
dntRwrong1 <- function(x, df, ncp, log = FALSE, M = 1000, check=FALSE, tol.check = 1e-7)
{
    ## R's source ~/R/D/r-devel/R/src/nmath/dnt.c  claims -- from 2003 till 2014 -- but *WRONGLY*
    ## *	   f(x, df, ncp) =
    ## *		df^(df/2) * exp(-.5*ncp^2) /
    ## *		(sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2)) *
    ## *		sum_{k=0}^Inf  gamma((df + k + df)/2)*ncp^k /
    ## *				prod(1:k)*(2*x^2/(df+x^2))^(k/2)
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              ncp >= 0, df > 0,
              is.numeric(M), M >= 1, M == round(M))

    if(check)
     fac <- df^(df/2) * exp(-.5*ncp^2) / (sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2))
    lfac <- (df/2)*log(df) -.5*ncp^2 - (log(pi)/2 + lgamma(df/2)  +(df+1)/2* log(df+x^2))
    if(check) stopifnot(all.equal(fac, exp(lfac), tol=tol.check))
    k <- 0:M
    if(check) suppressWarnings(
    terms  <-  gamma((df + k + df)/2) * ncp^k    /  factorial(k) * (2*x^2/(df+x^2)) ^ (k/2)  )
    lterms <- lgamma((df + k + df)/2)+k*log(ncp) - lfactorial(k) + (k/2) * log(2*x^2/(df+x^2))
    if(check) {
        ii <- which(!is.na(terms))
        stopifnot(all.equal(exp(lterms[ii]), terms[ii], tol=tol.check))
    }
    lf <- lfac + lsum(lterms)
    if(log) lf else exp(lf)
}
dntRwrong <- Vectorize(dntRwrong1, c("x", "df", "ncp"))
}##-- just for historical reference

### Johnson, Kotz and Balakrishnan (1995) [2nd ed.] have
###  (31.15) [p.516] and (31.15'), p.519 -- and they  contradict by a factor  (1 / j!)
###  in the sum
### ==> via the following: "prove" that (31.15) is correct, and  (31.15') is missing the "j!"
## was *old*  'dnt.1'
.dntJKBch1 <- function(x, df, ncp, log = FALSE, M = 1000, check=FALSE, tol.check = 1e-7)
{
    stopifnot(length(x) == 1, length(df) == 1, length(ncp) == 1, length(M) == 1,
              ncp >= 0, df > 0,
              is.numeric(M), M >= 1, M == round(M))
    if(check) ## this is the formula (31.15) unchanged:
     fac <- exp(-.5*ncp^2) * gamma((df+1)/2) / (  sqrt(pi*df)* gamma(df/2)) *    (df/(df+x^2))^((df+1)/2)
    lfac <-     -.5*ncp^2  +lgamma((df+1)/2) - (.5*log(pi*df)+lgamma(df/2)) + log(df/(df+x^2))*((df+1)/2)
    if(check) stopifnot(all.equal(fac, exp(lfac), tol=tol.check))
    j <- 0:M
    if(check)  suppressWarnings(## this is the formula (31.15) unchanged:
        ## [note that 1/gamma((df+1)/2) is indep. of j!]
    terms  <-  gamma((df + j + 1)/2) / ( factorial(j) * gamma((df + 1)/2))*    (x*ncp*sqrt(2)/sqrt(df+x^2))^ j )
    lterms <- lgamma((df + j + 1)/2) - (lfactorial(j) +lgamma((df + 1)/2))+ log(x*ncp*sqrt(2)/sqrt(df+x^2))* j
    if(check) {
        ii <- which(!is.na(terms))
        stopifnot(all.equal(exp(lterms[ii]), terms[ii], tol=tol.check))
    }
    lf <- lfac + lsum(lterms)
    if(log) lf else exp(lf)
}
.dntJKBch <- Vectorize(.dntJKBch1, c("x", "df", "ncp"))

## Q{MM}: is the [ a < x ] cutoff really exactly optimal?
logr <- function(x, a) { ## == log(x / (x + a)) -- but numerically smart; x >= 0, a > -x
    if(length(aS <- a < x) == 1L) {
        if(aS) -log1p(a/x) else log(x / (x + a))
    } else { # "vectors" : do ifelse(aS, .., ..) efficiently:
        r <- a+x # of correct (recycled) length and type (numeric, mpfr,  ..)
        r[ aS] <- -log1p((a/x)[aS])
        r[!aS] <- log((x / r)[!aS])
        r
    }
}

## New "optimized" and  "mpfr-aware" and *vectorized* (!) version:
## was 'dnt.1'
dntJKBf <- function(x, df, ncp, log = FALSE, M = 1000)
{
    stopifnot(length(M) == 1, df >= 0, is.numeric(M), M >= 1, M == round(M))
    ln2 <- log(2)
    ._1.1..M <- c(1L, seq_len(M)) # cumprod(.) = (0!, 1!, 2! ..) =  (1, 1, 2, 6, ...)
    isN <- is.numeric(x) && is.numeric(df) && is.numeric(ncp)
    isMpfr <- !isN && any_mpfr(x, df, ncp)
    if(!isN) {
        if(isMpfr) {
	    stopifnot(requireNamespace("Rmpfr"))
            mpfr <- Rmpfr::mpfr ; getPrec <- Rmpfr::getPrec
	    prec <- max(getPrec(x), getPrec(df), getPrec(ncp))
	    pi <- Rmpfr::Const("pi", prec = max(64, prec))
	    ## this *is* necessary / improving results!
	    if(!inherits(x,  "mpfr")) x   <- mpfr(x,  prec)
	    if(!inherits(df, "mpfr")) df  <- mpfr(df, prec)
	    if(!inherits(ncp,"mpfr")) ncp <- mpfr(ncp,prec)
	    ln2 <- log(mpfr(2, prec))
	    ._1.1..M <- mpfr(._1.1..M, prec)
        } else {
            warning(" Not 'numeric' but also not  'mpfr' -- untested, beware!!")
        }
        ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
	## .N <- if(isMpfr) Rmpfr::asNumeric else as.numeric
    }

    x2 <- x^2
    lfac <- -ncp^2/2 - (.5*log(pi*df)+lgamma(df/2)) + logr(df, x2)*(df+1)/2
    j <- 0:M
    lfact.j <- cumsum(log(._1.1..M)) ## == lfactorial(j)
    nd <- length(df)
    LogRt <- 2*log(abs(ncp)) + ln2 + logr(x2, df) # (full length)
    ## now vectorize "lSum(x,df,ncp)" :
    lSum <- dx <- ncp*x + 0*df # delta * x  [of full length],  (correct == 0 for dx == 0)
    if(any(negD <- dx < 0)) # if(negD[i]) "alternating sum"
	sigPM <- rep_len(c(1,-1), length(j))
    for(i in which(dx != 0)) { # --- compute lSum[i] ----------------
	## use abs(ncp) : if(ncp < 0) ncp <- -ncp
	##lterms <- lgamma((df+j + 1)/2) - lfact.j + log(x*ncp*sqrt(2)/sqrt(df+x^2))* j
	lterms <- lgamma((df[1L+ (i-1L)%% nd] + j + 1)/2) - lfact.j + LogRt[i] * j/2
        lSum[i] <-
            if(negD[i]) ## this is hard: even have *negative* sum {before log(.)} with mpfr,
                ## e.g. in  dnt.1(mpfr(-4, 128), 5, 10)   ???
                lssum(lterms, signs = sigPM, strict=FALSE)
            else
                lsum(lterms)
    }
    lf <- lfac + lSum
    if(log) lf else exp(lf)
}
## No longer, as have vectorized above!
## dntJKBf <- Vectorize(dntJKBf1, c("x", "df", "ncp"))
## instead, from 2019-10-04 :
dntJKBf1 <- function(x, df, ncp, log = FALSE, M = 1000) {
    .Deprecated("dntJKBf")
    dntJKBf(x=x, df=df, ncp=ncp, log=log, M=M)
}


## Orig: ~/R/MM/NUMERICS/dpq-functions/noncentral-t-density-approx_WV.R
##
## From: Wolfgang Viechtbauer <wviechtb@s.psych.uiuc.edu>
## To: <r-help@stat.math.ethz.ch>
## Subject: Re: [R] Non-central distributions
## Date: Fri, 18 Oct 2002 11:09:57 -0500 (CDT)
## .......
## This is an approximation based on Resnikoff & Lieberman (1957).
## .. quite accurate. .....
## .......
## .......
##
## MM: added 'log' argument and implemented  log=TRUE
## --- TODO: almost untested by MM [but see Wolfgang's notes in *_WV.R (s.above)]
dtWV <- function(x, df, ncp=0, log=FALSE) {
   dfx2 <- df + x^2 # = 'f+t^2' in Resnikoff+L.(1957), p.1 (by MM)
   y <- -ncp*x/sqrt(dfx2) # = 'y' in R.+L., p.1
   ## MM(FIXME): cancellation for y >> df  here :
   a <- (-y + sqrt(y^2 + 4*df)) / 2 # NB a = 't' in R.+L., p.25
   dfa2 <- df+a^2 ## << MM(2)
   if(log) {
       lHhmy <- df*log(a) + -0.5*(a+y)^2 +
           0.5*log(2*pi*a^2/dfa2) +
           log1p( - 3*df/(4*dfa2^2) + 5*df^2/(6*dfa2^3))
       lHhmy - (((df-1)/2)*log(2) + lgamma(df/2) + .5*log(pi*df)) +
           -0.5*df*ncp^2/dfx2 + ((df+1)/2)*log(df/dfx2)
   } else { ## MM: cancelled 1/f! = 1/gamma(df+1) in Hh_f(y) =: Hhmy : formula p.25
       Hhmy <- a^df * exp(-0.5*(a+y)^2) *
           sqrt(2*pi*a^2/dfa2) *
           (1 - 3*df/(4*dfa2^2) + 5*df^2/(6*dfa2^3))
       ## formula p.1:  h(f,δ,t) = (....) * Hh_f(-δ t / sqrt(f+t²)) = (....) * Hhmy
       Hhmy / (2^((df-1)/2) * gamma(df/2) * sqrt(pi*df)) *
           exp(-0.5*df*ncp^2/dfx2) * (df/dfx2)^((df+1)/2)
   }
}


###-- qnt() {C} i.e. qt(.., ncp=*) did not exist yet at the time I wrote this ...
##    ---
## was  qt.appr <-
qtAppr <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE,
                    method = c("a","b","c"))
{
  ## Purpose: Quantiles of approximate non-central t
  ##  using Johnson,Kotz,.. p.521, formula (31.26 a) (31.26 b) & (31.26 c)
  ## ----------------------------------------------------------------------
  ## Arguments: see  ?qt
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date:  6 Feb 99
    method <- match.arg(method)

  ##----------- NEED df >> 1 (this is from experiments below; what exactly??)
###___ TODO: Compare with the comment in qtR1() [ = R's src/nmath/qt.c ] about the df >= 20 qnorm() approx <<<<<<<<<<<<<
###___ ====  which is said to be  like  Abramowitz & Stegun 26.7.5 (p.949)
  z <- qnorm(p, lower.tail=lower.tail, log.p=log.p)
  if(method %in% c("a","c")) {
      b <- b_chi(df)
      b2 <- b*b
  }
  ## For huge `df';  b2 ~= 1 ---> method b) sets b = b2 = 1
  switch(method,
         "a" = {
             den <- b2 - z*z*(1-b2)
             (ncp*b + z*sqrt(den + ncp^2*(1-b2)))/den
         },
         "b" = {
             den <- 1 - z*z/(2*df)
             (ncp + z*sqrt(den + ncp^2/(2*df)))/den

         },
         "c" = {
             den <- b2 - z*z/(2*df)
             (ncp*b + z*sqrt(den + ncp^2/(2*df)))/den
         })
}

##
##' Pure R version of R's C function in  ..../src/nmath/qnt.c
##' additionally providing "tuning" consts (accu, eps)
##'
##' qntR1: for length-1 arguments

## These are the "old" versions,
qntRo1 <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE,
                  pnt = stats::pt, accu = 1e-13, eps = 1e-11)
{
    stopifnot(length(p) == 1L, length(df) == 1L, length(ncp) == 1L
            , accu >= 0, eps >= 0
            , is.function(pnt), names(formals(pnt)) >= 5
              )
    if(is.na(p) || is.na(df) || is.na(ncp))
        return(p + df + ncp)
    if (df <= 0.0) { ## ML_WARN_return_NAN;
        warning("Non-positive 'df': NaNs produced")
        return(NaN)
    }
    if(ncp == 0. && df >= 1.)
        return(qt(p, df, lower.tail, log.p))
    ## for  ncp=0 and df < 1  continue :

    ## R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF) :
    if(p == .D_0(log.p)) return(if(lower.tail) -Inf else  Inf)
    if(p == .D_1(log.p)) return(if(lower.tail)  Inf else -Inf)
    if(p < .D_0(log.p) ||
       p > .D_1(log.p)) { warning("p out of range"); return(NaN) }

    if (!is.finite(df)) # df = Inf ==> is limit = N(ncp,1) :
	return(qnorm(p, ncp, 1., lower.tail, log.p))

    ## MM: this *does* lose accuracy in case of very small p
    ## --- TODO: do *not* do this, but s/ (TRUE, FALSE) by (lower.tail, log.p) below
    ##           *and* "transform" (lx, ux) accordingly, for lower.tail swap "<" with ">"
    p <- .DT_qIv(p, lower.tail, log.p)

    ##* Invert pnt(.) :
    ##  -------------
    ##* 1. finding an upper and lower bound
    Mdeps <- .Machine$double.eps
    if(p > 1 - Mdeps) return(Inf)
    pp <- min(1 - Mdeps, p * (1 + eps))
    ux <- max(1., ncp)
    DBL_MAX <- .Machine$double.xmax
    while(ux < DBL_MAX && pnt(ux, df, ncp, TRUE, FALSE) < pp)
	ux <- ldexp(ux, 1L)
    pp <- p * (1 - eps)
    lx <- min(-1., -ncp)
    while(lx > -DBL_MAX && pnt(lx, df, ncp, TRUE, FALSE) > pp)
	lx <- ldexp(lx, 1L) ## * 2

    ##* 2. interval (lx,ux)  halving :
    repeat {
	nx <- 0.5 * (lx + ux) ## could be zero
	if(pnt(nx, df, ncp, TRUE, FALSE) > p) ux <- nx else lx <- nx
        ## while(...) :
        if(!((ux - lx) > accu * max(abs(lx), abs(ux)))) break
    }

    ## return
    ldexp(lx + ux, -1L) # = 0.5 * (lx + ux)
}

qntRo <- Vectorize(qntRo1, c("p", "df", "ncp"))

## These are new versions which do *not* transform p to "usual" scale:
qntR1 <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE,
                  pnt = stats::pt, accu = 1e-13, eps = 1e-11)
{
    stopifnot(length(p) == 1L, length(df) == 1L, length(ncp) == 1L
            , accu >= 0, eps >= 0
            , is.function(pnt), names(formals(pnt)) >= 5
              )
    if(is.na(p) || is.na(df) || is.na(ncp))
        return(p + df + ncp)
    if (df <= 0.0) { ## ML_WARN_return_NAN;
        warning("Non-positive 'df': NaNs produced")
        return(NaN)
    }
    if(ncp == 0. && df >= 1.)
        return(qt(p, df, lower.tail, log.p))
    ## for  ncp=0 and df < 1  continue :

    ## R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF) :
    if(p == .D_0(log.p)) return(if(lower.tail) -Inf else  Inf)
    if(p == .D_1(log.p)) return(if(lower.tail)  Inf else -Inf)
    if(p < .D_0(log.p) ||
       p > .D_1(log.p)) { warning("p out of range"); return(NaN) }

    if (!is.finite(df)) # df = Inf ==> is limit = N(ncp,1) :
	return(qnorm(p, ncp, 1., lower.tail, log.p))

    ## MM: this *does* lose accuracy in case of very small p
    ## --- TODO: do *not* do this, but s/ (TRUE, FALSE) by (lower.tail, log.p) below
    ##           *and* "transform" (lx, ux) accordingly, for lower.tail swap "<" with ">"
    ## p_ <- .DT_qIv(p, lower.tail, log.p)

    ##* Invert pnt(.) :
    ##  -------------
    ##* 1. finding an upper and lower bound
    Mdeps <- .Machine$double.eps
    ## if(p > 1 - Mdeps) return(Inf)
    pp <- if(log.p) min(log1p(-Mdeps), p + eps) else min(1 - Mdeps, p * (1 + eps))
    ux <- max(1., ncp)
    DBL_MAX <- .Machine$double.xmax
    MAX2 <- DBL_MAX / 2.
    while(ux < MAX2 && pnt(ux, df, ncp, lower.tail, log.p) < pp)
	ux <- ldexp(ux, 1L)
    pp <-  if(log.p) p - eps else p * (1 - eps)
    lx <- min(-1., -ncp)
    while(lx > -MAX2 && pnt(lx, df, ncp, lower.tail, log.p) > pp)
	lx <- ldexp(lx, 1L) ## * 2

    ##* 2. interval (lx,ux)  halving :
    repeat {
	nx <- 0.5 * (lx + ux) ## could be zero
	if(pnt(nx, df, ncp,  lower.tail, log.p) > p) ux <- nx else lx <- nx
        ## while(...) :
        if(!((ux - lx) > accu * max(abs(lx), abs(ux)))) break
    }

    ## return
    ldexp(lx + ux, -1L) # = 0.5 * (lx + ux)
}

qntR <- Vectorize(qntR1, c("p", "df", "ncp"))

## Invert pnt() via uniroot() -- is currently more reliable
qtU1 <- function(p, df, ncp, lower.tail=TRUE, log.p=FALSE,
                 interval = c(-10,10), ## FIXME! use qtAppr(), possibly simple 'method = "b"' ?
                 tol = 1e-5, verbose = FALSE, ...)
{
    ## Be smart when missing(ncp) :  "keep it" missing also in the pt() call
    stopifnot(length(p) == 1, length(df) == 1, missing(ncp) || length(ncp) == 1)
    ptFUN <- if(verbose) { # verbose: print(.)
                 function(q) {
                     p. <- pt(q, df, ncp, lower.tail, log.p)
                     d <- p. - p
                     cat(sprintf("q=%11g --> pt(q,*) = %11g, d = pt() - p = %g\n",
                                 q, p., d))
                     d
                 }
             } else if(missing(ncp)) {
                 function(q) pt(q, df,    , lower.tail, log.p) - p
             } else
                 function(q) pt(q, df, ncp, lower.tail, log.p) - p

    if(verbose && missing(ncp)) ## replace 'ncp' by <empty> in the pt() call:
        body(ptFUN)[[2:4]] <- alist(x=)$x
    r <- uniroot(ptFUN,
                 interval = interval,
                 extendInt = if(lower.tail) "upX" else "downX",
                 tol = tol, ...)
    if(verbose) { cat("uniroot(.):  "); str(r) }
    r$root
}

qtU <- Vectorize(qtU1, c("p","df","ncp"))
