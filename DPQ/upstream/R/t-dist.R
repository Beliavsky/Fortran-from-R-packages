#### t-distribution functionality
###
### NB: *non*-central t   is in  >>>>>  ./t-nonc-fn.R <<<<

## buglets in qt() ==> PR#18360 --  https://bugs.r-project.org/show_bug.cgi?id=18360
## --- look at a pure R version of R's own qt() in  src/nmath/qt.c

 ## Mathlib : A C Library of Special Functions
 ## Copyright (C) 2000-2023 The R Core Team
 ## Copyright (C) 2003-2022 The R Foundation
 ## Copyright (C) 1998 Ross Ihaka

 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 3 of the License, or
 ## (at your option) any later version.

 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.

 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, a copy is available at
 ## https://www.R-project.org/Licenses/

 ## DESCRIPTION

 ##     The "Student" t distribution quantile function.

 ## NOTES

 ##     This is a C translation of the Fortran routine given in:
 ##     Hill, G.W (1970) "Algorithm 396: Student's t-quantiles"
 ##     CACM 13(10), 619-620.

 ##     Supplemented by inversion for 0 < ndf < 1.

 ## ADDITIONS:
 ##     - lower_tail, log_p
 ##     - using	 expm1() : takes care of  Lozy (1979) "Remark on Algo.", TOMS
 ##     - Apply 2-term Taylor expansion as in
 ##       Hill, G.W (1981) "Remark on Algo.396", ACM TOMS 7, 250-1
 ##     - Improve the formula decision for 1 < df < 2

qtR1 <- function(p, df, lower.tail=TRUE, log.p=FALSE, eps = 1e-12,
                 d1_accu = 1e-13, d1_eps = 1e-11,
                 itNewt = 10L, epsNewt = 1e-14, logNewton = log.p,
                 verbose = FALSE)
{
    stopifnot(length(p) == 1, length(df) == 1, df > 0,
              (itNewt <- as.integer(itNewt)) >= 1L)

    if (is.na(p) || is.na(df))
	return(p + df)

    ## R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF) :
    if(p == .D_0(log.p)) return(if(lower.tail) -Inf else  Inf)
    if(p == .D_1(log.p)) return(if(lower.tail)  Inf else -Inf)
    if(p < .D_0(log.p) ||
       p > .D_1(log.p)) { warning("p out of range"); return(NaN) }

    if (df < 0.0) { ## ML_WARN_return_NAN;
        warning("Negative-positive 'df': NaNs produced")
        return(NaN)
    }

    R_FINITE <- is.finite

    if (df < 1) { ## based on qnt() --> see ./t-nonc-fn.R <<<<<<<<<<<<<
	## const static double d1_accu = 1e-13;
	## const static double d1_eps = 1e-11; ## must be > d1_accu */

	## double ux, lx, nx, pp;

	iter <- 0
        ##  ?? (not good)
        if(verbose) cat(sprintf("qt(p, df, *) -- case df < 1: p=%12g, ", p))
	p <- .DT_qIv(p, lower.tail, log.p) #  R_DT_qIv(p)
        if(verbose) cat(sprintf(" new p=%12g\n", p))

	## Invert pt(.) :
	##   1. finding an upper and lower bound */
	if(p > 1 - DBL_EPSILON) return (ML_POSINF);
	pp = min(1 - DBL_EPSILON, p * (1 + d1_eps))
        ux <- 1
	while(ux < DBL_MAX && pt(ux, df, , TRUE, FALSE) < pp) ux <- ux * 2
	pp = p * (1 - d1_eps)
	lx <- -1
        while(lx > -DBL_MAX && pt(lx, df, , TRUE, FALSE) > pp) lx <- lx * 2

        ## 2. interval (lx,ux)  halving
        ##   regula falsi failed on qt(0.1, 0.1)
        repeat {
	    nx = 0.5 * (lx + ux);
	    if (pt(nx, df, , TRUE, FALSE) > p) ux <- nx else lx <- nx
            if(! ((ux - lx) / abs(nx) > d1_accu && (iter <- iter+1L) < 1000))
                break
	}

	if(iter >= 1000) ML_WARNING(ME_PRECISION, "qt");

	return(0.5 * (lx + ux))
    }

    ## Old comment:
    ## FIXME: "This test should depend on  df  AND p  !!
    ## -----  and in fact should be replaced by
    ## something like Abramowitz & Stegun 26.7.5 (p.949)"
    ##
    ## That would say that if the qnorm value is x then
    ## the result is about x + (x^3+x)/4df + (5x^5+16x^3+3x)/96df^2
    ## The differences are tiny even if x ~ 1e5, and qnorm is not
    ## that accurate in the extreme tails.
    ##
    if (df > 1e20) return(qnorm(p, 0., 1., lower.tail, log.p))

    if(verbose) cat(sprintf("qt(p=%11g, df=%11g, *) -- general case\n", p, df))

    P <- .D_qIv(p, log.p) # R_D_qIv(p) -- if exp(p) underflows, we fix below */

    ## Rboolean
    neg = (!lower.tail || P < 0.5) && (lower.tail || P > 0.5)
    is_neg_lower = (lower.tail == neg) # both TRUE or FALSE == !xor */
    if(verbose) cat(sprintf(" -> P=%11g, neg=%s, is_neg_lower=%s;", P,
                            format(neg), format(is_neg_lower)))
    if(neg)
	P = 2 * if(log.p) (if(lower.tail) P else -expm1(p)) else .D_Lval(p, lower.tail)
    else
	P = 2 * if(log.p) (if(lower.tail) -expm1(p) else P) else .D_Cval(p, lower.tail)
    ## 0 <= P <= 1 ; P = 2*min(P', 1 - P')  in all cases */
    if(verbose) cat(sprintf(" -> final P=%11g\n", P))
    if (abs(df - 2) < eps) {	## df ~= 2 */
	if(P > DBL_MIN) {
	    if(3* P < DBL_EPSILON) ## P ~= 0 */
		q = 1 / sqrt(P)
	    else if (P > 0.9)	   ## P ~= 1 */
		q = (1 - P) * sqrt(2 /(P * (2 - P)))
	    else ## eps/3 <= P <= 0.9 */
		q = sqrt(2 / (P * (2 - P)) - 2)
	}
	else { ## P << 1, q = 1/sqrt(P) = ... */
	    if(log.p)
		q = if(is_neg_lower) exp(- p/2) / M_SQRT2  else  1/sqrt(-expm1(p))
	    else
		q = ML_POSINF
	}
    }
    else if (df < 1 + eps) { ## df ~= 1  (df < 1 excluded above): Cauchy */
	if(P == 1.)
            q = 0 ## some versions of tanpi give Inf, some NaN
	else if(P > 0)
	    q = 1/tanpi(P/2.) ## == - tan((P+1) * M_PI_2) -- suffers for P ~= 0 */
	else { ## P = 0, but maybe = 2*exp(p) ! */
	    if(log.p) ## 1/tan(e) ~ 1/e */
		q = if(is_neg_lower) M_1_PI * exp(-p) else -1./(M_PI * expm1(p))
	    else
		q = ML_POSINF
	}
    }
    else {		##-- usual case;  including, e.g.,  df = 1.1 */
        if(verbose) cat("  usual 'df' case: ")

	## double
        x = 0.
        log.p2 = 0.## -Wall */;
        a = 1 / (df - 0.5)
        b = 48 / (a * a)
        c = ((20700 * a / b - 98) * a - 16) * a + 96.36
        d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * df
	## Rboolean
        P_ok1 = P > DBL_MIN || !log.p
        P_ok  = P_ok1 ## when true (after check below), use "normal scale": log.p=FALSE
        if(verbose) cat(sprintf("  P_ok:= P_ok1 = %s, ", format(P_ok1)))
	if(P_ok1) {
	    y = (d * P) ^(2./df)
	    P_ok = (y >= DBL_EPSILON);
            if(verbose) cat(sprintf("  y=%11g, P_ok = %s, ", y, format(P_ok)))
	}
	if(!P_ok) {## log.p && P very.small  ||  (d*P)^(2/df) =: y < eps_c
	    log.p2 <- if(is_neg_lower)  .D_log(p, log.p) else .D_LExp(p, log.p); # == log(P / 2)
	    x = (log(d) + M_LN2 + log.p2) / df;
	    y = exp(2 * x);
            if(verbose) cat(sprintf("  !P_ok: log.p2=%11g, y=%11g\n", log.p2, y))
	}

	if ((df < 2.1 && P > 0.5) || y > 0.05 + a) { ## P > P0(df) */
            if(verbose) cat(" P > P0(df): Asymptotic inverse expansion about normal\n")
	    if(P_ok)
		x = qnorm(0.5 * P, 0., 1., TRUE,  FALSE)
	    else ## log.p && P underflowed */
		x = qnorm(log.p2,  0., 1., lower.tail, TRUE)

	    y = x * x;
	    if (df < 5)
		c <- c + 0.3 * (df - 4.5) * (x + 0.6);
	    c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
	    y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c
		  - y - 3) / b + 1) * x;
            if(verbose) cat(sprintf(" x=qnorm(..)=%12g, c=%12g, new y=%12g\n",
                                    x, c, y))
	    y = expm1(a * y * y);
	    q = sqrt(df * y);
	} else if(!P_ok && x < - M_LN2 * DBL_MANT_DIG) {## log(DBL_EPSILON) */
	    ## y above might have underflown */
	    q = sqrt(df) * exp(-x);
            if(verbose) cat(sprintf("!P_ok && x < log(epsC) = -36.04: q=%12g\n", q))
	}
	else { ## re-use 'y' from above */
	    y = ((1 / (((df + 6) / (df * y) - 0.089 * d - 0.822)
		       * (df + 2) * 3) + 0.5 / (df + 4))
		 * y - 1) * (df + 1) / (df + 2) + 1 / y;
	    q = sqrt(df * y);
            if(verbose) cat(sprintf("re-use y := .. = %12g, q = sqrt(df*y) = %12g\n", y,q))
	}


	## Now apply 2-term Taylor expansion improvement (1-term = Newton):
        ## as by Hill (1981) [ref.above] */

	## FIXME: This can be far from optimal when log.p = TRUE
        ##        but is still needed, e.g. for qt(-2, df=1.01, log=TRUE).
        ##    	  Probably also improvable when  lower.tail = FALSE */
	{ ## was if(P_ok1) {
	    it <- 0
            M <- abs(sqrt(DBL_MAX/2.) - df)
            if(logNewton && log.p) {
              if(verbose) cat("P_ok1: log-scale Taylor (iterated):\n")
	      while((it <- it+1L) <= itNewt && ## (y <- dt(q, df, log=FALSE)) > 0 &&
                    { lF <- pt(q, df, lower.tail=FALSE, log.p=TRUE)
                      R_FINITE(x <- exp(lF - dt(q, df, log=TRUE))*(lF - log(P/2))) } &&
                             ## FIXME:  directly from orig. p (lower.?) ^^^^^^^^ via log1mexp(.) !
                    abs(x) > epsNewt*abs(q)) {
                      ## Newton (=Taylor 1 term): q += x;
                      ## ___TODO___ (?) Taylor 2-term :....
                      if(verbose) cat(sprintf(
			"it=%3d, ... d{q}1=exp(lF - dt(q, df, log=TRUE))*(lF - log(P/2)) = %13g; ",
                                      it, x))
                      if(R_FINITE(qn <- q + x))
                          q <- qn
                      else ##/ FIXME??  if  q+x = +/-Inf is *better* than q should still use it
                          break ##; // cannot improve  q  with a Newton/Taylor step
                      if(verbose) cat(sprintf("new q=%12g\n", q))
                  }
            }
            else { ## log.p or logNewton is false
              if(verbose) cat("P_ok1: 2-step Taylor (iterated):\n")
	      while((it <- it+1L) <= itNewt && (y <- dt(q, df, log=FALSE)) > 0 &&
		  R_FINITE(x <- (pt(q, df, , FALSE, FALSE) - P/2) / y) &&
		  abs(x) > epsNewt*abs(q)) {
                      ## Newton (=Taylor 1 term):
                      ##  q += x;
                      ## Taylor 2-term : */
                      F <- if(abs(q) < M)
                               q * (df + 1) / (2 * (q * q + df))
                           else    (df + 1) / (2 * (q     + df/q))
                      del_q <- x * (1. + x * F)
                      if(verbose) cat(sprintf(
				"it=%3d, y=dt(*)=%13g, d{q}1=(pt(q,*) - P/2)/y = %13g; d{q}2 = %13g; ",
                                it, y, x, del_q))
                      if(R_FINITE(del_q) && R_FINITE(qn <- q + del_q))
                          q <- qn
                      else if(R_FINITE(qn <- q + x)) # have checked R_FINITE(x) already
                          q <- qn
                      else ##/ FIXME??  if  q+x = +/-Inf is *better* than q should still use it
                          break ##; // cannot improve  q  with a Newton/Taylor step
                      if(verbose) cat(sprintf("new q=%12g\n", q))
                  }
            }
            if(verbose && it <= 1L) cat(">> *no* Newton refinements <<\n")
	}
    }
    return (if(neg) -q else q)
}

qtR <- Vectorize(qtR1, c("p", "df"))
## when interactively:  assignInNamespace("qtR", qtR, ns=asNamespace("DPQ"))


## large df Approx. from comment above:
qtNappr <- function(p, df, lower.tail=TRUE, log.p=FALSE, k = 2) {
    stopifnot(k == (k. <- as.integer(k)), 0 <= (k <- k.), k <= 4)
    ## something like Abramowitz & Stegun 26.7.5 (p.949)"
    ##
    ## That would say that if the qnorm value is x then
    ## the result is about x + (x^3+x)/4df + (5x^5+16x^3+3x)/96df^2
    ## The differences are tiny even if x ~ 1e5, and qnorm is not
    ## that accurate in the extreme tails.
    x <- qnorm(p, lower.tail=lower.tail, log.p=log.p)
    switch(k+1
         , x # k=0
         , x*(1+ (x^2+1)/(4*df)) # == x + (x^3+x)/4df  --- k=1

         ##{ x2 <- x^2; x*(1 +(x2+1)/(4*df) + ((5*x2+16)*x2+3)/(96*df^2)) } # --- k=2
         , { x2 <- x^2; x*(1 +(x2+1 + ((5*x2+16)*x2+3)/(24*df))/(4*df)) }   # --- k=2
           ##         = x + (x^3+x)/4df + (5x^5+16x^3+3x)/96df^2

         , { x2 <- x^2; n4 <- 4*df
            ## x*(1 + (x2+1 + ((5*x2+16)*x2+3 + (((3*x2+19)*x2+17)*x2-15)/(4*df))/(24*df)))/(4*df) } # --- k=3
             x * (1 + (x2+1 + ((5*x2+16)*x2+3 + (((3*x2+19)*x2+17)*x2-15)/n4) / (24*df)) / n4) } # --- k=3
           ## = x + (x^3+x)/4df + (5x^5+16x^3+3x)/96df^2 + (3x^7+19x^5+17x^3-15x)/384df^3


         , { x2 <- x^2; n4 <- 4*df
         ## x * (1 + (x2+1 + ((5*x2+16)*x2+3 + (((3*x2+19)*x2+17)*x2-15)/n4)/(24*df))/n4) +
         ## (79*x^9+776*x^7+1482*x^5-1920*x^3-945*x)/(92160*df^4) }
         ## x * (1 + (x2+1 + ((5*x2+16)*x2+3 + (((3*x2+19)*x2+17)*x2-15)/n4)/(24*df))/n4 +
         ##      ((((79*x2+776)*x2+1482)*x2-1920)*x2-945)/(92160*df^4)) }
             x * (1 + (x2+1 + ((5*x2+16)*x2+3 + (((3*x2+19)*x2+17)*x2-15 +
                                                 ((((79*x2+776)*x2+1482)*x2-1920)*x2-945)/(240*df)) /
                               n4) / (24*df)) / n4) }
           ## = x + (x^3+x)/4df + (*)/96df^2 + (*)/*df^3 + (79x^9+776x^7+1482x^5-1920x^3-945x)/92160df^4
           )
}

