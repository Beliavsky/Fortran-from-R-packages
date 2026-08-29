#### qnorm(), pnorm() etc
#### --------------------

### Bounds on  1-Phi(.) = pnorm(*, lower.tail=FALSE) -- typically generalized Mill's ratios

### From Lutz Dümbgen (2010)'s arXiv
invisible("https://arxiv.org/abs/1012.2063")

## Sampford(1953)'s upper bound for 1-Phi(x)  --- Dumbgen's (5), p.2 :
pnormU_S53 <- function(x, lower.tail=FALSE, log.p=FALSE) {
    stopifnot(x >= 0)
    ## upper tail and *not* log.p :
    ##  4*dnorm(x) / (sqrt(8+x^2) + 3*x)

    ## log.p=TRUE and upper tail, i.e.  !lower.tail :
    ## log(4) + dnorm(x, log=TRUE) - log(sqrt(8+x^2) + 3*x)
    ## with less quick overflow for very large 'x' :
    ## sqrt(8+x^2) + 3x = |x| sqrt(1 + 8/x^2) + 3x = x(sqrt(1 + 8/x^2) + 3)
    ##
    ## dnorm(x, log=TRUE) - log(x) + (log(4) - log(sqrt(1 + 8/x^2) + 3))
    r <- dnorm(x, log=TRUE) - log(x) + log(4 / (3 + sqrt(1 + (8/x)/x)))
    if(log.p) {
        if(lower.tail) ## log(1 - exp(r)) = log1mexp(-r)
            log1mexp(-r)
        else ## upper tail: log(1 - (1 - exp(r))) = r
            r
    } else {
        if(lower.tail) -expm1(r) else exp(r)
    }
}

## Duembgen's lower bound (10), p.6
## { which is strictly better than Komatu(1955)'s lower bound (3) }
pnormL_LD10 <- function(x, lower.tail=FALSE, log.p=FALSE) {
    stopifnot(x > 0)
    ## non-log, upper tail :
    ## 1-Phi(x) >  ~=~ pi*dnorm(x) / ((pi-1)*x + sqrt(2*pi + x^2))
    ## log.p=TRUE and upper tail, i.e.  !lower.tail :
    r <- dnorm(x, log=TRUE) - log(x) + log(pi / (pi + sqrt(1 + (2*pi/x)/x) -1))
    if(log.p) {
        if(lower.tail) ## log(1 - exp(r)) = log1mexp(-r)
            log1mexp(-r)
        else
            r
    } else {
        if(lower.tail) -expm1(r) else exp(r)
    }
}

## "Upper tail" asymptotic approximation of Q(x) --
## Using the asymptotic series of  Abramowitz & Stegun, 26.2.13, p.932

pnormAsymp <- function(x, k, lower.tail=FALSE, log.p=FALSE) {
    stopifnot(k == round(k), 0 <= (k <- as.integer(k)), k <= 5)
    ##
    ## FIXME: Want to keep vectorized; TODO: lower.tail as *vector* ==> can *swap* lower.tail *where* x < 0
    stopifnot(x >= 0)
    ## if(x < 0) { ## swap tail
    ##     x <- -x
    ##     lower.tail <- ! lower.tail
    ## }
    r <- dnorm(x, log=TRUE) - log(x) ## <==> exp(r) = phi(x)/x
    if(k > 0) {
        xsq <- x*x
        del <-
            switch(k,
                   1/(xsq + 2.) , # k = 1
                   (1 - 1/(xsq + 4.))/(xsq + 2.) , # k = 2
                   ## k = 3:
                   (1 - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.),
                   ## k = 4:
                   (1. - (1. - (5. - 9/(xsq+8.)) / (xsq+6.)) / (xsq+4.)) / (xsq+2.),
                   ## k = 5:
                   (1. - (1. - (5. - (9 - 129/(xsq+10))/ (xsq+8)) / (xsq+6)) / (xsq+4)) / (xsq+2),
                   stop("invalid 'k': ", k)) # should never happen
        r <- r + log1p(-del)
    }
    if(log.p) {
        if(lower.tail) ## log(1 - exp(r)) = log1mexp(-r)
            log1mexp(-r)
        else ## upper tail: log(1 - (1 - exp(r))) = r
            r
    } else {
        if(lower.tail) -expm1(r) else exp(r)
    }
}


## TODO: if 'k' is NA / NULL ==> find optimal 'k' for given 'x'
## ===   ==> but that needs *separate* function which we need  Vectorize()




###----------- NB>  qnormAppr(), qnormUappr() and qnormUappr6() -->>>> ./beta-fns.R <<<<<<
##                  =========    ==========       ===========            ~~~~~~~~~~


### R version of   ~/R/D/r-devel/R/src/nmath/qnorm.c

## Mathlib : A C Library of Special Functions
## Copyright (C) 2000--2020 The R Core Team
## Copyright (C) 1998       Ross Ihaka
## based on AS 111 (C) 1977 Royal Statistical Society
## and   on AS 241 (C) 1988 Royal Statistical Society

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/

## SYNOPSIS

##     double qnorm5(double p, double mu, double sigma,
##     	      int lower_tail, int log_p)
##           {qnorm (..) is synonymous and preferred inside R}

## DESCRIPTION

##     Compute the quantile function for the normal distribution.

##     For small to moderate probabilities, algorithm referenced
##     below is used to obtain an initial approximation which is
##     polished with a final Newton step.

##     For very large arguments, an algorithm of Wichura is used.

## REFERENCE

##     Beasley, J. D. and S. G. Springer (1977).
##     Algorithm AS 111: The percentage points of the normal distribution,
##     Applied Statistics, 26, 118-121.

##     Wichura, M.J. (1988).
##     Algorithm AS 241: The Percentage Points of the Normal Distribution.
##     Applied Statistics, 37, 477-484.


## <--> ../man/qnormR.Rd

qnormR1 <- function(p, mu=0, sd=1, lower.tail=TRUE, log.p=FALSE,
                    trace = 0,
                    version = c("4.0.x", "1.0.x", "1.0_noN", "2020-10-17", "2022-08-04"))
{
    stopifnot(length(p) == 1, length(mu) == 1, length(sd) == 1,
              length(lower.tail) == 1, length(log.p) == 1,
              !is.na(lower.tail), !is.na(log.p))
    version <- match.arg(version)
    if(is.na(p) || is.na(mu) || is.na(sd)) return(p+mu+sd)
    ## R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF) :
    if(p == .D_0(log.p)) return(if(lower.tail) -Inf else  Inf)
    if(p == .D_1(log.p)) return(if(lower.tail)  Inf else -Inf)
    if(p < .D_0(log.p) ||
       p > .D_1(log.p)) { warning("p out of range"); return(NaN) }

    if(sd < 0) { warning("sd < 0"); return(NaN) }
    if(sd == 0)	return(mu)

    p. <- .DT_qIv(p, lower.tail, log.p=log.p) # = lower_tail prob (in any case)
    q = p. - 0.5;

    if(trace)
        cat(sprintf(
            "qnormR1(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d, version=\"%s\"): q = %g\n",
            p,mu,sd, lower.tail, log.p,  version,  q))

    if(startsWith(version, "1.0")) {
        final_Newton <- (version == "1.0.x")
        if (abs(q) <= 0.42) {
            ## 0.08 <= p <= 0.92 */
            r <- q * q;
            val <- q * (((-25.44106049637 * r + 41.39119773534) * r
                - 18.61500062529) * r + 2.50662823884
            ) / ((((3.13082909833 * r - 21.06224101826) * r
                + 23.08336743743) * r + -8.47351093090) * r + 1.)
        }
        else {
            ## p < 0.08 or p > 0.92, set r = min(p, 1 - p) */

            if (q > 0)
                r = .DT_CIv(p, lower.tail, log.p) # 1-p
            else
                r = p. # = R_DT_Iv(p) ^=  p
            if(trace)
                cat(sprintf("\t 'middle p': r = %7g\n", r))

            if(r > DBL_EPSILON) {
                r = sqrt(- if(log.p &&
                              ((lower.tail && q <= 0) || (!lower.tail && q > 0))) p else log(r))
                if(trace)
                    cat(sprintf("\t new r = %7g ( =? sqrt(- log(r)) )\n", r))

                val = (((2.32121276858 * r + 4.85014127135) * r
                        - 2.29796479134) * r - 2.78718931138
                      ) / ((1.63706781897 * r + 3.54388924762) * r + 1.)
	    if (q < 0)
		val = -val;
            }
            else if(r >= DBL_MIN) { ## r = p <= eps : Use Wichura */
                val <- -2 * if(log.p) .D_Lval(p, lower.tail)
                            else  log(.D_Lval(p, lower.tail))
                r = log(2 * M_PI * val);
                if(trace)
                    cat(sprintf("\t DBL_MIN <= r <= DBL_EPS: val = %g, new r = %g\n",
                                val, r))
                p = val * val;
                r = r/val + (2 - r)/p + (-14 + 6 * r - r * r)/(2 * p * val);
                val = sqrt(val * (1 - r));
                if(q < 0)
                    val = -val;
                final_Newton <- FALSE
            }
            else {
                ## if(trace)
                ##     cat("\t r < DBL_MIN : giving up (-> +- Inf \n")
                cat("r < DBL_MIN (ME_RANGE error) -- infinite result\n") # ML_ERROR(ME_RANGE);
                return(if(q < 0) ML_NEGINF else ML_POSINF)
            }
        }
        if(final_Newton) {
            ## FIXME: This could be improved when log.p or !lower.tail ?
            ## 	  (using p, not p. , and a different derivative )
            if(trace)
                cat(sprintf("\t before final step: val = %7g\n", val))
            ## Final Newton step: */
            val = val - (pnorm(val) - p.) / dnorm(val);
        }
    }

    ## otherwise -- use AS 241 ---

##  double ppnd16_(double *p, long *ifault)*/
##       ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

##       Produces the normal deviate Z corresponding to a given lower
##       tail area of P; Z is accurate to about 1 part in 10**16.

##       (original fortran code used PARAMETER(..) for the coefficients
##        and provided hash codes for checking them...)

    else if (abs(q) <= .425) { ## 0.075 <= p <= 0.925 */
        r  <- .180625 - q * q; ## History: AS241 has SPLIT1 = 0.425, CONST1 = 0.180625
        if(trace) cat(sprintf(" --> usual rational form 1, r=%12.7g\n", r))
	val <-
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608
            ) /
            (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.)
    }
    else { ## closer than 0.075 from {0,1} boundary */
           ## "r" := lp := log(p~);  p~ = min(p, 1-p) < 0.075 :
	lp <-
            if(log.p && ((lower.tail && q <= 0) || (!lower.tail && q > 0))) {
                p
            } else {
                log( if(q > 0) .DT_CIv(p, lower.tail, log.p) ## 1-p
                     else p. ) # p. = R_DT_Iv(p) ^=  p
            }
	###  r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 ) :
        r  <- sqrt(-lp)
        if(trace) cat(sprintf("\t somewhat close to 0 or 1: r := sqrt(-lp) = %11g\n", r))
        if(is.na(r)) # safety-check, should never happen
            warning("r = sqrt( - log(min(p,1-p)) )  is NA -- \"we have lost\"")

        if (!is.na(r) && r <= 5.) { ## <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r <- r + -1.6; ## History: AS241 has SPLIT2 = 5.0, CONST2 = 1.6
            if(trace) cat(sprintf("\t as r <= 5, using rational form 2, for (updated) r=%12.7g\n", r))
            val <- (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                  1.42343711074968357734
            ) / (((((((r *
                       1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                      r + .0151986665636164571966) * r +
                     .14810397642748007459) * r + .68976733498510000455) *
                   r + 1.6763848301838038494) * r +
                  2.05319162663775882187) * r + 1.)
        }
        else { ## p is very close to  0 or 1:  r > 5 <==> min(p,1-p) < exp(-25) = 1.3888..e-11
          if(version == "2020-10-17" && r >= 816) { # p is *extremly* close to 0 or 1 - only possibly when log_p =TRUE
              if(trace) cat("\t *extremely* close to 0 or 1; ==> r large: asymptotic sqrt(2*s) = r*sqrt(2)\n")
              ## Using the asymptotical formula -- is *not* optimal but uniformly better than branch below
              val = r * M_SQRT2;
          }
          else if(version == "2022-08-04" && r > 27) {
              ## p is *really* close to 0 or 1 .. practically only when log_p =TRUE
            if(trace) cat("\t *really* close to 0 or 1; ==> using an asymptotic formula; ")
	    if(r >= 6.4e8) { ## p is *very extremly* close to 0 or 1
		## Using the asymptotical formula ("0-th order"): qn = sqrt(2*s)
                if(trace) cat("  r large:  sqrt(2*s) = r*sqrt(2)\n")
		val = r * M_SQRT2;
	    } else {
		s2 <- -ldexp(lp, 1) ## = -2*lp = 2s
                x2 <- s2 - log(M_2PI * s2); ## = xs_1
                if(trace) cat(sprintf("  1st order x2=%11g \n", x2))
		## if(r >= 36000.)  # use x2 = xs_1 above
		if(r < 36000.) { ## <==> s < 36000^2
		    x2 = s2 - log(M_2PI * x2) - 2./(2. + x2); ## == xs_2
                    if(trace) cat(sprintf("  2nd order x2=%11g \n", x2))
		    if(r < 840.) { ## 27 < r < 840
			x2 = s2 - log(M_2PI * x2) + 2*log1p(- (1 - 1/(4 + x2))/(2. + x2)); ## == xs_3
                        if(trace) cat(sprintf("  3rd order x2=%11g \n", x2))
			if(r < 109.) { ## 27 < r < 109
			  x2 = s2 - log(M_2PI * x2) +
			      2*log1p(- (1 - (1 - 5/(6 + x2))/(4. + x2))/(2. + x2)); ## == xs_4
                          if(trace) cat(sprintf("  4-th order x2=%11g \n", x2))
			  if(r < 55.) { ## 27 < r < 55
                            if(trace) cat("27 < r < 55: using 5-th order x2\n")
			    x2 = s2 - log(M_2PI * x2) +
			      2*log1p(- (1 - (1 - (5 - 9/(8. + x2))/(6. + x2))/(4. + x2))/(2. + x2)); ## == xs_5
			  }
			}
		    }
		}
                val = sqrt(x2);
	    }
          }
          else { ## r > 5 and depending on version  r < **
	    ## Wichura, p.478: minimax rational approx R_3(t) is for 5 <= t <= 27  (t :== r)
            r <- r + -5.;
            if(trace) cat(sprintf("\t r > 5, using rational form R_3(t), for t=%12.7g\n", r))
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772
            ) / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
          }
        }
	if(q < 0)
	    val <- -val;
    }
    ## return
    mu + sd * val
} # qnormR1()

## FIXME:    Rmpfr / mpfrize  as bd0() in  --> ./dgamma.R
## FIXME(2): In that case, need N(.) for sprintf(.. %g ..)  -- or even just format(<num>)
qnormR <- Vectorize(qnormR1, c("p", "mu", "sd"))
##==

### Partly modeled after qnormUappr()  from ./beta-fns.R

qnormAsymp <- function(p, lp = .DT_Clog(p, lower.tail=lower.tail, log.p=log.p),
                           # ~= log(1-p) -- independent of lower.tail, log.p
                       order,
                       lower.tail = TRUE, log.p = missing(p))
{
    stopifnot(length(order) == 1L, order == (ord <- as.integer(order)), 0L <= ord, ord <= 5L)
    if(missing(p)) { ## log.p unused;  lp = log(1-p)  <==>  e^lp = 1-p  <==>  p = 1 - e^lp
        p. <- -expm1(lp)
        ## swap p <--> 1-p -- so we are where approximation is better
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector
    } else {
        p. <- .D_qIv(p, log.p)
        ## swap p <--> 1-p   (as above)
        swap <- if(lower.tail) p. < 1/2 else p. > 1/2 # logical vector
        p[swap] <- if(log.p) log1mexp(-p[swap]) else 1 - p[swap]
    }
    iFin <- is.finite(lp)
    ##  r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 )
    x2 <- s2 <- -ldexp(lp[iFin], 1) ## = -2*lp = 2s =: xs_0
    if(ord >= 1L) {
        x2 <- s2 - log(M_2PI * x2); ## = xs_1
        if(ord >= 2L) { ## need for (r < 36000.)  <==> s < 36000^2
            x2 = s2 - log(M_2PI * x2) - 2./(2. + x2); ## == xs_2
            if(ord >= 3L) { ## need for (r < 840)
                x2 = s2 - log(M_2PI * x2) + 2*log1p(- 1./(2. + x2)*(1 - 1/(4 + x2))); ## == xs_3
                if(ord >= 4L) { ## need for (r < 109)
                    x2 = s2 - log(M_2PI * x2) + 2*log1p(- 1./(2. + x2)*(1 - 1/(4. + x2)*(1 - 5/(6 + x2)))); ## == xs_4
                    if(ord >= 5L) { ## need for (r < 55)
                        x2 = s2 - log(M_2PI * x2) + 2*log1p( - 1./(2. + x2)*
                                                            (1 - 1/(4. + x2)*
                                                             (1 - 1/(6. + x2)*
                                                              (5 - 9/(8. + x2))))); ## == xs_5
                    }
                }
            }
        }
    }
    R <- -lp # incl. Inf
    R[iFin] <- sqrt(x2)
    R[swap] <- -R[swap]
    ## R[p. == 1/2] <- 0
    R
}

## An approximation close to p=1/2, i.e., x ~= 0
## Using Upper Bound P_1(x)  of  Abramowitz & Stegun, 26.2.24, p.933
## and *inverting it
##
qnormU1 <- function(p) {
    neg <- p < 1/2
    r <- sqrt(pi/2 * (-log(4*p*(1-p))))
    r[neg] <- -r[neg]
    r
}

## "CAppr" := [C]enter [Appr]oximation
## TODO: 1) lower.tail=TRUE/FALSE, log.p=FALSE/TRUE
##       2) stable  log(4*p*(1-p))  for log.p=TRUE ==> explored below
## Using Upper Bound P_1(x)  of  Abramowitz & Stegun, 26.2.24, p.933
## and               P_3(x)  of       "         "   , 26.2.25
## *and* iterative plugin
qnormCappr <- function(p, k=1L) {
    stopifnot(length(k <- as.integer(k)) == 1, k >= 1L)
    neg <- p < 1/2
    pih <- pi/2
    x2 <- -pih * (l4p1p <- log(4*p*(1-p)))
    if(k >= 2L) {
        fpi <- 2*(pi-3)/(3*pi^2)
        for(j in 1:(k-1L))
            x2 <- -pih * (l4p1p - fpi*x2^2*exp(-x2/2))
    }
    r <- sqrt(x2)
    r[neg] <- -r[neg]
    r
}

##' Compute log(4*p'*(1-p')) for p' \in [0,1];  where
##' p' = p (log.p=FALSE) -- lower.tail does not matter as we need both p,1-p !
##' p' = exp(p) or exp(1-p) for log.p=TRUE
log4p1p <- function(p, log.p=FALSE) {
    if(log.p) {
        ## log(4) + p + log1mexp(-p)
        ## +/- equivalent, but slightly faster; around p=log(1/2) even slightly less problematic
        log(-4*expm1(p)) + p
    } else { # p' = p or 1-p -- does not matter: p'(1-p') == p(1-p)
        log(4*p*(1-p)) # is the best from those tried (see Rmpfr below)
        ## == log1p(-4*(p-1/2)^2)
        ## log(4)+log(p)+log1p(-p)
    }
}
