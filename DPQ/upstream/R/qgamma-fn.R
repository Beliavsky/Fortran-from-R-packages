#### R simulation of C's qgamma()  -*- delete-old-versions: never -*-
#### ---------------------------------
### MM: ~/R/D/r-devel/R/src/nmath/qgamma.c

## using this, see ../tests/qgamma-ex.R , for testing ../tests/chisq-nonc-ex.R
##	                    ~~~~~~~~~~~                        ~~~~~~~~~~~~~~~

### Exports:
### -------
###  qchisqAppr.R (p, nu, g=.., lower.tail,log.p, tol, verbose, kind)
###
###  qchisq.appr.Kind (p, nu, g, lower.tail,log.p, tol, verbose)
###
###  qgamma.R      (p, alpha, scale, lower.tail,log.p, EPS1,EPS2,.......)

## source("/u/maechler/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/beta-fns.R")
## if(FALSE) ## lgamma1p(), ... ----- also loads
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")
## --> .DT_qIv() etc etc


### This is the R equivalent of the .C() calling
### qgammaAppr() in ../qchisqAppr.R
##                  ---------------
qchisqAppr.R <- function(p, df, lower.tail = TRUE, log.p = FALSE,
                         tol = 5e-7, maxit = 1000, verbose = getOption('verbose'),
                         kind = NULL)

{
    ## Purpose: Cheap fast "initial" approximation to qgamma()
    ## --- also to explore the different kinds / cutoffs..
    ## ----------------------------------------------------------------------
    ## Arguments: tol: tolerance with default = EPS2 of qgamma()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 23 Mar 2004, 16:16

    if(length(df) != 1 || length(lower.tail) != 1 || length(log.p) != 1)
        stop("arguments 'df', 'lower.tail', and 'log.p' must have length 1")
    if(log.p) stopifnot(p <= 0) else stopifnot(0 <= p, p <= 1)
    n <- length(p)
    if (df <= 0) return(rep(NaN,n))

    alpha <- 0.5 * df ##/* = [pq]gamma() shape */
    c <- alpha-1
    g <- lgamma(df/2)
    if(!is.null(kind))
        kind <- match.arg(kind,
                          c("chi.small", "WH", "WHchk", "p1WH", "df.small"))

    Cat <- function(...) if(verbose > 0) cat(...)

    ## vectorizing the "rest"  ``in p'' (including choice of 'kind' which depends on p):
    vapply(p, FUN.VALUE = 1, function(p) {

        ## test arguments and initialise
        if (is.na(p) || is.na(df))
            return(p + df)
        p1 <- .DT_log(p, lower.tail, log.p)
        if(is.null(kind)) ## default
            kind <- {
                if(df < -1.24 * p1) "chi.small"
                else if(df > 0.32) {
                    ## 'ch' not known! if(ch > 2.2*df + 6) "p1WH" else "WH"
                    "WHchk"
                }
                else "df.small"
            }
        switch(kind,
               "chi.small" = { ##/* for small chi-squared */

		   ## log(alpha) + g = log(alpha) + log(gamma(alpha)) =
		   ## = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
		   ##  catastrophic cancellation when alpha << 1
		   lgam1pa <- if(alpha < 0.25) lgamma1p(alpha) else log(alpha) + g
		   ch <- exp((lgam1pa + p1)/alpha + M_LN2)
		   Cat(sprintf("(p1,df)=(%g,%g) ==> small chi-sq., ch0 = %g\n",
			       p1,df,ch))
               },
               "WH" = ,
               "WHchk" = ,
               "p1WH" = { ##/*  using Wilson and Hilferty estimate */

                   x <- qnorm(p, 0, 1, lower.tail, log.p)
                   p1 <- 2./(9*df)
                   ch <- df* (x*sqrt(p1) + 1-p1)^3

                   Cat(sprintf(" df=%g: Wilson-Hilferty; x = %.12g\n",df,x))

                   ##/* approximation for p tending to 1: */
                   if(kind == "p1WH" || (kind == "WHchk" && ch > 2.2*df + 6))
                       ch <- -2*(.DT_Clog(p, lower.tail,log.p)
                                 - c*log(0.5*ch)+ g)
               },
               "df.small" = { ##/* small  df : 1.24*(-log(p)) <= df <= 0.32 */

                   C7 <- 4.67
                   C8 <- 6.66
                   C9 <- 6.73
                   C10 <- 13.32

                   ch <- 0.4
                   a <- .DT_Clog(p, lower.tail, log.p) + g + c*M_LN2

                   Cat(sprintf("'df.small', df=%g, a(p,df) = %19.13g ", df,a))

                   it <- 0; converged <- FALSE
                   while(!converged && it < maxit) {
                       ## q <- ch
                       p1 <- 1. / (1+ch*(C7+ch))
                       p2 <- ch*(C9+ch*(C8+ch))
                       t <- -0.5 +(C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2
                       Del <-  - (1- exp(a+0.5*ch)*p2*p1)/t
                       if(is.na(Del))
                           break
                       ch <- ch + Del
                       it <- it + 1
                       converged <- (abs(Del) <= tol * abs(ch))
                   }
                   Cat(sprintf("%5.0f iterations --> %s\n",it,
                               paste(if(!converged)"NOT", "converged")))

               }) ## end{ switch( kind ) }
        ch
    })# vapply(p, *)
} ## end qchisqAppr.R()

## qchisq.appr.Kind() is currently in ../qchisqAppr.R
## ==================                   ~~~~~~~~~~~~~

qgamma.R <- function(p, alpha, scale = 1, lower.tail = TRUE, log.p = FALSE,
                     EPS1 = 1e-2,
                     EPS2 = 5e-7,##/* final precision of AS 91 */
                     epsN= 1e-15, ##/* precision of Newton step / iterations */
                     maxit = 1000,
                     pMin = 1e-100,
                     pMax = (1-1e-14),
                     verbose = getOption('verbose')
                     )

{
    ## test arguments

    if(length(p) != 1 || length(alpha) != 1 || length(scale) != 1 ||
       length(lower.tail) != 1 || length(log.p) != 1)
        stop("all arguments but 'p'  must have length 1")

    if (is.na(p) || is.na(alpha) || is.na(scale))
	return(p + alpha + scale)
    ## R.Q.P01.check(p):
    if(log.p) stopifnot(p <= 0) else stopifnot(0 <= p && p <= 1)
    if (alpha < 0) stop("alpha < 0")
    if (scale <= 0) stop("scale <= 0")

    if (alpha == 0) return(0)

    ## initialize

    i420  <- 1./ 420
    i2520 <- 1./ 2520
    i5040 <- 1./ 5040

    max.it.Newton <- 1

    Cat <- function(...) if(verbose > 0) cat(...)
    p. <- .DT_qIv(p, lower.tail, log.p)

    Cat(sprintf("qgamma(p=%7g, alpha=%7g, scale=%7g, l.t.=%2d, log.p=%2d): ",
                p,alpha,scale, lower.tail, log.p))

    g <- lgamma(alpha)##/* log Gamma(v/2) */

    ##/*----- Phase I : Starting Approximation */
    ch <- qchisqAppr.R(p, df = 2*alpha, ## g = g,##/* = lgamma(nu/2) */
                       lower.tail, log.p, tol= EPS1, verbose=verbose)
    do.phaseII <- TRUE
    if(!is.finite(ch)) {
	##/* forget about all iterations! */
	max.it.Newton <- 0; do.phaseII <- FALSE
        warning("ch =", formatC(ch)," is not finite -- bug in qgamma() ?")
    }
    else if(ch < EPS2) {##/* Corrected according to AS 91; MM, May 25, 1999 */
        max.it.Newton <- 20;
	do.phaseII <- FALSE##/* and do Newton steps */
    }
    else if(p. > pMax || p. < pMin) {
        ##/* FIXME: This (cutoff to {0, +Inf}) is far from optimal
        ## * -----  when log.p or !lower.tail : */

	##/* did return ML.POSINF or 0.;	much better: */
	max.it.Newton <- 20;
	do.phaseII <- FALSE##/* and do Newton steps */
    }

    if(do.phaseII) {

        Cat(sprintf("\n==> ch = %10g:", ch))

        ##/*----- Phase II: Iteration
        ## *	Call pgamma() [AS 239]	and calculate seven term taylor series
        ## */
        c <- alpha-1;
        s6 <- (120+c*(346+127*c)) * i5040; ##/* used below, is "const" */

        do.break <- FALSE
        ch0 <- ch ##/* save initial approx. */
        for(i in 1:maxit) {
            q <- ch
            p1 <- 0.5*ch
            p2 <- p. - pgamma(p1, alpha, 1, lower.tail=TRUE, log.p=FALSE)

            if(i == 1) Cat(sprintf(" Ph.II iter; ch=%g, p2=%g\n", ch, p2))
            if(i >= 2) Cat(sprintf("     it=%d,  ch=%g, p2=%g\n", i, ch, p2))

            if(!is.finite(p2)) {
                Cat("--> non finite p2\n")
                ch <- ch0; max.it.Newton <- 27
                do.break <- TRUE ; break
            }

            t <- p2*exp(alpha*M_LN2 + g + p1 - c*log(ch));
            b <- t/ch;
            a <- 0.5*t - b*c;
            s1 <- (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
            s2 <- (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
            s3 <- (210+ a*(462+a*(707+932*a))) * i2520;
            s4 <- (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
            s5 <- (84+2264*a + c*(1175+606*a)) * i2520;

            Del <- t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))))

            if(is.na(Del) || abs(Del) < EPS2*(ch <- ch + Del)) {
                if(is.na(Del)) {
                    Cat("--> Del became NA\n")
                    ch <- ch0; max.it.Newton <- 27
                }
                do.break <- TRUE ; break
            }
        }
        ##/* no convergence in maxit iterations -- but we add Newton now... */
                                        #ifdef DEBUG.q
        if(!do.break)
        warning(sprintf("qgamma(%g) not converged in %d iterations; rel.ch=%g\n",
                        p, maxit, ch/abs(q - ch)))
                                        #endif
        ##/* was
        ## *    ML.ERROR(ME.PRECISION);
        ## * does nothing in R !*/
    }

##/* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
##  --------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision
##
##  * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
##  *
##  * Improved (MM): - only if rel.Err > EPS.N (= 1e-15);
##  *		    - also for lower.tail = FALSE	 or log.p = TRUE
##  * 		    - optionally *iterate* Newton
##  */
    x <- 0.5*scale*ch

    if(max.it.Newton) {
	if (!log.p) {
	    p <- log(p)
	    log.p <- TRUE
	}
	if(x == 0) {
	    x <- .Machine$double.xmin
	    p. <- pgamma(x, alpha,
			 scale=scale, lower.tail=lower.tail, log.p=log.p)
	    Cat(sprintf(" x == 0: p. := pgamma(XMIN,*) = %g; p. - p = %g\n",
			p., p. - p))
	    if(( lower.tail && p. > p *(1 + 1e-7)) ||
               (!lower.tail && p. < p *(1 - 1e-7)))
		return(0.)
	    ## else: continue, using x = DBL_MIN instead of  0
	}
	else
	    p. <- pgamma(x, alpha,
			 scale=scale, lower.tail=lower.tail, log.p=log.p)
    }
    for(i in seq_len(max.it.Newton)) {
	p1 <- p. - p
	if(i == 1)
            Cat(sprintf("\n it=%d: p=%g, x = %g, p.=%g; p1:=D{p}=%g\n",
                        i, p, x, p., p1))
	if(i >= 2) Cat(sprintf("         it=%d,  d{p}=%g\n",    i, p1))

	if(abs(p1) < abs(epsN * p) ||
	   (g <- dgamma(x, alpha, scale=scale, log = log.p)) == .D_0(log.p)) {
            if(i == 1 && g == .D_0(log.p))
                warning("no Newton step done because dgamma(*) = 0 !")
	    break
	}
	else {
	    ##/* delta x = f(x)/f'(x);
            ## * if(log.p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
            ## * ==> f(x)/f'(x) = f*P / P' = f*exp(p.) / P' (since p. = log P(x))
            ## */
	    t <- ifelse(log.p, p1*exp(p. - g), p1/g)##/* = "delta x" */
	    t <- ifelse(lower.tail, x - t, x + t)
            p. <- pgamma(t, alpha,
                         scale=scale, lower.tail=lower.tail, log.p=log.p)
            Cat(sprintf("new t= %15.9g,  p.=%15.9g ", t, p.))
	    if (abs(p. - p) >= abs(p1)) { ##/* no improvement */
                if(i == 1)
                warning("no Newton step done since delta{p} >= last delta")
		break
            }
	    else x <- t
	}
    }
    return (x)
}## qgamma.R()
