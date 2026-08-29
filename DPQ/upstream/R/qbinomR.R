### Was part of ./qnbinom-slow-ex.R
###               ~~~~~~~~~~~~~~~~~

## MM_FIXME_: solve  "small 'size'"  cases: do_search() is really "too cheap" sometimes, there
## ---

if(FALSE) ## C source:
invisible( "~/R/D/r-devel/R/src/nmath/qbinom.c" )

##======= TODO:  this is practically identical to  doSearch_pnb() in  ./qnbinomR.R === use ONE function
##                                             and doSearch_poi() in  ./qpoisR.R   ...  "  "    "
##  with  (<pDist>, <pDname>)  arguments, here (pnbinom, "pnbinom")
## ==> can MPFR-ify it by passing pnbinomMpfr()..

##' @title For P(x) := pbinom(x, n, pr),  find p-quantile  y = y(p)   :<==>   P(y-1) < p <= P(y)
##'   @param y    current guess
##'   @param z   == pbinom(y, ..)
##'   @param p    target probability
##'   @param n, pr : (n = size); parameters of the  binomial
##'   @param incr increment, integer valued >= 1.
##'   @param lower.tail
##'   @param log.p
##'   @param trace
##' @return  root 'y'  and  z[0] = pbinom(y, ..)
##' @author Martin Maechler
doSearch_bin <- function(y, z, p, n, pr, incr,
                         lower.tail=TRUE, log.p=FALSE, trace = 0)
{
    formatI <- function(x) format(x, scientific = 16)
    if(y < 0) y <- 0
    iStr <- if(incr != 1) sprintf(" incr = %s", formatI(incr)) else ""
    left <- if(lower.tail) (z >= p) else (z < p)
    it <- 0L
    if(left) { ## search to the left
        if(trace) cat(sprintf("doSearch() to the left (pbinom(y=%g, *) = z = %g  %s  p=%g):%s\n",
                              y,z, if(lower.tail) ">=" else "<", p, iStr))
	repeat {
            it <- it + 1L
	    if(y > 0)
                newz <- pbinom(y - incr, n, pr, lower.tail=lower.tail, log.p=log.p)
            else if(y < 0)
                y <- 0
	    if(y == 0 || is.na(newz) || ## pbinom() gives NaN rarely because of pbeta(*, log.p=TRUE) bug
               if(lower.tail) (newz < p) else (newz >= p)) {
                if(trace >= 2) cat(sprintf(
			"  new y=%.15g, pbinom(y-incr,*) %s p; ==> doSearch() returns prev z=%g after %d iter.\n",
                        y, if(lower.tail) "<" else ">=", z, it))
                ## pbinom(y-incr, *) < p <= prev.z = pbinom(y, *)
		return(list(y=y, bin.y = z, iter = it))
            }
	    y <- max(0, y - incr)
            z <- newz
	}
    }
    else { ## !left : search to the right
        if(trace) cat(sprintf("doSearch() to the right (pbinom(y=%g, *) = z = %g  %s  p=%g):%s\n",
                              y,z, if(lower.tail) "<" else ">=", p, iStr))
	repeat {
            it <- it + 1L
	    y <- y + incr
	    z <- pbinom(y, n, pr, lower.tail=lower.tail, log.p=log.p)
	    if(is.na(z) || ## pbinom() returning NaN  {pbeta() bug}
               if(lower.tail) z >= p else z < p) {
                if(trace >= 2) cat(sprintf(
			"  new y=%.15g, z=%g = pbinom(y,*) %s p ==> doSearch() returns after %d iter.\n",
                        y, z, if(lower.tail) ">=" else "<", it))
                ## pbinom(y, *) >= p > prev.z = pbinom(prev.y, *)= pbinom(y-incr, *)  <===>
                ## pbinom(y-incr, *) = prev.z < p <= pbinom(y, *)
		return(list(y=y, bin.y = z, iter = it))
            }
	}
    }
}

qbinomR1 <- function(p, size, prob, lower.tail=TRUE, log.p=FALSE,
                     yLarge = 4096, # was hard wired to 1e5
                     incF = 1/64,   # was hard wired to .001
                     iShrink = 8,   # was hard wired to 100
                     relTol = 1e-15,# was hard wired to 1e-15
                     pfEps.n = 8,   # was hard wired to 64: "fuzz to ensure left continuity"
                     pfEps.L = 2,   # was hard wired to 64:   "   "   ..
                     fpf = 4, # *MUST* be >= 1 (did not exist previously)
                     trace = 0)
{
    stopifnot(length(p) == 1, length(size) == 1, length(prob) == 1)
                                        # #ifdef IEEE_754
    if (is.na(p) || is.na(size) || is.na(prob))
	return(p + size + prob)
                                        # #endif

    if (prob <= 0 || prob > 1 || size < 0) {
        warning("returning NaN because of (prob,size) out of range")
        return(NaN)
    }
    if (prob == 1 || size == 0) return(0)

    ## R_Q_P01_boundaries(p, 0, size); (n == size)
    ## {TODO? simplify via if(p == .D_0(log.p)) .. as in ./t-dist.R }
    if(log.p) {
        if(p == -Inf) return(if(lower.tail) 0 else size)
        if(p ==   0 ) return(if(lower.tail) size else 0)
        if(p > 0) { warning("p > 0 -> returning NaN"); return(NaN) }
        ## now  p \in (-Inf, 0)  (open interval)
    } else {
        if(p == 0) return(if(lower.tail) 0 else size)
        if(p == 1) return(if(lower.tail) size else 0)
        if(p < 0 || p > 1) {
            warning("p out of [0,1] -> returning NaN"); return(NaN) }
        ## now  p \in (0, 1)     (open interval)
    }

    q <- 1 - prob;
    mu <- size * prob;
    sigma <- sqrt(size * prob * q);
    gamma <- (q - prob) / sigma;
    if(trace) {
        cat(sprintf("qbinom(p=%.12g, size=%.15g, prob=%g, l.t.=%d, log=%d):\n
                  mu=%g, sigma=%g, gamma=%g;\n",
                  p, size, prob, lower.tail, log.p, mu, sigma, gamma))
    }

    ## /* Note : "same" code in qpois.c, qbinom.c, qbinom.c --
    ##  * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if(!lower.tail || log.p) {
	p_n <- .DT_qIv(p, lower.tail, log.p) ## need check again (cancellation!):
 	if(trace) cat(sprintf("  upper tail or log_p: (p_n, 1-p_n) = (%.15g, %.15g)\n", p_n, 1-p_n));
	if (p_n == 0) message("p_n=0: NO LONGER return(0)\n")
	if (p_n == 1) message("p_n=1: NO LONGER return(Inf)\n")
    } else
        p_n <- p
    ## temporary hack --- FIXME --- */
    if (p_n + 1.01*.Machine$double.eps >= 1.) {
        if(trace)
            cat("p_n + 1.01 * c_eps >= 1 ; (1-p = ",format(1-p),
                ")   ==> NOW LONGER returning  Inf  (Hack -- FIXME ?)\n", sep="")
        ## return(Inf)
    }

    ## y := approx.value (Cornish-Fisher expansion) :  */
    z <- qnorm(p, lower.tail=lower.tail, log.p=log.p)
    y <- round(mu + sigma * (z + gamma * (z*z - 1) / 6)) # R_forceint(.)
    if(trace) cat(sprintf("Cornish-Fisher: initial z=%g, y=%g; ", z,y))
    if(y > size) # way off
        y <- size
    else if(y < 0) ## (new check; (when) does it happen ?? )
        y <- 0.
    z <- pbinom(y, size, prob, lower.tail=lower.tail, log.p=log.p)
    if(trace) cat(sprintf(" then: y=%g, z= pbinom(y,*) = %g;", y,z))

    ## fuzz to ensure left continuity: */
    ##p <- p*( 1 - pEps.f *.Machine$double.eps ) ## pEps.f was hard wired to 64
    ## fuzz to ensure left continuity: do not loose too much (=> error in upper tail)
    c.eps <- .Machine$double.eps
    if(log.p) { ## <==> p \in [-Inf, 0]  different adjustment: "other sign"
        e <- pfEps.L * c.eps
    	if(lower.tail && p > -.Machine$x.max)# prevent underflow to -Inf
    	    p <- p * (1 + e)
    	else ## if(p < - .Machine$x.min) # not too close to -0
    	    p <- p * (1 - e)
    } else { ##// not log scale
        e <- pfEps.n * c.eps
	if(lower.tail)
	    p <- p * (1 - e)
	else if(1 - p > fpf*e)# otherwise get p > 1
	    p <- p * (1 + e)
    }
    if(trace) cat(sprintf(" left-cont. fuzz => p = %15.11g\n", p))

    ## If the C-F value is not too large a simple search is OK */
    if(y < yLarge) {

        doSearch_bin(y=y, z=z, p=p, n=size, pr=prob, incr = 1,
                     lower.tail=lower.tail, log.p=log.p, trace=trace)$ y

        ## pbinom(y-1, *) < p <= pbinom(y, *)   ==> y is solution

    } else { ## y >= yLarge, be a bit cleverer in the search:
        ## use larger increments, notably initially:
	incr <- floor(y * incF) ## incF was hard-wired to 0.001
        if(trace) cat(sprintf("large y: --> use larger increments than 1: incr=%s\n",
                              format(incr, scientific=16)))
        qIt <- totIt <- 0L
	repeat {
	    oldincr <- incr
	    yz <- doSearch_bin(y=y, z=z, p=p, n=size, pr=prob, incr = incr,
                               lower.tail=lower.tail, log.p=log.p, trace=trace)
            ## We know (for both, left or right search !) :
            ## pbinom(y-oldincr, *) < p <= pbinom(y, *)
            y <- yz$y
            z <- yz$bin.y
            totIt <- totIt + yz$iter
            qIt   <- qIt   + 1L

	    incr <- max(1, floor(incr/iShrink)); ## was (incr/100)
            ## while(oldincr > 1 && incr > y*1e-15);
            if(oldincr <= 1 || incr <= y * relTol)
                break
	}
        if(trace) cat(sprintf("  \\--> %s; needed %d doSearch() calls, total %d iterations\n",
                              if(oldincr == 1) "oldincr=1"
                              else sprintf("oldincr = %s, incr = %s < %.11g = y * relTol",
                                           format(oldincr, scientific = 9),
                                           format(incr,    scientific = 9), y*relTol),
                              qIt, totIt))
        ## return
	y
    }
} ## qbinomR1()

## A version vectorized in 'p' :
qbinomR <- Vectorize(qbinomR1, c("p", "size", "prob"))
                                        # but not {lower.tail, log.p, iShrink, ....}

