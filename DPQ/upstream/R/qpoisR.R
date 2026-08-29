
## MM: _FIXME_solve  "small 'lambda'"  cases: do_search() is really "too cheap" sometimes, there
## ---

if(FALSE) ## C source:
invisible( "~/R/D/r-devel/R/src/nmath/qpois.c" )

##======= TODO(?):  this is almost identical to  doSearch_bin() in  ./qbinomR.R === use ONE function
##  with  (<pDist>, <pDname>)  arguments, here (ppois, "ppois")
## ==> can MPFR-ify it by passing ...
## "Extra problem":  here 1 param  lambda <-->  n, pr (2 param) for {n}binom()

##' @title For P(x) := ppois(x, n, pr),  find p-quantile  y = y(p)   :<==>   P(y-1) < p <= P(y)
##'   @param y    current guess
##'   @param z   == ppois(y, ..)
##'   @param p    target probability
##'   @param lambda  parameter of Poisson
##'   @param incr increment, integer valued >= 1.
##'   @param lower.tail
##'   @param log.p
##'   @param trace
##' @return  root 'y'  and  z[0] = ppois(y, ..)
##' @author Martin Maechler
doSearch_poi <- function(y, z, p, lambda, incr,
                         lower.tail=TRUE, log.p=FALSE, trace = 0)
{
    formatI <- function(x) format(x, scientific = 16)
    if(y < 0) y <- 0
    iStr <- if(incr != 1) sprintf(" incr = %s", formatI(incr)) else ""
    left <- if(lower.tail) (z >= p) else (z < p)
    it <- 0L
    if(left) { ## search to the left
        if(trace) cat(sprintf("doSearch() to the left (ppois(y=%g, *) = z = %g  %s  p=%g):%s\n",
                              y,z, if(lower.tail) ">=" else "<", p, iStr))
	repeat {
            it <- it + 1L
	    if(y > 0)
                newz <- ppois(y - incr, lambda, lower.tail=lower.tail, log.p=log.p)
            else if(y < 0)
                y <- 0
	    if(y == 0 || is.na(newz) || ## ppois() giving NaN -- does this happen??
               if(lower.tail) (newz < p) else (newz >= p)) {
                if(trace >= 2) cat(sprintf(
			"  new y=%.15g, ppois(y-incr,*) %s p; ==> doSearch() returns prev z=%g after %d iter.\n",
                        y, if(lower.tail) "<" else ">=", z, it))
                ## ppois(y-incr, *) < p <= prev.z = ppois(y, *)
		return(list(y=y, poi.y = z, iter = it))
            }
	    y <- max(0, y - incr)
            z <- newz
	}
    }
    else { ## !left : search to the right
        if(trace) cat(sprintf("doSearch() to the right (ppois(y=%g, *) = z = %g  %s  p=%g):%s\n",
                              y,z, if(lower.tail) "<" else ">=", p, iStr))
	repeat {
            it <- it + 1L
	    y <- y + incr
	    z <- ppois(y, lambda, lower.tail=lower.tail, log.p=log.p)
	    if(is.na(z) || ## ppois() giving NaN -- does this happen??
               if(lower.tail) z >= p else z < p) {
                if(trace >= 2) cat(sprintf(
			"  new y=%.15g, z=%g = ppois(y,*) %s p ==> doSearch() returns after %d iter.\n",
                        y, z, if(lower.tail) ">=" else "<", it))
                ## ppois(y, *) >= p > prev.z = ppois(prev.y, *)= ppois(y-incr, *)  <===>
                ## ppois(y-incr, *) = prev.z < p <= ppois(y, *)
		return(list(y=y, poi.y = z, iter = it))
            }
	}
    }
}

qpoisR1 <- function(p, lambda, lower.tail=TRUE, log.p=FALSE,
                    yLarge = 4096, # was hard wired to 1e5
                    incF = 1/64,   # was hard wired to .001
                    iShrink = 8,   # was hard wired to 100
                    relTol = 1e-15,# was hard wired to 1e-15
                    pfEps.n = 8,   # was hard wired to 64: "fuzz to ensure left continuity"
                    pfEps.L = 2,   # was hard wired to 64:   "   "   ..
                    fpf = 4, # *MUST* be >= 1 (did not exist previously)
                    trace = 0)
{
    stopifnot(length(p) == 1, length(lambda) == 1)
    if (is.na(p) || is.na(lambda))
	return(p + lambda)

    if (lambda == 0) return(0)
    if (lambda < 0) {
        warning("returning NaN because of lambda < 0 is out of range")
        return(NaN)
    }

    ## R_Q_P01_boundaries(p, 0, ML_POSINF);
    ## {TODO? simplify via if(p == .D_0(log.p)) .. as in ./t-dist.R }
    if(log.p) {
        if(p == -Inf) return(if(lower.tail) 0 else Inf)
        if(p ==   0 ) return(if(lower.tail) Inf else 0)
        if(p > 0) { warning("p > 0 -> returning NaN"); return(NaN) }
        ## now  p \in (-Inf, 0)  (open interval)
    } else {
        if(p == 0) return(if(lower.tail) 0 else Inf)
        if(p == 1) return(if(lower.tail) Inf else 0)
        if(p < 0 || p > 1) {
            warning("p out of [0,1] -> returning NaN"); return(NaN) }
        ## now  p \in (0, 1)     (open interval)
    }

    mu <- lambda
    sigma <- sqrt(lambda);
    gamma <- 1.0/sigma;

    if(trace) {
        cat(sprintf("qpois(p=%.12g, lambda=%.15g, l.t.=%d, log=%d):\n
                  mu=%g, sigma=%g, gamma=%g;\n",
                  p, lambda, lower.tail, log.p, mu, sigma, gamma))
    }

    ## /* Note : "same" code in qpois.c, qbinom.c, qpois.c --
    ##  * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if(!lower.tail || log.p) {
	p_n <- .DT_qIv(p, lower.tail, log.p) ## need check again (cancellation!):
	if (p_n == 0) message("p_n=0: NO LONGER return(0)\n")
	if (p_n == 1) message("p_n=1: NO LONGER return(Inf)\n")
    } else
        p_n <- p
    ## temporary hack --- FIXME --- */
    if (p_n + 1.01*.Machine$double.eps >= 1.) {
        if(trace)
            cat("p__n + 1.01 * c_eps >= 1 ; (1-p = ",format(1-p),
                ")   ==> NOW LONGER returning  Inf  (Hack -- FIXME ?)\n", sep="")
        ## return(Inf)
    }

    ## y := approx.value (Cornish-Fisher expansion) :  */
    z <- qnorm(p, lower.tail=lower.tail, log.p=log.p)
    y <- round(mu + sigma * (z + gamma * (z*z - 1) / 6)) # R_forceint(.)
    if(trace) cat(sprintf("Cornish-Fisher: initial z=%g, y=%g; ", z,y))
    if(y < 0) y <- 0.; ## e.g., for qpois(0.5, mu = 3, size = 1e-10)

    z <- ppois(y, lambda, lower.tail=lower.tail, log.p=log.p)
    if(trace) cat(sprintf(" then: y=%g, z= ppois(y,*) = %g\n", y,z))

    ## fuzz to ensure left continuity: do not loose too much (=> error in upper tail):
    ## p <- p*( 1 - pEps.f *.Machine$double.eps ) ## pEps.f was hard wired to 64
    c.eps <- .Machine$double.eps
    if(log.p) { ## <==> p \in [-Inf, 0]  different adjustment: "other sign"
        e <- pfEps.L * c.eps
    	if(lower.tail && p > -.Machine$x.max)# prevent underflow to -Inf
    	    p <- p * (1 + e)
    	else if(p < - .Machine$x.min) # not too close to -0
    	    p <- p * (1 - e)
    } else { ##// not log scale
        e <- pfEps.n * c.eps
	if(lower.tail)
	    p <- p * (1 - e)
	else if(1 - p > fpf*e)# otherwise get p > 1
	    p <- p * (1 + e)
    }


    ## If the C-F value is not too large a simple search is OK */
    if(y < yLarge) { ## yLarge was hard-wired to 1e5

        doSearch_poi(y=y, z=z, p=p, lambda=lambda, incr = 1,
                     lower.tail=lower.tail, log.p=log.p, trace=trace)$ y

        ## ppois(y-1, *) < p <= ppois(y, *)   ==> y is solution

    } else { ## y >= yLarge, be a bit cleverer in the search:
        ## use larger increments, notably initially:
	incr <- floor(y * incF) ## incF was hard-wired to 0.001
        if(trace) cat(sprintf("large y: --> use larger increments than 1: incr=%s\n",
                              format(incr, scientific=16)))
        qIt <- totIt <- 0L
	repeat {
	    oldincr <- incr
	    yz <- doSearch_poi(y=y, z=z, p=p, lambda=lambda, incr = incr,
                               lower.tail=lower.tail, log.p=log.p, trace=trace)
            ## We know (for both, left or right search !) :
            ## ppois(y-oldincr, *) < p <= ppois(y, *)
            y <- yz$y
            z <- yz$poi.y
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
} ## qpoisR1()

## A version vectorized in (p, lambda) :
qpoisR <- Vectorize(qpoisR1, c("p", "lambda"))
                                        # but not {lower.tail, log.p, iShrink, ....}

