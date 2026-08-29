### A pure-R version of pbeta(), << of R <= 2.2.x -- before we used bratio() ! --
##
##  /usr/local.nfs/app/R/R-2.2.1/src/nmath/pbeta.c
##			 ------- ~~~~~~~~~~~~~~~~~
## For consistency within R, we are using
## R's argument names    (q, shape1, shape2)
## instead of C code's   (x, pin,    qin   )
##
pbetaRv1 <- function(q, shape1, shape2, lower.tail = TRUE,
		     eps = .5 * .Machine$double.eps, ## = 1.11e-16
		     sml = .Machine$double.xmin,     ## = 2.22e-308
		     verbose = 0)
{
  ## Purpose: emulate the C function pbeta_raw() in R -- for diagnosing..
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  3 Apr 2004, 16:28

    if(length(q) != 1 || length(shape1) != 1 || length(shape2) != 1)
	stop("arguments must have length 1 !")
    x <- q; pin <- shape1; qin <- shape2 # back to "C-code compatible" argument names
    isN <- is.numeric(x) && is.numeric(pin) && is.numeric(qin)
    isMpfr <- !isN && any_mpfr(x, pin, qin)
    if(isMpfr) isMpfr <- requireNamespace("Rmpfr")
    ## needed for printing mpfr numbers {-> pkg Rmpfr}, e.g.
    .N <- if(isMpfr) Rmpfr::asNumeric else as.numeric
    if(isMpfr) {
        lbeta <- Rmpfr::lbeta
        ## beta <- Rmpfr::beta
        if(verbose) format <- Rmpfr::format
    }
    Cat <- function(...) if(verbose > 0) cat(..., sep="")

    lneps <- log(eps) #  -36.7368  ---	 but may be lower in mpfr case
    lnsml <- log(sml) # -708.3964  ---	 but may be lower in mpfr case

    ## /* swap tails if x is greater than the mean */
    swap.tail <- (pin / (pin + qin) < x)
    if(swap.tail) {
	y <- 1 - x
	p <- qin
	q <- pin
    }
    else {
	y <- x
	p <- pin
	q <- qin
    }

    if ((p + q) * y / (p + 1) < eps) {

	## /* tail approximation */
	Cat("pbetaR(): _tail approximation_ swap.tail : ", swap.tail,"\n")

	xb <- p * log(max(y, sml)) - log(p) - lbeta(p, q)
	if (xb > lnsml && y != 0) {
	    ans	 <- if(swap.tail == lower.tail) -expm1(xb) else exp(xb)
	} else {
	    ans	 <- if(swap.tail == lower.tail) 1. else 0
	}
    }
    else {
	## /* MM: __ FIXME __ : This takes forever (or ends wrongly)
	##	  when (one or) both p & q  are huge
	##
	##  ./pf.c  now has a cheap fix -- can use that here, but better
	##  "get it right"  (PD to R-core on 20 Feb 2000)

	## /* evaluate the infinite sum first.	term will equal */
	## /* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */

	Ly <- if(swap.tail) log1p(-x) else log(y)
	ps <- q - floor(q)
	xb <- p * Ly
	if (ps == 0) {
	    ps <- 1 # in this case: lbeta(ps,p)= log Beta(1,p) = log(1/p)= -log(p)
	}
	else xb <- xb - (lbeta(ps, p) + log(p))

	Cat("pbetaR() 'normal', swap.tail= ", swap.tail,
	    "; xb=", format(xb), ", lnsml=", format(lnsml), ", ")
	if (xb >= lnsml) {
	    ans	 <- exp(xb)
	    term <- ans * p
            Cat("xb >= lnsml: ans_1 := exp(xb) =", format(ans), ", term = ", format(term), " ")
            if(term == 0)
                Cat("term = 0 => no summation ")
            else if (ps != 1) {
		n <- floor(max(lneps/Ly, 4.0))
		## Now if 'n' is really too large, we should give up!
		if(n > 1e9)
		    stop("'n' := max(4, lneps/Ly) =",format(n),
			 " is way too large")
		n <- as.integer(n)# if it was "mpfr"
		Cat(", n=", n,", ")
		for(i in 1:n) {
		    xi	<- i
		    term <- term * ((xi - ps) * y / xi)
		    ans <- ans + term / (p + xi)
		}
	    }
	    Cat("1st ans= ", format(ans), "\n")
	} else {
            ans <- 0
	    Cat("xb < lnsml ==> no first sum; ans := 0\n")
        }

	## /* now evaluate the finite sum, maybe. */

	if (q > 1) {

	    if(swap.tail) {
		c <- 1./x;##/* == 1/(1 - y) */
		liy <- log(x)
	    }
	    else {
		c <- 1./(1. - y)
		liy <- log1p(-y)
	    }

	    xb	<- p * Ly + q * liy - lbeta(p, q) - log(q)
	    ib	<- floor(max(xb / lnsml, 0.))
	    term <- exp(xb - ib * lnsml)
	    p1	<- q * c / (p + q - 1)

	    Cat(" q > 1: xb=",format(xb),", ib=",.N(ib),", term=",format(term),
		", p1=",format(p1), if(p1 > 1)" > 1 !!!!" else '', "\n")
	    finsum <- 0
	    n <- floor(q)
	    if (q == n)
		n <- n-1
	    Cat("2nd n := n(q) = ", n, ":\n")
	    ##L	 seq() can fail for really large 'n' !!
	    ##L for(i in seq(length=n)) {
	    i <- 1
	    while(i <= n) {
		if (p1 <= 1 && term / eps <= finsum)
		    break
		xi <- i
		term <- (q - xi + 1) * c * term / (p + q - xi)
		## do not: when n = 2e9, almost kills emacs...
		if(verbose) 
		   cat(if(verbose >= 2 && i < 100) sprintf("term=%g ",term) else ".", if(!(i %% 100))"\n")
		if (term > 1) {
		    ib <- ib-1
		    term <- term * sml
		}
		if (ib == 0)
		    finsum <- finsum + term
		##L :
		i <- i+1
	    }
            if(verbose && i <= n) Cat("finite sum converged early: i=",i,"\n")
	    ans <- ans + finsum
	}
	if (swap.tail == lower.tail)
	    ans <- 1 - ans
	ans <- max(min(ans, 1.), 0.)
    }

    ans
} ## /* pbeta_raw() */
