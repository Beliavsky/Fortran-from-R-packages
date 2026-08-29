#### R version of  dbinom_raw()          from ~/R/D/r-devel/R/src/nmath/dbinom.c  <<< now (Dec 30, 2023) improved!
####        and dnbinom() / dnbinom_mu()  "   ~/R/D/r-devel/R/src/nmath/dnbinom.c
##
##  (orig. was ~/R/MM/NUMERICS/dpq-functions/dbinom_raw.R )

##' Compute  (1+x)^y  accurately also for |x| << 1

pow1p.1 <- function(x, y,  # C code in ...../dbinom.c
                    ## Use naive/direct (1+x)^y  in two cases: (1) when 1+x is exact
                    ## and (2) when |x| > 1/2 and we have no better algorithm.
                    pow = ((x + 1) - 1) == x || abs(x) > 0.5 || is.na(x))
{   ## comment apart from when developing:
    stopifnot(length(x) == 1, length(y) == 1)
    if(is.na(y))
	if(x == 0) 1. else y # (0+1)^NaN := 1  by standards
    else if(0 <= y && y == trunc(y) && y <= 4.) {
        if(!y) 1
        else switch(as.integer(y),
		    x + 1,  # ^1
		    x*(x + 2) + 1,			# (x+1)^2
		    x*(x*(x + 3.) + 3.) + 1,		# ( . )^3
		    x*(x*(x*(x + 4.) + 6.) + 4.) + 1	# ( . )^4
		    )
    } else {
        if(pow)
            (1 + x)^y
        else # not perfect, e.g., for small |x|, non-huge y, use binom expansion 1 + y*x + y(y-1)/2 x^2 + ..
            exp(y * log1p(x))
    }
}
pow1p  <- Vectorize(pow1p.1, c("x", "y"))

dbinom_raw <- function(x, n, p, q=1-p, log = FALSE # >> ../man/dbinom_raw.Rd <<<
                       , version = c("2008", "R4.4")
                       , verbose = getOption("verbose")
                       )
{
  ## Purpose: R version of dbinom_raw()   { from .../R/src/nmath/dbinom.c }
  ## ----------------------------------------------------------------------
  ## Arguments:  p + q == 1

    ##__________ for Rmpfr mpfr() numbers, need a more accurate stirlerr()
    ##--- otherwise cannot get more than double prec.! >>> ./dgamma.R  (using 'DPQmpfr::stirlerrM()')

    stopifnot(is.logical(log), length(log) == 1)
    ## Recycle to common length
    M <- max(length(x), length(n), length(p), length(q))
    r <- double(M) # result
    if(M == 0) return(r)
    ## M >= 1 :
    x <- rep_len(x,  M)
    n <- rep_len(n,  M)
    p <- rep_len(p,  M)
    q <- rep_len(q,  M)

    version <- match.arg(version)
    if(any(B <- p == 0)) r[B] <- ifelse(x[B] ==  0  , .D_1(log), .D_0(log))
    if(any(B <- q == 0)) r[B] <- ifelse(x[B] == n[B], .D_1(log), .D_0(log))
    BB <- !B & (p != 0)
    if(verbose) cat(sprintf("dbinom_raw(*, version='%s'): #(bndr, rglr): (%d,%d)\n",
                            version, sum(!BB), sum(BB)))
    verb1 <- pmax(0L, verbose - 1L)

    if(any(B <- BB & x == 0)) {
        ii <- which(B)
        if(verbose) cat(sprintf("x=0 for i = %s\n", deparse(ii, control="S")))
	if(any(i0 <- n[ii] == 0)) {
            r[ii[i0]] <- .D_1(log)
            ii <- ii[!i0] # for the rest of this clause
        }
        n. <- n[ii]
        p. <- p[ii]
        q. <- q[ii]
        r[ii] <-
         switch(version,
                "2008" = {
                    lc <- ifelse(p. < 0.1,
                                 -bd0(n., n.*q., verbose=verb1) - n.*p.,
                                 n.*log(q.))
                    .D_exp(lc, log)
                },
                "R4.4" = {
                    ifelse(p. > q.,
                           if(log) n. * log(q.)    else q. ^ n.,
                           ## else   0 < p. <= 1/2
                           if(log) n. * log1p(-p.) else pow1p(-p., n.))
                },
                stop("invalid 'version':", version))
        BB <- BB & !B
    }
    if(any(B <- BB & x == n)) {
        ii <- which(B)
        if(verbose) cat(sprintf("x=n for i = %s\n", deparse(ii, control="S")))
        n. <- n[ii]
        p. <- p[ii]
        q. <- q[ii]
        r[ii] <-
         switch(version,
                "2008" = {
                   lc <- ifelse(q. < 0.1,
                                -bd0(n.,n.*p., verbose=verb1) - n.*q.,
                                n.*log(p.))

                   .D_exp(lc, log)
               },
               "R4.4" = {
                    ifelse(p. > q.,
                           if(log) n. * log1p(-q.) else pow1p(-q., n.),
                           ## else   0 < p. <= 1/2
                           if(log) n. * log  (p.)  else p. ^ n.)
               },
               stop("invalid 'version':", version))

        BB <- BB & !B
    }
    if(any(B <- BB & x < 0 | x > n)) {
        if(verbose) cat(sprintf("have %d x values outside [0,n]\n", sum(B)))
        r[B] <- .D_0(log)
        BB <- BB & !B
    }

    if(any(BB)) {
        ii <- which(BB)
        n <- n[ii]
        x <- x[ii]
        if(verbose) cat(sprintf("main, i = %s\n", deparse(ii, control="S")))
        ##  n*p or n*q can underflow to zero if n and p or q are small.  This
        ##  used to occur in dbeta, and gives NaN as from R 2.3.0.
        lc <- { stirlerr(n, verbose=verb1) - stirlerr(x, verbose=verb1) - stirlerr(n-x, verbose=verb1) -
                    bd0( x , n*p[ii], verbose=verb1) - bd0(n-x, n*q[ii], verbose=verb1)
        }

        ## f = (M.2PI*x*(n-x))/n; could overflow or underflow */
        ##lf <- log(2*pi) + log(x) + log(n-x) - log(n)
        ##                          ---------------- = log((n-x)/n)=log(1 - x/n)
        lf  <- log(2*pi) + log(x) + log1p(-x/n)

        if(verbose) cat(sprintf("  lc=%g, lf=%g ==> lc - 0.5*lf = %g\n", lc,lf,lc - 0.5*lf))

        r[BB] <- .D_exp(lc - 0.5*lf, log)
    }
    r
}

dnbinomR <- function (x, size, prob, log = FALSE, eps = 1e-10)
{
  ## Purpose: R version of R'C level dnbinom() in .../R/src/nmath/dnbinom.c

    x <- floor(x + 1e-7)
    stopifnot(is.logical(log), length(log) == 1,
              0 < prob, prob <= 1, size >= 0, x >= 0)
    M <- max(length(x), length(size), length(prob))
    r <- double(M)
    if(M == 0) return(r)
    x    <- rep_len(x,     M)
    size <- rep_len(size,  M)
    prob <- rep_len(prob,  M)

    ## This is  ** in addition ** to the C code .. and is part of x < eps * size below
    if(any(i0 <- x == 0)) {
        if(any(i0s0 <- i0 & size == 0)) {
            r[i0s0] <- .D_1(log)
        }
        if(any(i0P <- i0 & size > 0)) {  ## x = 0, size > 0
            ## pr(x,...) = pr^n :
            r[i0P] <- if(log) size[i0P]*log(prob[i0P]) else prob[i0P] ^ size[i0P]
        }
           x <-    x[!i0]
        size <- size[!i0]
        prob <- prob[!i0]
    }
    if(any(B <- x < eps * size)) { ## don't use dbinom_raw() but MM's formula
        i <- which(B)
        x. <-    x[i]
        n. <- size[i]
        pr <- prob[i]
        r[!i0][i] <- .D_exp(n. * log(pr) + x. * (log(n.) + log1p(-pr))
                            -lgamma1p(x.) + log1p(x.*(x.-1)/(2*n.)),
                            log)
    }
    if(any(!B)) {
        i <- which(!B)
        x. <- x   [i]
        n. <- size[i]
        pr <- prob[i]
        ans <- dbinom_raw(x = n., n = x.+n., p = pr, q = 1-pr, log = log)
        ## p <- n./(n.+x) ## == 1 if |x| << n.;
        ## better in log case when x < n: log(n/(n+x)) = log(1 - x/(n+x))
        r[!i0][i] <- if(log) ifelse(x. < n., log1p(-x./(n.+x.)), log(n./(n.+x.))) + ans
                     else  n./(n.+x.) * ans
    }
    r
}

dnbinom.mu <- function(x, size, mu, log = FALSE, eps = 1e-10)
{
  ## Purpose: R version of dnbinom_mu() { in .../R/src/nmath/dnbinom.c }
  ## ----------------------------------------------------------------------
  ## Arguments: as  dbinom_mu()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Jul 2008, 21:54

    x <- floor(x + 1e-7)
    stopifnot(is.logical(log), length(log) == 1,
              mu >= 0, size >= 0, x >= 0)
    M <- max(length(x), length(size), length(mu))
    r <- double(M)
    if(M == 0) return(r)
    x    <- rep_len(x,    M)
    size <- rep_len(size, M)
    mu   <- rep_len(mu,   M)

    if(any(i0 <- x == 0)) {
        i <- which(i0)
        x.  <-   x[i]
        mu. <-  mu[i]
        n. <- size[i]
        ## pr(x,...) = p^n = (n / (n+mu))^n  -- carefully evaluated:
	r[i] <- .D_exp(n. * ifelse(n. < mu.,
                                    log(n./(n.+mu.)),
                                    log1p( - mu./(n.+mu.))),
			 log)
	size <- size[!i0]
	   x <-    x[!i0]
	  mu <-   mu[!i0]
    }
    if(length(i <- which(B <- x < eps * size))) { ## don't use dbinom_raw() but MM's formula
        ## log p__r =
        x.  <-   x[i]
        mu. <-  mu[i]
        n. <- size[i]
        r[!i0][i] <- .D_exp(x. * log(n.*mu./(n.+mu.)) - mu.
                            -lgamma1p(x.) + log1p(x.*(x.-1)/(2*n.)),
                            log)
    }
    if(length(i <- which(!B))) {
        x    <-    x[i]
        mu   <-   mu[i]
        size <- size[i]
        ans <- dbinom_raw(x= size, n= x+size,
                          p= size/(size+mu), q= mu/(size+mu),
                          log = log)
        ## p <- size/(size+x) ## == 1 if  |x| << size; better in log case,
        ## log(n/(n+x)) = log(1 - x/(n+x))  if (x < size):
        r[!i0][i] <- if(log) ifelse(x < size, log1p(-x/(size+x)), log(size/(size+x))) + ans
                     else  size/(size+x) * ans
    }
    r
}
