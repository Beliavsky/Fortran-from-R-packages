#### Hypergeometric Distribution
#### ===========================
#### R ( / S) functions :
####	dhyper(x, m, n, k)
####	phyper(q, m, n, k)
###
### Parametrization / Nomenclature :
###
### JohKK  R/S        C    AS152
### ~~~~~  ~~~       ~~
###  Np   = m     =  NR   = nn   : # { WHITE (Red) balls   in the urn }
###  N-Np = n     =  NB   = mm-nn: # { BLACK       balls   in the urn }
###   n   = k     =   n   = kk   : # { balls drawn from the urn }
###         x     =   x   = ll   : # { WHITE (Red) balls AMONG the `n' drawn}
###  ----  ---      -----
###   N   =m+n    = NR+NB = mm   : the TOTAL # { balls   in the urn }
###   p   =m/(m+n)= NR/N  = nn/mm: probability of WHITE(Red)

### where
###  JohKK := Johnson, Kotz & Kemp (1992)
###           Univariate Discrete Distributions, 2nd Ed.
###           [Chapter 6 -- Hypergeometric Distributions]
###  R/S   := the [dpqr]hyper() functions defined in R {and already S}

### Some tests (no longer problematic):
###  ./dpq-functions/Knuesel-splus-ex.R  [ -->  phyper() ]

#### Used to be part of /u/maechler/R/MM/NUMERICS/hyper-dist.R -- till Jan.24 2020

#### First version committed is the code "as in Apr 1999":
####  file | 11782 | Apr 22 1999 | hyper-dist.R

phyperApprAS152 <- function(q, m, n, k)
{
  ## Purpose: Normal Approximation to cumulative Hyperbolic
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 Apr 99, 23:39
  kk <- n
  nn <- m
  mm <- m+n
  ll <- q
  mean <- kk * nn / mm
  sig <- sqrt(mean * (mm - nn) / mm * (mm - kk)/(mm - 1) )
  pnorm(ll + 1/2, mean=mean, sd = sig)
}

phyperIbeta <- function(q, m, n, k)
{
  ## Purpose: Pearson [Incompl.Beta] Approximation to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.90), p.260 -- Bol'shev (1964)
  ## ----------------------------------------------------------------------

  Np <- m; N <- n + m; n <- k; x <- q
  p <- Np/N ; np <- n*p
  xi <- (n + Np -1 - 2*np) / (N-2)

  d.c <- (N-n)*(1-p) + np - 1
  cc <- n*(n-1)*p*(Np-1) / ((N-1)*d.c)
  lam <- (N-2)^2 *np *(N-n)*(1-p) / ((N-1)*d.c*(n + Np-1 - 2*np))
  pbeta(1 - xi,
        lam - x + cc,
        x - cc + 1)
}

phyper1molenaar <- function(q, m, n, k)
{
  ## Purpose: Normal Approximation to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.91), p.261
  ## ----------------------------------------------------------------------

  Np <- m; N <- n + m; n <- k; x <- q
  pnorm(2/sqrt(N-1) *
        (sqrt((x + 1)*(N - Np - n + x + 1)) -
         sqrt((n - x)*(Np - x)) ))
}

phyper2molenaar <- function(q, m, n, k)
{
  ## Purpose: Normal Approximation to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.92), p.261
  ## ----------------------------------------------------------------------
  Np <- m; N <- n + m; n <- k; x <- q
  pnorm(2/sqrt(N) *
        (sqrt((x + .75)*(N - Np - n + x + .75)) -
         sqrt((n - x -.25)*(Np - x - .25)) ))
}

phyperPeizer <- function(q, m, n, k)
{
  ## Purpose: Peizer's extremely good Normal Approx. to cumulative Hyperbolic
  ##  Johnson, Kotz & Kemp (1992):  (6.93) & (6.94), p.261 __CORRECTED__
  ## ----------------------------------------------------------------------
  Np <- m; N <- n + m; n <- k; x <- q
  ## (6.94) -- in proper order!
  nn <- Np			;  n. <- Np     + 1/6
  mm <- N - Np                  ;  m. <- N - Np + 1/6
  r <- n                        ;  r. <- n      + 1/6
  s <- N - n                    ;  s. <- N - n  + 1/6
                                   N. <- N      - 1/6
  A <- x + 1/2                  ;  A. <- x      + 2/3
  B <- Np - x - 1/2             ;  B. <- Np - x - 1/3
  C <- n  - x - 1/2             ;  C. <- n  - x - 1/3
  D <- N - Np - n + x + 1/2     ;  D. <- N - Np - n + x + 2/3

  n <- nn
  m <- mm
  ## After (6.93):
  L <- # NB:  for all 4 log(r_j) :  have r_j ~=~ 1 -- would probably *benefit* from log1p(.) -- cancellation?
    A * log((A*N)/(n*r)) +
    B * log((B*N)/(n*s)) +
    C * log((C*N)/(m*r)) +
    D * log((D*N)/(m*s))
  ## (6.93) :
  pnorm((A.*D. - B.*C.) / abs(A*D - B*C) *
        sqrt(2*L* (m* n* r* s* N.)/
                  (m.*n.*r.*s.*N )))
  # The book wrongly has an extra "2*" before `m* ' (after "2*L* (" ) above
}

### Binomial Approximation(s) ==========
### ======== ~~~~~~~~~~~~~~~~ ==========

##  JohKK  R     C
##  ~~~~~  ~~    ~~
##  Np   = m  =  NR  : # { WHITE (Red) balls   in the urn }
##  N-Np = n  =  NB  : # { BLACK       balls   in the urn }
##   n   = k  =   n  : # { balls drawn from the urn }
##         x  =   x  : # { WHITE (Red) balls AMONG the `n' drawn}

##' Molenaar(1970a) -- formula (6.80) in JohKK -- improved p for Bin(n,p) approximation:
hyper2binomP <- function(x, m,n,k) {
    N <- m+n
    p <- m / N # "Np / N"
    N.n <- N - (k-1)/2
    ## p^{+} =
    (m - x/2)/N.n - k*(x - k*p - 1/2) / (6 * N.n^2)
}

### NOTA BENE/TODO: (this paragraph has not __yet__ been applied below at all !!)
### ---------
### JohKK, p.258 (top) mention the *four* different binomial approximations
### for a given hypergeometric, and then
### "Brunk et al. (1968) .. support the opinion of ...(1961) that it is
### best to use the binomial with smallest power parameter, that is,
### n' = min(n, Np, N-Np, N-n)
### --------------------------
###   ( x ; (  n , p= Np/N) )  <--> ( x  ; (Np , p= n /N)) <-->
###   (n-x; (N-Np, p= n /N) )  <--> (Np-x; (N-n, p= Np/N))
### translated to R's notation:
###   ( q ; ( k ,  p= m/(m+n)) <--> ( q  ; (m,     p= k/(m+n))) <-->
###   (k-q; ( n ,  p= k/(m+n)) <--> (m-q ; (m+n-k, p= m/(m+n)))

##' The support of the hypergeometric distrib. as a function of its parameters
.suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)

##' The two symmetries <---> the four different ways to compute :
phypers <- function(m,n,k, q = .suppHyper(m,n,k), tol = sqrt(.Machine$double.eps)) {
    N <- m+n
    pm <- cbind(ph = phyper(q,     m,  n , k), # 1 = orig.
                p2 = phyper(q,     k, N-k, m), # swap m <-> k (keep N = m+n)
                ## "lower.tail = FALSE"  <==>  1 - p..(..)
                Ip2= phyper(m-1-q, N-k, k, m, lower.tail=FALSE),
                Ip1= phyper(k-1-q, n,   m, k, lower.tail=FALSE))

    ## check that all are (approximately) the same :
    stopifnot(all.equal(pm[,1], pm[,2], tolerance=tol),
              all.equal(pm[,2], pm[,3], tolerance=tol),
              all.equal(pm[,3], pm[,4], tolerance=tol))
    list(q = q, phyp = pm)
}

phyperBinMolenaar <-
function(q, m, n, k, lower.tail=TRUE, log.p=FALSE) {
    .Deprecated("phyperBinMolenaar.1")
    phyperBinMolenaar.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p)
}

phyperBinMolenaar.1 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(q, size = k, prob = hyper2binomP(q, m,n,k),
           lower.tail=lower.tail, log.p=log.p)
phyperBinMolenaar.2 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    ## swap k ('n') with m ('Np') -- but with R's notation n=N-Np changes too:
    pbinom(q, size = m, prob = hyper2binomP(q, k, n-k+m, m),
           lower.tail=lower.tail, log.p=log.p)
phyperBinMolenaar.3 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE) {
    ## "Ip2"
    pbinom(m-1-q, size = m, prob = hyper2binomP(m-1-q, m+n-k, k, m),
           lower.tail = !lower.tail, log.p=log.p)
    ##                 ===
}
phyperBinMolenaar.4 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE) {
    ## "Ip1"
    pbinom(k-1-q, size = k, prob = hyper2binomP(k-1-q, n, m, k),
           lower.tail = !lower.tail, log.p=log.p)
    ##                 ===
}

## Now, for completeness, also the simple binomial approximations:
phyperBin.1 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(q, size = k, prob = m/(m+n), lower.tail=lower.tail, log.p=log.p)
phyperBin.2 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    ## swap k ('n') with m ('Np') -- but with R's notation n=N-Np changes too:
    pbinom(q, size = m, prob = k/(m+n), lower.tail=lower.tail, log.p=log.p)
phyperBin.3 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(m-1-q, size = m, prob = (m+n-k)/(m+n), lower.tail = !lower.tail, log.p=log.p)
phyperBin.4 <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
    pbinom(k-1-q, size = k, prob = n/(m+n), lower.tail = !lower.tail, log.p=log.p)


phyperAllBinM <- function(m, n, k, q = .suppHyper(m,n,k),
                            lower.tail=TRUE, log.p=FALSE)
{
    cbind(pM1 = phyperBinMolenaar.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM2 = phyperBinMolenaar.2(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM3 = phyperBinMolenaar.3(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM4 = phyperBinMolenaar.4(q, m, n, k, lower.tail=lower.tail, log.p=log.p))
}

phyperAllBin <- function(m, n, k, q = .suppHyper(m,n,k),
                           lower.tail=TRUE, log.p=FALSE)
{
    cbind(
          p1 = phyperBin.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          p2 = phyperBin.2(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          p3 = phyperBin.3(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          p4 = phyperBin.4(q, m, n, k, lower.tail=lower.tail, log.p=log.p),

          pM1 = phyperBinMolenaar.1(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM2 = phyperBinMolenaar.2(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM3 = phyperBinMolenaar.3(q, m, n, k, lower.tail=lower.tail, log.p=log.p),
          pM4 = phyperBinMolenaar.4(q, m, n, k, lower.tail=lower.tail, log.p=log.p))
}

dhyperBinMolenaar <- function(x, m, n, k, log=FALSE)
    dbinom(x, size=k, prob = hyper2binomP(x, m,n,k), log=log)

### ----------- ----------- lfastchoose() etc -----------------------------

lfastchoose <- function(n,k) lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

f05lchoose <- function(n,k)
    lfastchoose(n = floor(n + .5),
                k = floor(k + .5))

## in  math/gamma.c :
p1 <- c(0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2)
## (nowhere used ??)

## The  Bernoulli Numbers:--------------------------------------------

.bernoulliEnv <- new.env(parent = emptyenv(), hash = FALSE)

Bern <- function(n, verbose = getOption("verbose", FALSE))
{
  ## Purpose: n-th Bernoulli number -- exercise in cashing
  ## ----------------------------------------------------------------------
  ## Arguments: n >= 0   B0 = 1, B1 = +1/2,  B2 = 1/6,  B4 = -1/30,...
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 26 Apr 1997  (using B1 = -1/2, back then)
  n <- as.integer(n)
  if(n < 0) stop("'n'  must not be negative")
  if(n== 0) 1 else # n == 1, now have +1/2, being compatible with Rmpfr::Bernoulli() :
  if(n==1) +1/2 else
  if(n %% 2 == 1) 0 else { ## use 'cache': .Bernoulli[n] == Bern(2 * n)
    n2 <- n %/% 2
    if(do.new <- is.null(.bernoulliEnv$.Bern))
                         .bernoulliEnv$.Bern <- numeric(0)
    if(do.new || length(.bernoulliEnv$.Bern) < n2) { # Compute  Bernoulli(n)
        if(verbose) cat("n=",n,": computing", sep='', "\n")
        Bk <- k0 <- seq(length=n2-1)
        if(n2 > 1) {
            for(k in k0) Bk[k] <- Bern(2*k)
            k0 <- 2 * k0
        }
        .bernoulliEnv$.Bern[n2] <-
            B <- - sum( choose(n+1, c(0,1,k0)) * c(1,-1/2, Bk)) / (n+1)
        B
    }
    else
        .bernoulliEnv$.Bern[n2]
  }
}


### lgammaAsymp() --- Asymptotic log gamma function :

lgammaAsymp <- function(x, n)
{
  ## Purpose: asymptotic log gamma function
  ## ----------------------------------------------------------------------
  ## Arguments: x >~ 3 ; n: number of terms in sum
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 27 Apr 97, 22:18
  stopifnot(length(n) == 1, n >= 0)
  s <- (x-1/2)*log(x) -x + log(2*pi)/2
  if(n >= 1) {
    Ix2 <- 1/(x*x)
    k <- 1:n
    Bern(2*n) #-> assigning .bernoulliEnv$.Bern[1:n]
    Bf <- rev(.bernoulliEnv$.Bern[k] / (2*k*(2*k-1)))
    bsum <- Bf[1]
    for(i in k[-1])
      bsum <- Bf[i] + bsum*Ix2
    s + bsum/x
  } else s
}


### R version of  {old} C version in  <R>/src/nmath/phyper.c
##  before Morten Welinder's (15 Apr 2004 '[Rd] phyper accuracy and efficiency (PR#6772)')
## from Gnumeric: The last of this was for R 1.9.1 i.e., R-1.9.1/src/nmath/phyper.c

## NB: Now *vectorized* in all four arguments
## --
## was     function(x, NR, NB, n)
phyperR <- function(q, m,  n, k, lower.tail=TRUE, log.p=FALSE)
{
    q <- floor(q)
    NR <- floor(m + 0.5)
    NB <- floor(n + 0.5)
    N <- NR + NB
    k <- floor(k + 0.5)
    stopifnot(NR >= 0, NB >= 0,  0 <= k, k <= N)
    xstart <- pmax(0, k - NB)
    xend   <- pmin(k, NR)
    inside <- (xstart <= q & q < xend) # << is of correct recycled length
    len <- length(inside)
    ## result
    r <- numeric(len) # = 0
    ## if(q < xstart) return(0.0)
    r[q >= xend] <- 1
    qI <- q[inside]
    ## recycle others *before* [inside]:
    NB <- rep_len(NB,    len)[inside]
    NR <- rep_len(NR,    len)[inside]
    N  <- rep_len(N ,    len)[inside]
    xr <- rep_len(xstart,len)[inside]
    k  <- rep_len(k,     len)[inside]
    xb <- k - xr
    ltrm <- lchoose(NR, xr) + lchoose(NB, xb) - lchoose(N, k)
    ##o term <- exp(ltrm)
    NR <- NR - xr
    NB <- NB - xb
    ## Sm <-
    s2 <- rep(0, length(qI))
    while(any(I <- xr <= qI)) {
        ##o Sm <- Sm + term
        s2[I] <- s2[I] + exp(ltrm[I])
        xr <- xr+1
        NB <- NB+1
        ff <- (NR * xb / (xr * NB))[I] ## (NR / xr) * (xb / NB)
        ltrm[I] <- ltrm[I] + log(ff)
        ##o term <- term * ff
        xb <- xb-1
        NR <- NR-1
    }
    r[inside] <- s2 ##o return(Sm,s2)
    ## return
    if(lower.tail)
         .D_val (r, log.p)
    else .D_Clog(r, log.p)

}

###------- pure R version of "new" Morten_Welinder--phyper() ----

## NB: Try to write this in a way to be used easily via  Rmpfr
## --  ( --> a version of this to become part of package 'DPQmpfr')

## (Taken from R (devel)source <R>/src/main/phyper.c, 2020-06-15 :

 ## * From: Morten Welinder <terra@gnome.org>
 ## * Cc: R-bugs@biostat.ku.dk
 ## * Subject: [Rd] phyper accuracy and efficiency (PR#6772)
 ## * Date: Thu, 15 Apr 2004 18:06:37 +0200 (CEST)
 ## ......

 ## The current version has very serious cancellation issues.  For example,
 ## if you ask for a small right-tail you are likely to get total cancellation.
 ## For example,  phyper(59, 150, 150, 60, FALSE, FALSE) gives 6.372680161e-14.
 ## The right answer is dhyper(0, 150, 150, 60, FALSE) which is 5.111204798e-22.

 ## phyper is also really slow for large arguments.

 ## Therefore, I suggest using the code below. This is a sniplet from Gnumeric ...
 ## The code isn't perfect.  In fact, if  x*(NR+NB)  is close to	n*NR,
 ## then this code can take a while. Not longer than the old code, though.

 ## -- Thanks to Ian Smith for ideas.


## #include "nmath.h"
## #include "dpq.h"

## C code args:     x, NR, NB, n,
pdhyper <- function(q,  m,  n, k, log.p = FALSE,
                     epsC = .Machine$double.eps, verbose = getOption("verbose")) {
## Calculate
##
##        phyper (q, m, n, k, TRUE, FALSE)
## [log]  ----------------------------------
##           dhyper (q, m, n, k, FALSE)
##
## without actually calling phyper.  This assumes that
##
##    q * (m + n) <= k * m   <==>   q/k  <=  m / (m+n)

### For now:
    stopifnot(length(q) == 1,
              q == floor(q),
              length(m) == 1, length(n) == 1, length(k) == 1)
    ## C code uses LDOUBLE (= long double) which we can't in R.
    sum  <- 0 # LDOUBLE sum = 0;
    term <- 1 # LDOUBLE term = 1;
    if(verbose) q0 <- q
    while (q > 0 && term >= epsC * sum) {
	term <- term * (q * (n - k + q) / (k + 1 - q) / (m + 1 - q))
	sum  <- sum + term
	q <- q-1
    }
    if(verbose)
        message("pdhyper(q=",q0,"): used q0-q = ", q0-q, " while(.) iterations")

    if(log.p) log1p(sum) else 1 + sum
}


## C code args:      x, NR, NB, n,
phyperR2 <- function(q,  m,  n, k, lower.tail=TRUE, log.p=FALSE, ...)
{
## Sample of  k balls from  m red  and	 n black ones;	 q are red

### For now:
    stopifnot(length(q) == 1,
              ## q == floor(q),
              length(m) == 1, length(n) == 1, length(k) == 1)

    if(is.na(q) || is.na(m) || is.na(n) || is.na(k))
	return(q + m + n + k)

    q  <- floor (q + 1e-7)
    m <- round(m) ## round(.)  was "nmath.h"s  R_forceint()
    n <- round(n)
    k  <- round(k)

    if(m < 0 || n < 0 || !is.finite(m + n) || k < 0 || k > m + n) {
	warning("Invalid values for (m, n, k)")
        return(NaN)
    } else if(q * (m + n) > k * m) { ##  Swap tails.
	oldn <- n; n <- m; m <- oldn
	q <- k - q - 1
	lower.tail <- !lower.tail
    }

    ## support of dhyper() as a function of its parameters
    ##   .suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)
    if (q < 0 || q < k - n)
	return(.DT_0(lower.tail, log.p))
    if (q >= m || q >= k)
	return(.DT_1(lower.tail, log.p))

    d <- dhyper (q, m, n, k, log.p)
   ## dhyper(.., log.p=FALSE) > 0 mathematically, but not always numerically :
    if((!log.p && d == 0.) ||
        (log.p && d == -Inf))
	return(.DT_0(lower.tail, log.p))

    pd <- pdhyper(q, m, n, k, log.p, ...)
    ## return
    if(log.p)
        .DT_Log (d + pd, lower.tail)
    else .D_Lval(d * pd, lower.tail)
}


