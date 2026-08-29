### find (double precision accurate expm1x(x) :=  e^x - 1 - x   for small |x|
### --                                  using taylor e^x = \sum_{k=0}^n x^n / n!
##  expm1x(x) :=  e^x - 1 - x = expm1(x) - x =
##    = x^2/2! + x^3/3! + x^4/4! + ....
##    = x^2/2 (1 + x/3(1 + x/4(1 + x/5(1 + x/6(....) ))))
expm1x.0 <- function(x) exp(x) -1 - x
expm1x.1 <- function(x) expm1(x)  - x
expm1xTser <- function(x, k) { ## perfectly vectorized in x
    stopifnot(length(k) == 1L, k == (k <- as.integer(k)), 1L <= k)
    if(k <= 9)
       switch(k   # k = 1, 2, 3 :
         , x^2/2
         , x^2/2*(1 + x/3)
         , x^2/2*(1 + x/3*(1 + x/4))
           ## k = 4, 5, 6 :
         , x^2/2*(1 + x/3*(1 + x/4*(1 + x/5)))
         , x^2/2*(1 + x/3*(1 + x/4*(1 + x/5*(1 + x/6))))
         , x^2/2*(1 + x/3*(1 + x/4*(1 + x/5*(1 + x/6*(1 + x/7)))))
           ## k = 7, 8, 9 :
         , x^2/2*(1 + x/3*(1 + x/4*(1 + x/5*(1 + x/6*(1 + x/7*(1 + x/8))))))
         , x^2/2*(1 + x/3*(1 + x/4*(1 + x/5*(1 + x/6*(1 + x/7*(1 + x/8*(1 + x/9)))))))
         , x^2/2*(1 + x/3*(1 + x/4*(1 + x/5*(1 + x/6*(1 + x/7*(1 + x/8*(1 + x/9*(1 + x/10))))))))
           )
    else { # k >= 10 .. typically when x is mpfr, in general when x has *higher* than double precision
        p <- x / (k+1L)
        for(j in k:2L)
            p <- (1 + p)*x/j
        x * p
    }
}


##' Accurate expm1x(x) :=  e^x - 1 - x   (notably for small |x| )
##'                     = x^2/2! + x^3/3! + x^4/4! + ....
##'                     = x^2/2 (1 + x/3(1 + x/4(1 + x/5(1 + x/6(....) ))))
### Experimental version w/  'eps'   --------------  "production" version for double prec ==> expm1x() below
expm1xXp <- function(x, eps = .Machine$double.eps,
                   verbose = FALSE,
                   x1bnd = 1.05) # see experiments below (w/ x1bnd = 9 to *see* where default should be)
{
    stopifnot(length(x) == 1L, length(eps) == 1L, 0 < eps, eps < 1)
    if((ax <- abs(x)) >= x1bnd)  return(expm1(x) - x)
    P <- -log2(eps) # > 0; typically >= 52
    ## Using Taylor series order (k), relative error is <=  |(x^(k+1) / (k+1)!) / (x^2/2)| =  2 |x|^{k-1} / (k-1)!  <= eps = 2^-P
    ##   2 |x|^(k-1) / (k-1)!  <=  2^-P
    ## ln2 + (k-1)*log(|x|)    <=  - P * ln2 + log[(k-1)!]
    ##      (k-1)*log(|x|)    <=  -(P+1)*ln2 + log[(k-1)!]
    ##      log(|x|)    <=  (-(P+1)*ln2 + log[(k-1)!])/(k-1) ;  log((k-1)!) == log(Gamma(k)) = lgamma(k)
    ln2 <- log(2) ## log(if(inherits(x, "mpfr")) mpfr(2, .getPrec(x)) else 2)
    C <- -ln2*(P+1)
    lx <- log(ax)
    f0 <- function(k) (C + lgamma(k))/(k-1) - lx
    ur <- uniroot(f0, c(1.25, 100))# "FIXME": better interval
    ##    -------
    if(verbose) str(ur, digits.d = 12)
    k <- ceiling(k. <- ur$root)
    structure(expm1xTser(x, k = k), k = k.)
}

if(FALSE) {
    op <- options(digits = 4) ; require(Rmpfr)
    x <- .5 ; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x)) / eM - 1) #  3.182e-16 -5.511e-17
    x <- .75; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x)) / eM - 1) #  3.153e-16 1.277e-17
    x <- 1.0; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# -2.013e-16 -4.670e-17
    x <- 1.05;eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# -7.881e-17  5.866e-17
    ##----------------------------   about from here, expm1x.1(x) := expm1(x) - x  is fully accurate
    x <- 1.1; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# 1.76e-17    1.76e-17
    x <- 1.25;eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# 3.712e-17   3.712e-17
    x <- 1.5; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# 7.028e-17  -4.177e-17
    x <- 2.0; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# 4.095e-17   4.095e-17
    x <- 3.0; eM <- expm1x.1(mpfr(x, 256)); asNumeric(c(expm1x.1(x), expm1xXp(x, x1bnd=9))/eM -1)# 1.136e-17   1.136e-17
}

## *hidden* --- more experiments in ../tests/expm1x-tst.R <<<<<<<<<<<<<

##' Find the x_k boundaries (for each k = 2, ... k.max)
##' for a given binary precision P = -log2(eps) <==> eps = 2^-P
##'  ==================> k = *order* of Taylor approx;
##' "k" = order 2  <===> expm1xTser(*, k = 1)
expm1xBnds <- function(ord, P = 52L, tol = 1e-7,
                           interval = c(-100, 5)) # FIXME?  smarter interval
{
    stopifnot(length(ord) == 1L, ord == (ord <- as.integer(ord)), ord >= 2L,
              length(P) == 1L, (P <- as.integer(P)) >= 16L, length(interval) == 2L)
    C <- -log(2)*(P+1L)
    Ckp <- (C + lgamma(ord))/(ord-1L)
    f0 <- function(lx) Ckp - lx # lx = log(abs(x))
    uniroot(f0, interval=interval, extendInt = "yes", tol = tol)
}

if(FALSE)  {
    op <- options(digits = 4) ; require(Rmpfr)
  if(FALSE)
      op <- c(op, options(error=recover))
xBnds <- lapply(2:18, expm1xBnds)
## op <- c(op, options(digits = 5)) # op has 'digits' already
dBnds <- do.call(rbind, lapply(xBnds, data.frame)) # <- nice data.frame from the  uniroot() lists
dBnds <- cbind(ord = 2:18, within(dBnds, { rm(init.it); x. <- exp(root) }))
dBnds #
##    ord     root     f.root iter estim.prec        x.   "empirically sufficient" (relErr..) -- see "k = 1" below
## 1    2 -36.7368  7.105e-15    2  5.000e-08 1.110e-16    5 *  10^{-16}
## 2    3 -18.0218 -3.553e-15    2  5.000e-08 1.490e-08    4.4* 10^{-8}
## 3    4 -11.6483  1.776e-15    2  5.000e-08 8.733e-06    2.02*10^{-5}
## 4    5  -8.3897 -1.776e-15    2  5.000e-08 2.272e-04    3.95*10^{-4}
## 5    6  -6.3899  1.776e-15    2  5.000e-08 1.678e-03    3 *  10^{-3}
## 6    7  -5.0263 -8.882e-16    2  5.000e-08 6.563e-03     0.0111
## 7    8  -4.0302 -2.665e-15    2  5.000e-08 1.777e-02     0.031
## 8    9  -3.2665  8.882e-16    2  5.000e-08 3.814e-02     0.06
## 9   10  -2.6594  8.882e-16    2  5.000e-08 6.999e-02     0.10
## 10  11  -2.1632  1.332e-15    2  5.000e-08 1.150e-01     0.185
## 11  12  -1.7486 -2.220e-16    2  5.000e-08 1.740e-01     0.27
## 12  13  -1.3958 -2.220e-16    2  5.000e-08 2.476e-01     0.385
## 13  14  -1.0911 -1.110e-15    2  5.000e-08 3.358e-01     0.49
## 14  15  -0.8247  0.000e+00    1  9.918e+01 4.384e-01     0.63
## 15  16  -0.5892  2.220e-16    2  5.000e-08 5.548e-01     0.79
## 16  17  -0.3791  2.220e-16    2  5.000e-08 6.845e-01     1.12
## 17  18  -0.1901  3.053e-16    2  5.000e-08 8.269e-01     1.25
##							   ========== take (some of) these, as cutoffs 'cutx' below

} # if(FALSE)

expm1x <- function(x, cutx = c(  4.4e-8, 0.10, 0.385, 1.1, 2),  # cutoff x[k]
                         k = c( 2,      9,   12,   17))
{
    nk <- length(k)
    stopifnot(is.numeric(cutx), is.numeric(k), k == (k <- as.integer(k)), nk >= 2, nk + 1 == length(cutx),
              !is.unsorted(k), !is.unsorted(cutx), 0 < cutx, is.finite(cutx))
    ## vectorizing using findInterval():
    in.x <- findInterval(abs(x), c(0, cutx, Inf), all.inside = TRUE)
    iLst <- split(seq_along(x), in.x)
    r <- x
    for(nmi in names(iLst)) {
        jk <- as.integer(nmi)
        i <- iLst[[nmi]]
        xi <- x[i]
        r[i] <-
            if(jk <= nk)
                expm1xTser(xi, k = k[jk])
            else if(jk == nk+1L)
                expm1(xi) - xi
            else # if(jk == nk+2L)
                exp(xi) -1 -xi
    }
    ## return
    r
}

## Experiments changing cutoff etc -- see  "Older Experiments for finding cutoffs ..
## ===>  ../tests/expm1x-tst.R
##       -- ----- ^^^^^^^^^^^^

