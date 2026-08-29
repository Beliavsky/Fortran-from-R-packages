### dgamma() --  Catherine Loader's code in R ------

## for now: just the plan
## ~/R/D/r-devel/R/src/nmath/dgamma.c :

## R_D__0 : <==>  (if(log) -Inf else 0)
## R_D__1 : <==>  (if(log) 0 else 1) <==>  !log

dgamma.R <- function(x, shape, scale = 1, log = FALSE, dpois_r_args = list())
{
    ## vectorized in 'x' only
    stopifnot(length(shape) == 1, length(scale) == 1, length(log) == 1)
    if (is.na(shape) || is.na(scale))
        return (x + shape + scale)
    if (shape < 0 || scale <= 0)
        stop("invalid 'shape < 0' or 'scale <= 0'")
    r <- x # result
    if (N <- anyNA(x))
        r[N] <- x
    ok <- !N
    if (any(neg <- ok & (x < 0)))
        r[neg] <- if(log) -Inf else 0
    if (shape == 0) { ##/* point mass at 0 */
	r[ok & x == 0] <- Inf
	r[ok & x != 0] <- if(log) -Inf else 0
        return(r)
    }
    ## else: scale, shape > 0
    ## simplified --- "hoping" dpois_raw() would work ...
    ## } else if (x == 0) {
    ##     if (shape < 1)
    ##         Inf
    ##     else if (shape > 1) {
    ##         if(log) -Inf else 0
    ##     } else
    ##         if(log) -log(scale) else 1 / scale
    ## } else
    if(length(iOk <- which(ok & !neg))) {
        x <- x[iOk]
        r[iOk] <-
            if (shape < 1) {
                pr <- if(length(dpois_r_args))
                          do.call(dpois_raw, c(list(x = shape, lambda = x/scale, log = log),
                                               dpois_r_args))
                      else
                          dpois_raw(shape, x/scale, log)
                ## NB: currently *always*  shape/x > 0  if shape < 1:
                ## -- overflow to Inf happens, but underflow to 0 does NOT :
                if(log) pr + ifelse(shape/x == Inf, log(shape)-log(x), log(shape/x))
                else pr*shape/x
            }
            else {##  shape >= 1
                pr <- if(length(dpois_r_args))
                          do.call(dpois_raw, c(list(x = shape-1, lambda = x/scale, log = log),
                                               dpois_r_args))
                      else
                          dpois_raw(shape-1, x/scale, log)
                if(log) pr - log(scale) else pr/scale
            }
    }
    ## return
    r
}

## not clear if we should *ever* use this {OTOH, factorial(x) can be 100% accurate, e.g., for MPFR
##          ( <==>  Rmpfr's  factorialMpfr() and gmp factorialZ() )
dpois_simpl0 <- function(x, lambda, log=FALSE)
    .D_val(exp(-lambda) * lambda^x / factorial(x), log)

## when  x & lambda  are different orders of magnitude, this may be more accurate than anything:
dpois_simpl <- function(x, lambda, log=FALSE)
    .D_exp(-lambda + (x*log(lambda) - lgamma(x+1)), log) # use lgamma1pC(x) for x << 1  (?)

## ~/R/D/r-devel/R/src/nmath/dpois.c   --> dpois_raw()
## // called also from dgamma.c, pgamma.c, dnbeta.c, dnbinom.c, dnchisq.c :
dpois_raw <- function(x, lambda, log=FALSE,
                      version = c("bd0_v1", "bd0_p1l1d", "bd0_p1l1d1", "bd0_p1l1", "bd0_l1pm",
                                  ## NB:      ebd0_v1 R version *broken*
                                  "ebd0_C1", "ebd0_v1"),
                      ## future ?! version = c("ebd0_v1", "bd0_v1"),
                      small.x__lambda = .Machine$double.eps,
                      bd0.delta = 0.1,
                      ## optional arguments of log1pmx() :
                      tol_logcf = 1e-14, eps2 = 0.01, minL1 = -0.79149064, trace.lcf = verbose,
                      logCF = if (is.numeric(x)) logcf else logcfR,
                      verbose = FALSE)
{
    ##  x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
    ##     lambda >= 0
    stopifnot(is.logical(log), length(log) == 1)
    R_D__0 <- if(log) -Inf else 0
    R_D__1 <- !log #  (if(log) 0 else 1)
    if(!(nx <- length(x)) || !(nl <- length(lambda))) # 0-length result if one {x, lambda} is:
        return(x[0])
    M <- max(nx, nl)
    r <- rep_len(x+0*lambda, M) # result of same number class as 'x+lambda' (e.g. "mpfr")

    if(verbose) {
        cat(sprintf("dpois_raw(): max length M = %d\n", M))
        inds <- function(i)
            if((n <- length(i)) <= 4)
                paste(i, collapse=", ")
            else paste(paste(i[1:3], collapse=", "), "..", i[n])
    }
    if(M == 0) return(r)
    stopifnot(length(small.x__lambda) == 1, small.x__lambda >= DBL_MIN)

    verb1 <- pmax(0L, verbose - 1L)

    if(M > 1L) { # do recycle:
        x      <- rep_len(x,      M)
        lambda <- rep_len(lambda, M)
    }
    if(any(B <- lambda == 0))
        r[B] <- ifelse(x[B] == 0, R_D__1, R_D__0)
    BB <- !B # BB is true "for remaining cases"
    if(any(B <- BB & !is.finite(lambda))) {
        r[B] <- R_D__0 ## including for the case where  x = lambda = +Inf
        BB <- BB & !B
    }
    if(any(B <- BB & x < 0)) {
        r[B] <- R_D__0
        BB <- BB & !B
    }
    if(any(BB & x <= lambda * small.x__lambda)) {
        if(any(B <- BB & x <= lambda * DBL_MIN)) {
            if(verbose & any(i <- B & x > 0))
                cat(sprintf(" extremely small x/lambda in [%g,%g]\n for x[%s]\n",
                            min(x[i]), max(x[i]), inds(which(i))))
            r[B] <- .D_exp(-lambda[B], log)
        } else { ## x << lambda (less extreme)  x/lambda  in (DBL_MIN, small.x__lambda]
            ## ebd0(x, lambda) not good, and there's no bad cancellation here
            B <- !B
            if(verbose)
                cat(sprintf(" very small x/lambda in [%g,%g]\n for x[%s]\n",
                            min(x[B]), max(x[B]), inds(which(B))))
	    r[B] <- .D_exp(-lambda[B] + (x[B]*log(lambda[B]) -lgamma(x[B]+1)), log)
        }
        BB <- BB & !B
    }
    if(any(B <- BB & lambda < x * DBL_MIN)) {
	if(any(B. <- !is.finite(x[B]))) ## lambda < x = +Inf
	    r[B][B.] <- R_D__0
        lamb <- lambda[B][!B.]
        x.   <- x     [B][!B.]
        if(verbose && length(lamb))
           cat(sprintf(" very small lambda/x in [%g,%g]\n for x[%s]\n",
                       min(lamb/x.), max(lamb/x.), inds(which(B)[!B])))
        r[B][!B.] <- .D_exp(-lamb + x.*log(lamb) -lgamma(x.+1), log)
        BB <- BB & !B
    }
    if(any(BB)) { ## else remaining "main branch"
        version <- match.arg(version)
        x <- x[BB]
        lambda <- lambda[BB]
        BL <- x >= 2^1023 / pi # really large x ( <==> 2*pi*x  overflows to Inf)
        D_fE <- function(x, T0, log) {
            r <- x
            r[ BL] <- .D_rtxp(sqrt(2*pi)*sqrt(x[BL]), -T0[ BL], log) # avoid overflow for really large x
            r[!BL] <- .D_fexp(     2*pi *    x[!BL],  -T0[!BL], log)
            r
        }
        del.x <- stirlerr(x, verbose=verb1) # = \delta(x)
        ## version = c("bd0_v1", "bd0_p1l1d", "bd0_p1l1d1", "bd0_l1pm", "ebd0_v1"),
        r[BB] <- switch(version,
           "bd0_v1"    = D_fE(x, del.x + bd0(x, lambda, delta=bd0.delta, verbose=verb1), log),
           "bd0_p1l1d" = D_fE(x, del.x + bd0_p1l1d (x, lambda, tol_logcf=tol_logcf, eps2=eps2,
                                                    minL1=minL1, trace.lcf=trace.lcf, logCF=logCF),
                              log),
           "bd0_p1l1d1"= D_fE(x, del.x + bd0_p1l1d1(x, lambda, tol_logcf=tol_logcf, eps2=eps2,
                                                    minL1=minL1, trace.lcf=trace.lcf, logCF=logCF),
                              log),
           "bd0_p1l1"  = D_fE(x, del.x + bd0_p1l1  (x, lambda, tol_logcf=tol_logcf, eps2=eps2,
                                                    minL1=minL1, trace.lcf=trace.lcf, logCF=logCF),
                              log),
           "bd0_l1pm"  = D_fE(x, del.x + bd0_l1pm  (x, lambda, tol_logcf=tol_logcf, eps2=eps2,
                                                    minL1=minL1, trace.lcf=trace.lcf, logCF=logCF),
                              log),
           "ebd0_v1" = , "ebd0_C1" = {
               yM <- switch(version,
                            "ebd0_v1" = ebd0 (x, lambda, verbose=verb1,
                                              tol_logcf=tol_logcf, eps2=eps2, minL1=minL1,
                                              trace.lcf=trace.lcf, logCF=logCF),
                            "ebd0_C1" = ebd0C(x, lambda, verbose=verb1),
                            stop("internal version error:", version))
               yl <- yM[["yl"]] + del.x
               ## return
               if(log) ## "FIXME" ifelse() is inefficient (and fails for MPFR!)
                   -yl - yM[["yh"]] - ifelse(BL,
                                            0.5 * (log(2*pi)+log(x)),
                                            0.5 * log(2*pi*x))
               else
                   exp(-yl) * exp(-yM[["yh"]]) / ifelse(BL, sqrt(2*pi)*sqrt(x), sqrt(2*pi*x))
           },
           stop("invalid 'version'": version))
    }
    r
}


## ~/R/D/r-devel/R/src/nmath/bd0.c     --> bd0()
## /* Martin Plummer (in priv. e-mail):   t := (x-M)/M  ( <==> 1+t = x/M  ==>
##  *
##  * bd0 = M*[ x/M * log(x/M) + 1 - (x/M) ] = M*[ (1+t)*log1p(t) + 1 - (1+t) ]
##  *     = M*[ (1+t)*log1p(t) - t ] =: M * p1log1pm(t) =: M * p1l1(t)
##  *
## /* Morten Welinder -- when providing ebd0() in 2014:
##  *
##  * Compute x * log (x / M) + (M - x)
##  * aka    -x * log1pmx ((M - x) / x)
##  */
## double p1log1pm(double x) {
## }

bd0_p1l1d1 <- function(x, M, tol_logcf = 1e-14, ...) {
    t <- (x-M)/M
    t <- M * (log1pmx(t, tol_logcf=tol_logcf, ...) + t*log1p(t))# == M * p1l1p(t, ...)
    t[is.infinite(x) & is.finite(M)] <- Inf
    t
}

bd0_p1l1d <- function(x, M, tol_logcf = 1e-14, ...) {
    d <- x-M
    t <- d/M
    fi <- is.finite(dlg1 <- d*log1p(t))
    ## slightly better(?) , *not* computing M*t = M*((x-M)/M) for 2nd term:
    d[ fi] <- (M *  log1pmx(t, tol_logcf=tol_logcf, ...))[fi] + dlg1[fi]
    d[!fi] <- (M * (log1pmx(t, tol_logcf=tol_logcf, ...)  +  t*log1p(t)))[!fi]
    d[is.infinite(x) & is.finite(M)] <- Inf
    d
}

bd0_p1l1 <- function(x, M, ...) {
    t <- (x-M)/M
    t <- M * p1l1(t, ...) # dbl.prec "optimized" p1l1(t) == log1pmx(t) + t*log1p(t)  outside [-.07,.07]
    ##       ------- see below
    t[is.infinite(x) & is.finite(M)] <- Inf
    t
}


##' Version mentioned by Morten Welinder in PR#15628  --- ?? (not really ?? somewhere else?
bd0_l1pm <- function(x, M, tol_logcf = 1e-14, ...) {
    r <- s <- (M-x)/x
    I <- is.infinite(x) & is.finite(M) # MM, 2025-05-14
    r[I] <- s[I] <- -1
    ## L: |s| is large  -- NB, for x = 0, s=Inf (or NaN) ==> use  M * D~(x/M) there
    L <- !is.na(s) & abs(s) > 1e10
    r[L] <- M*(1 + {a <- u <- (x/M)[L]; p <- u > 0; a[p] <- u[p] * (log(u[p]) - 1); a })
    ok <- !L   # vvvvv (recycling x)
    r[ok]  <- - (x+0*M)[ok] * log1pmx(s[ok], tol_logcf=tol_logcf, ...)
    r
}

bd0 <- function(x, np,
                delta = 0.1, maxit = max(1L, as.integer(-1100 / log2(delta))),
                s0 = .Machine$double.xmin, # = DBL_MIN
                verbose = getOption("verbose"))
{
    stopifnot(0 <= delta, delta < 1, maxit >= 1L)
    if(min(n <- length(x), n2 <- length(np)) == 0) # length 0 result
        return(x[0])
    if(n > n2)
        np <- rep_len(np, n)
    else if(n < n2)
        x <- rep_len(x, n2)
    N <- as.numeric # to be used in verbose / warning printing -- also for <mpfr> !
    if(useMpfr <- inherits(x, "mpfr") || inherits(np, "mpfr")) {
        if(!requireNamespace("Rmpfr"))
            stop("need to  install.packages(\"Rmpfr\") ")
        ldexp <- function(f, E)
            if(inherits(f, "mpfr")) Rmpfr::ldexpMpfr(f, E) else DPQ::ldexp(f, E)
    }
    bd0.1 <- function(x, np) {
        if(!is.finite(x) || !is.finite(np) || np == 0.0) {
            if((isTRUE(np < x) && x == Inf) || # << no NA/NaN in {np,x}
               (isTRUE(x < np) && np== Inf))
                return(Inf)
            ## ML_ERR_return_NAN;
            warning(sprintf("invalid argument(s) in x=%g, np=%g", x, np))
            return(NaN)
        }
        d <- x-np
        ##       ____
        if(abs(d) <= delta * (x+np)) { #  '<='  was  '<' ;  '<=' allows to use delta=0
            ##   ~~~~
            v <- d/(x+np)
            if(v == 0 && d != 0) { # had underflow --> do as  bd0C(.., version = "R_2025_0510")
                if(verbose)
                    cat(sprintf(
                        "bd0(%g, %g): Initial v == 0 --> rescaling with F = %g\n",
                        N(x), N(np), N(F)))
                ## could do better: replace  v * N  with  sv * N * sv  where sv := sqrt(v)
                x. <- ldexp(x, -2L)
                n. <- ldexp(np,-2L)
                v <- (x. - n.)/(x. + n.)
            }
            s <- d*v
            if(abs(s) < s0) {
                if(verbose)
                    cat(sprintf(
                        "bd0(%g, %g): Very small initial |s| = |%g| < s0=%g -> bd0:=s.\n", 
                        N(x), N(np), N(s), N(s0)))
                return(s)
            }
            s <- ldexp(s, -1L) # such that ej = v*x (instead of 2*v*x) does *not* overflow
            ej <- x*v
            v  <- v*v # // "v = v^2"
            for (j in seq_len(maxit-1L)) { #/* Taylor series; 1000: no infinite loop
                                        # as |v| < .1,  v^2000 is "zero" */
                ej <- ej* v ##/ = x v^(2j+1)
                s_ <- s
                s  <- s + ej/(2*j+1)
                if (s == s_) { ##/* last term was effectively 0 */
                    if(verbose)
                        cat(sprintf("bd0(%g, %g): T.series w/ %d terms -> bd0=%g\n",
                                    N(x), N(np), j, N(ldexp(s, 1L))))
                    return(ldexp(s, 1L))
                }
            }
            warning(gettextf(
                "bd0(%g, %g): T.series failed to converge in %d it.; s=%g, ej/(2j+1)=%g\n -- should *NOT* happen, please report!",
                N(x), N(np), maxit, N(ldexp(s, 1L)), N(ej/(2*maxit+1))),
                domain=NA)
        }
        ## else   |x - np|  is not too small (or the iterations failed !)
        if(verbose) cat(sprintf("bd0(x=%g, np=%g): |x-np| not too \"small\" -> direct formula\n",
                                N(x), N(np)))
        l.x.n <-
            if(x/np == Inf) # overflow
                log(x) - log(np)
            else log(x/np)

        if(x > np) x*(l.x.n -1.) + np
        ##         ==================
        else       x* l.x.n + np -x
        ##         ================
    } # {bd0.1}

    if(useMpfr) # Vectorize(.) returns <mpfr1> (not "mpfr")
        Rmpfr::mpfr(Vectorize(bd0.1, c("x","np"))(x, np))
    else            Vectorize(bd0.1, c("x","np"))(x, np)
} ## {bd0}

##' The version calling into C (to the "same" code R's Rmathlib bd0() uses
bd0C <- function(x, np,
                 delta = 0.1, maxit = 1000L,
                 ##s0 = .Machine$double.xmin, # = DBL_MIN
                 version = c("R4.0", "R_2025_0510"),
                 verbose = getOption("verbose"))
{
    iVer <- pmatch(match.arg(version), eval(formals()$version))
    .Call(C_dpq_bd0, x, np, delta, maxit, iVer, verbose)# >> ../src/bd0.c
}

ebd0C <- function(x, M, verbose = getOption("verbose"))
    data.frame(.Call(C_dpq_ebd0, x, M, verbose))# >> ../src/bd0.c  {now returns *list* {case of __long__ vectors}}



##' bd0(x,M) = M * p1l1(t),  t := (x-M)/M
##' -------        =======
##'
##' p1l1(t) = (t+1)*log(1+t) - t   -- the "naive" :
p1l1. <- function(t) (t+1)*log1p(t) - t

##' As we found, this is *MUCH* better -- actually almost perfect
##' @param ... notably  'tol_logcf = 1e-14' is default in log1pmx()
p1l1p <- function(t, ...) log1pmx(t, ...) + t*log1p(t)


##' p1l1 Taylor series approximation (there must be a faster converging one with quadratic terms !???!
##  ~~~~        ~~~
.p1l1ser <- function(t, k, F = t^2/2) {
    two <- if(is.numeric(t)) 2 else 0*t[1] + 2
    if(k == 1L)
        return(F + t*0) # incl recycling
    ## else k >= 2 :
    a <- two/((k+1)*k)
    for(n in k:2)
        a <- two/(n*(n-1)) - t*a
    F * a
}

p1l1ser <- function(t, k, F = t^2/2) {
    if(!length(t)) return(t) # else can use t[1]
    stopifnot(k == (k <- as.integer(k)), k >= 1)
    if(k <= 12) { # i(.) for "inverse"
        i <- if(is.numeric(t))
                   function(N) 1/N
               else { # one := 1, but should inherit 'numeric class', e.g., "mpfr" from t
                   one <- t[1]*0 + 1
                   function(N) one/N
               }
        switch(k
               , F			# k = 1
               , F*(1 - t/3)		# k = 2
               , F*(1 - t*(i(3) - t/6))	# k = 3
               , F*(1 - t*(i(3) - t*(i(6) - t/10)))			# k = 4
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t/15))))		# k = 5
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t/21)))))	# k = 6
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t*(i(21) - t/28)))))) # k = 7
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t*(i(21) - t*(i(28) - t/36))))))) # k = 8
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t*(i(21) - t*(i(28) - t*(i(36) - t/45)))))))) # k = 9
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t*(i(21) - t*(i(28) - t*(i(36) - t*(i(45) - t/55))))))))) # k = 10
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t*(i(21) - t*(i(28) - t*(i(36) - t*(i(45) - t*(i(55) - t/66)))))))))) # k = 11
               , F*(1 - t*(i(3) - t*(i(6) - t*(i(10) - t*(i(15) - t*(i(21) - t*(i(28) - t*(i(36) - t*(i(45) - t*(i(55) - t*(i(66) - t/78))))))))))) # k = 12
               )
    } else {
        .p1l1ser(t, k, F)
    }
}

p1l1 <- function(t, F = t^2/2, ...)  {
    r <- t
    if(!is.numeric(t)) warning("non-numeric 't' -- approximation intervals are for double precision numeric")
    for(i in seq_along(r)) { # (manual vectorization)
        t_ <- t[i] # with t_ := t[i], compute r[i] := F[i] * p1l1(t_) / (t_^2/2)  = p1l1(t_) {for default F[i]}
        r[i] <-
            if(-.0724 <= t_ & t_ <= .0718) ## NB: depends on last case dealt with below
              F[i]*( ## the cutoffs are  t.0  values see the corrDig[] table computations
                if(-6.60e-16 <= t_ & t_ <= 6.69e-16)
                    1				# k = 1
                else if(-3.65e-8 <= t_ & t_ <= 3.65e-8)
                    1 - t_/3			# k = 2
                else if(-1.30e-5 <= t_ & t_ <= 1.32e-5)
                    1 - t_*(1/3 - t_/6) 	# k = 3
                else if(-2.39e-4 <= t_ & t_ <= 2.42e-4)
                    1 - t_*(1/3 - t_*(1/6 - t_/10)) # k = 4
                else if(-1.35e-3 <= t_ & t_ <= 1.38e-3)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_/15))) # k = 5
                else if(-4.27e-3 <= t_ & t_ <= 4.34e-3)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_/21)))) # k = 6
                else if(-9.60e-3 <= t_ & t_ <= 9.78e-3)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_*(1/21 - t_/28))))) # k = 7
                else if(-.0178 <= t_ & t_ <= .0180)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_*(1/21 - t_*(1/28 - t_/36)))))) # k = 8
                else if(-.0285 <= t_ & t_ <= .0285)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_*(1/21 - t_*(1/28 - t_*(1/36 - t_/45))))))) # k = 9
                else if(-.0413 <= t_ & t_ <= .0414)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_*(1/21 - t_*(1/28 - t_*(1/36 - t_*(1/45 - t_/55)))))))) # k = 10
                else if(-.0562 <= t_ & t_ <= .0564)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_*(1/21 - t_*(1/28 - t_*(1/36 - t_*(1/45 - t_*(1/55 - t_/66))))))))) # k = 11
                else if(-.0724 <= t_ & t_ <= .0718)
                    1 - t_*(1/3 - t_*(1/6 - t_*(1/10 - t_*(1/15 - t_*(1/21 - t_*(1/28 - t_*(1/36 - t_*(1/45 - t_*(1/55 - t_*(1/66 - t_/78)))))))))) # k = 12
              )
            else if(missing(F)) ## default F = t^2 / 2, i.e. "direct formula" directly
                              log1pmx(t_, ...) + t_ * log1p(t_)
            else # F not missing, relatively large |t_|
                2*F[i]/t_^2 * log1pmx(t_, ...) + t_ * log1p(t_)
        ## was:  (t_+1)*log1p(t_) - t_
    }
    r
}


## ebd0(): R Bugzilla PR#15628 -- proposed accuracy improvement by Morten Welinder

## /*
##  * A table of logs for scaling purposes.  Each value has four parts with
##  * 23 bits in each.  That means each part can be multiplied by a double
##  * with at most 30 bits set and not have any rounding error.  Note, that
##  * the first entry is log(2).
##  *
##  * Entry i is associated with the value r = 0.5 + i / 256.0.  The
##  * argument to log is p/q where q=1024 and p=floor(q / r + 0.5).
##  * Thus r*p/q is close to 1.
##  */
logf_mat <- ## 'bd0_scale' in C
    matrix(nrow = 4,
       c(
	 +0x1.62e430p-1, -0x1.05c610p-29, -0x1.950d88p-54, +0x1.d9cc02p-79 , # /* 128: log(2048/1024.) */
	 +0x1.5ee02cp-1, -0x1.6dbe98p-25, -0x1.51e540p-50, +0x1.2bfa48p-74 , # /* 129: log(2032/1024.) */
	 +0x1.5ad404p-1, +0x1.86b3e4p-26, +0x1.9f6534p-50, +0x1.54be04p-74 , # /* 130: log(2016/1024.) */
	 +0x1.570124p-1, -0x1.9ed750p-25, -0x1.f37dd0p-51, +0x1.10b770p-77 , # /* 131: log(2001/1024.) */
	 +0x1.5326e4p-1, -0x1.9b9874p-25, -0x1.378194p-49, +0x1.56feb2p-74 , # /* 132: log(1986/1024.) */
	 +0x1.4f4528p-1, +0x1.aca70cp-28, +0x1.103e74p-53, +0x1.9c410ap-81 , # /* 133: log(1971/1024.) */
	 +0x1.4b5bd8p-1, -0x1.6a91d8p-25, -0x1.8e43d0p-50, -0x1.afba9ep-77 , # /* 134: log(1956/1024.) */
	 +0x1.47ae54p-1, -0x1.abb51cp-25, +0x1.19b798p-51, +0x1.45e09cp-76 , # /* 135: log(1942/1024.) */
	 +0x1.43fa00p-1, -0x1.d06318p-25, -0x1.8858d8p-49, -0x1.1927c4p-75 , # /* 136: log(1928/1024.) */
	 +0x1.3ffa40p-1, +0x1.1a427cp-25, +0x1.151640p-53, -0x1.4f5606p-77 , # /* 137: log(1913/1024.) */
	 +0x1.3c7c80p-1, -0x1.19bf48p-34, +0x1.05fc94p-58, -0x1.c096fcp-82 , # /* 138: log(1900/1024.) */
	 +0x1.38b320p-1, +0x1.6b5778p-25, +0x1.be38d0p-50, -0x1.075e96p-74 , # /* 139: log(1886/1024.) */
	 +0x1.34e288p-1, +0x1.d9ce1cp-25, +0x1.316eb8p-49, +0x1.2d885cp-73 , # /* 140: log(1872/1024.) */
	 +0x1.315124p-1, +0x1.c2fc60p-29, -0x1.4396fcp-53, +0x1.acf376p-78 , # /* 141: log(1859/1024.) */
	 +0x1.2db954p-1, +0x1.720de4p-25, -0x1.d39b04p-49, -0x1.f11176p-76 , # /* 142: log(1846/1024.) */
	 +0x1.2a1b08p-1, -0x1.562494p-25, +0x1.a7863cp-49, +0x1.85dd64p-73 , # /* 143: log(1833/1024.) */
	 +0x1.267620p-1, +0x1.3430e0p-29, -0x1.96a958p-56, +0x1.f8e636p-82 , # /* 144: log(1820/1024.) */
	 +0x1.23130cp-1, +0x1.7bebf4p-25, +0x1.416f1cp-52, -0x1.78dd36p-77 , # /* 145: log(1808/1024.) */
	 +0x1.1faa34p-1, +0x1.70e128p-26, +0x1.81817cp-50, -0x1.c2179cp-76 , # /* 146: log(1796/1024.) */
	 +0x1.1bf204p-1, +0x1.3a9620p-28, +0x1.2f94c0p-52, +0x1.9096c0p-76 , # /* 147: log(1783/1024.) */
	 +0x1.187ce4p-1, -0x1.077870p-27, +0x1.655a80p-51, +0x1.eaafd6p-78 , # /* 148: log(1771/1024.) */
	 +0x1.1501c0p-1, -0x1.406cacp-25, -0x1.e72290p-49, +0x1.5dd800p-73 , # /* 149: log(1759/1024.) */
	 +0x1.11cb80p-1, +0x1.787cd0p-25, -0x1.efdc78p-51, -0x1.5380cep-77 , # /* 150: log(1748/1024.) */
	 +0x1.0e4498p-1, +0x1.747324p-27, -0x1.024548p-51, +0x1.77a5a6p-75 , # /* 151: log(1736/1024.) */
	 +0x1.0b036cp-1, +0x1.690c74p-25, +0x1.5d0cc4p-50, -0x1.c0e23cp-76 , # /* 152: log(1725/1024.) */
	 +0x1.077070p-1, -0x1.a769bcp-27, +0x1.452234p-52, +0x1.6ba668p-76 , # /* 153: log(1713/1024.) */
	 +0x1.04240cp-1, -0x1.a686acp-27, -0x1.ef46b0p-52, -0x1.5ce10cp-76 , # /* 154: log(1702/1024.) */
	 +0x1.00d22cp-1, +0x1.fc0e10p-25, +0x1.6ee034p-50, -0x1.19a2ccp-74 , # /* 155: log(1691/1024.) */
	 +0x1.faf588p-2, +0x1.ef1e64p-27, -0x1.26504cp-54, -0x1.b15792p-82 , # /* 156: log(1680/1024.) */
	 +0x1.f4d87cp-2, +0x1.d7b980p-26, -0x1.a114d8p-50, +0x1.9758c6p-75 , # /* 157: log(1670/1024.) */
	 +0x1.ee1414p-2, +0x1.2ec060p-26, +0x1.dc00fcp-52, +0x1.f8833cp-76 , # /* 158: log(1659/1024.) */
	 +0x1.e7e32cp-2, -0x1.ac796cp-27, -0x1.a68818p-54, +0x1.235d02p-78 , # /* 159: log(1649/1024.) */
	 +0x1.e108a0p-2, -0x1.768ba4p-28, -0x1.f050a8p-52, +0x1.00d632p-82 , # /* 160: log(1638/1024.) */
	 +0x1.dac354p-2, -0x1.d3a6acp-30, +0x1.18734cp-57, -0x1.f97902p-83 , # /* 161: log(1628/1024.) */
	 +0x1.d47424p-2, +0x1.7dbbacp-31, -0x1.d5ada4p-56, +0x1.56fcaap-81 , # /* 162: log(1618/1024.) */
	 +0x1.ce1af0p-2, +0x1.70be7cp-27, +0x1.6f6fa4p-51, +0x1.7955a2p-75 , # /* 163: log(1608/1024.) */
	 +0x1.c7b798p-2, +0x1.ec36ecp-26, -0x1.07e294p-50, -0x1.ca183cp-75 , # /* 164: log(1598/1024.) */
	 +0x1.c1ef04p-2, +0x1.c1dfd4p-26, +0x1.888eecp-50, -0x1.fd6b86p-75 , # /* 165: log(1589/1024.) */
	 +0x1.bb7810p-2, +0x1.478bfcp-26, +0x1.245b8cp-50, +0x1.ea9d52p-74 , # /* 166: log(1579/1024.) */
	 +0x1.b59da0p-2, -0x1.882b08p-27, +0x1.31573cp-53, -0x1.8c249ap-77 , # /* 167: log(1570/1024.) */
	 +0x1.af1294p-2, -0x1.b710f4p-27, +0x1.622670p-51, +0x1.128578p-76 , # /* 168: log(1560/1024.) */
	 +0x1.a925d4p-2, -0x1.0ae750p-27, +0x1.574ed4p-51, +0x1.084996p-75 , # /* 169: log(1551/1024.) */
	 +0x1.a33040p-2, +0x1.027d30p-29, +0x1.b9a550p-53, -0x1.b2e38ap-78 , # /* 170: log(1542/1024.) */
	 +0x1.9d31c0p-2, -0x1.5ec12cp-26, -0x1.5245e0p-52, +0x1.2522d0p-79 , # /* 171: log(1533/1024.) */
	 +0x1.972a34p-2, +0x1.135158p-30, +0x1.a5c09cp-56, +0x1.24b70ep-80 , # /* 172: log(1524/1024.) */
	 +0x1.911984p-2, +0x1.0995d4p-26, +0x1.3bfb5cp-50, +0x1.2c9dd6p-75 , # /* 173: log(1515/1024.) */
	 +0x1.8bad98p-2, -0x1.1d6144p-29, +0x1.5b9208p-53, +0x1.1ec158p-77 , # /* 174: log(1507/1024.) */
	 +0x1.858b58p-2, -0x1.1b4678p-27, +0x1.56cab4p-53, -0x1.2fdc0cp-78 , # /* 175: log(1498/1024.) */
	 +0x1.7f5fa0p-2, +0x1.3aaf48p-27, +0x1.461964p-51, +0x1.4ae476p-75 , # /* 176: log(1489/1024.) */
	 +0x1.79db68p-2, -0x1.7e5054p-26, +0x1.673750p-51, -0x1.a11f7ap-76 , # /* 177: log(1481/1024.) */
	 +0x1.744f88p-2, -0x1.cc0e18p-26, -0x1.1e9d18p-50, -0x1.6c06bcp-78 , # /* 178: log(1473/1024.) */
	 +0x1.6e08ecp-2, -0x1.5d45e0p-26, -0x1.c73ec8p-50, +0x1.318d72p-74 , # /* 179: log(1464/1024.) */
	 +0x1.686c80p-2, +0x1.e9b14cp-26, -0x1.13bbd4p-50, -0x1.efeb1cp-78 , # /* 180: log(1456/1024.) */
	 +0x1.62c830p-2, -0x1.a8c70cp-27, -0x1.5a1214p-51, -0x1.bab3fcp-79 , # /* 181: log(1448/1024.) */
	 +0x1.5d1bdcp-2, -0x1.4fec6cp-31, +0x1.423638p-56, +0x1.ee3feep-83 , # /* 182: log(1440/1024.) */
	 +0x1.576770p-2, +0x1.7455a8p-26, -0x1.3ab654p-50, -0x1.26be4cp-75 , # /* 183: log(1432/1024.) */
	 +0x1.5262e0p-2, -0x1.146778p-26, -0x1.b9f708p-52, -0x1.294018p-77 , # /* 184: log(1425/1024.) */
	 +0x1.4c9f08p-2, +0x1.e152c4p-26, -0x1.dde710p-53, +0x1.fd2208p-77 , # /* 185: log(1417/1024.) */
	 +0x1.46d2d8p-2, +0x1.c28058p-26, -0x1.936284p-50, +0x1.9fdd68p-74 , # /* 186: log(1409/1024.) */
	 +0x1.41b940p-2, +0x1.cce0c0p-26, -0x1.1a4050p-50, +0x1.bc0376p-76 , # /* 187: log(1402/1024.) */
	 +0x1.3bdd24p-2, +0x1.d6296cp-27, +0x1.425b48p-51, -0x1.cddb2cp-77 , # /* 188: log(1394/1024.) */
	 +0x1.36b578p-2, -0x1.287ddcp-27, -0x1.2d0f4cp-51, +0x1.38447ep-75 , # /* 189: log(1387/1024.) */
	 +0x1.31871cp-2, +0x1.2a8830p-27, +0x1.3eae54p-52, -0x1.898136p-77 , # /* 190: log(1380/1024.) */
	 +0x1.2b9304p-2, -0x1.51d8b8p-28, +0x1.27694cp-52, -0x1.fd852ap-76 , # /* 191: log(1372/1024.) */
	 +0x1.265620p-2, -0x1.d98f3cp-27, +0x1.a44338p-51, -0x1.56e85ep-78 , # /* 192: log(1365/1024.) */
	 +0x1.211254p-2, +0x1.986160p-26, +0x1.73c5d0p-51, +0x1.4a861ep-75 , # /* 193: log(1358/1024.) */
	 +0x1.1bc794p-2, +0x1.fa3918p-27, +0x1.879c5cp-51, +0x1.16107cp-78 , # /* 194: log(1351/1024.) */
	 +0x1.1675ccp-2, -0x1.4545a0p-26, +0x1.c07398p-51, +0x1.f55c42p-76 , # /* 195: log(1344/1024.) */
	 +0x1.111ce4p-2, +0x1.f72670p-37, -0x1.b84b5cp-61, +0x1.a4a4dcp-85 , # /* 196: log(1337/1024.) */
	 +0x1.0c81d4p-2, +0x1.0c150cp-27, +0x1.218600p-51, -0x1.d17312p-76 , # /* 197: log(1331/1024.) */
	 +0x1.071b84p-2, +0x1.fcd590p-26, +0x1.a3a2e0p-51, +0x1.fe5ef8p-76 , # /* 198: log(1324/1024.) */
	 +0x1.01ade4p-2, -0x1.bb1844p-28, +0x1.db3cccp-52, +0x1.1f56fcp-77 , # /* 199: log(1317/1024.) */
	 +0x1.fa01c4p-3, -0x1.12a0d0p-29, -0x1.f71fb0p-54, +0x1.e287a4p-78 , # /* 200: log(1311/1024.) */
	 +0x1.ef0adcp-3, +0x1.7b8b28p-28, -0x1.35bce4p-52, -0x1.abc8f8p-79 , # /* 201: log(1304/1024.) */
	 +0x1.e598ecp-3, +0x1.5a87e4p-27, -0x1.134bd0p-51, +0x1.c2cebep-76 , # /* 202: log(1298/1024.) */
	 +0x1.da85d8p-3, -0x1.df31b0p-27, +0x1.94c16cp-57, +0x1.8fd7eap-82 , # /* 203: log(1291/1024.) */
	 +0x1.d0fb80p-3, -0x1.bb5434p-28, -0x1.ea5640p-52, -0x1.8ceca4p-77 , # /* 204: log(1285/1024.) */
	 +0x1.c765b8p-3, +0x1.e4d68cp-27, +0x1.5b59b4p-51, +0x1.76f6c4p-76 , # /* 205: log(1279/1024.) */
	 +0x1.bdc46cp-3, -0x1.1cbb50p-27, +0x1.2da010p-51, +0x1.eb282cp-75 , # /* 206: log(1273/1024.) */
	 +0x1.b27980p-3, -0x1.1b9ce0p-27, +0x1.7756f8p-52, +0x1.2ff572p-76 , # /* 207: log(1266/1024.) */
	 +0x1.a8bed0p-3, -0x1.bbe874p-30, +0x1.85cf20p-56, +0x1.b9cf18p-80 , # /* 208: log(1260/1024.) */
	 +0x1.9ef83cp-3, +0x1.2769a4p-27, -0x1.85bda0p-52, +0x1.8c8018p-79 , # /* 209: log(1254/1024.) */
	 +0x1.9525a8p-3, +0x1.cf456cp-27, -0x1.7137d8p-52, -0x1.f158e8p-76 , # /* 210: log(1248/1024.) */
	 +0x1.8b46f8p-3, +0x1.11b12cp-30, +0x1.9f2104p-54, -0x1.22836ep-78 , # /* 211: log(1242/1024.) */
	 +0x1.83040cp-3, +0x1.2379e4p-28, +0x1.b71c70p-52, -0x1.990cdep-76 , # /* 212: log(1237/1024.) */
	 +0x1.790ed4p-3, +0x1.dc4c68p-28, -0x1.910ac8p-52, +0x1.dd1bd6p-76 , # /* 213: log(1231/1024.) */
	 +0x1.6f0d28p-3, +0x1.5cad68p-28, +0x1.737c94p-52, -0x1.9184bap-77 , # /* 214: log(1225/1024.) */
	 +0x1.64fee8p-3, +0x1.04bf88p-28, +0x1.6fca28p-52, +0x1.8884a8p-76 , # /* 215: log(1219/1024.) */
	 +0x1.5c9400p-3, +0x1.d65cb0p-29, -0x1.b2919cp-53, +0x1.b99bcep-77 , # /* 216: log(1214/1024.) */
	 +0x1.526e60p-3, -0x1.c5e4bcp-27, -0x1.0ba380p-52, +0x1.d6e3ccp-79 , # /* 217: log(1208/1024.) */
	 +0x1.483bccp-3, +0x1.9cdc7cp-28, -0x1.5ad8dcp-54, -0x1.392d3cp-83 , # /* 218: log(1202/1024.) */
	 +0x1.3fb25cp-3, -0x1.a6ad74p-27, +0x1.5be6b4p-52, -0x1.4e0114p-77 , # /* 219: log(1197/1024.) */
	 +0x1.371fc4p-3, -0x1.fe1708p-27, -0x1.78864cp-52, -0x1.27543ap-76 , # /* 220: log(1192/1024.) */
	 +0x1.2cca10p-3, -0x1.4141b4p-28, -0x1.ef191cp-52, +0x1.00ee08p-76 , # /* 221: log(1186/1024.) */
	 +0x1.242310p-3, +0x1.3ba510p-27, -0x1.d003c8p-51, +0x1.162640p-76 , # /* 222: log(1181/1024.) */
	 +0x1.1b72acp-3, +0x1.52f67cp-27, -0x1.fd6fa0p-51, +0x1.1a3966p-77 , # /* 223: log(1176/1024.) */
	 +0x1.10f8e4p-3, +0x1.129cd8p-30, +0x1.31ef30p-55, +0x1.a73e38p-79 , # /* 224: log(1170/1024.) */
	 +0x1.08338cp-3, -0x1.005d7cp-27, -0x1.661a9cp-51, +0x1.1f138ap-79 , # /* 225: log(1165/1024.) */
	 +0x1.fec914p-4, -0x1.c482a8p-29, -0x1.55746cp-54, +0x1.99f932p-80 , # /* 226: log(1160/1024.) */
	 +0x1.ed1794p-4, +0x1.d06f00p-29, +0x1.75e45cp-53, -0x1.d0483ep-78 , # /* 227: log(1155/1024.) */
	 +0x1.db5270p-4, +0x1.87d928p-32, -0x1.0f52a4p-57, +0x1.81f4a6p-84 , # /* 228: log(1150/1024.) */
	 +0x1.c97978p-4, +0x1.af1d24p-29, -0x1.0977d0p-60, -0x1.8839d0p-84 , # /* 229: log(1145/1024.) */
	 +0x1.b78c84p-4, -0x1.44f124p-28, -0x1.ef7bc4p-52, +0x1.9e0650p-78 , # /* 230: log(1140/1024.) */
	 +0x1.a58b60p-4, +0x1.856464p-29, +0x1.c651d0p-55, +0x1.b06b0cp-79 , # /* 231: log(1135/1024.) */
	 +0x1.9375e4p-4, +0x1.5595ecp-28, +0x1.dc3738p-52, +0x1.86c89ap-81 , # /* 232: log(1130/1024.) */
	 +0x1.814be4p-4, -0x1.c073fcp-28, -0x1.371f88p-53, -0x1.5f4080p-77 , # /* 233: log(1125/1024.) */
	 +0x1.6f0d28p-4, +0x1.5cad68p-29, +0x1.737c94p-53, -0x1.9184bap-78 , # /* 234: log(1120/1024.) */
	 +0x1.60658cp-4, -0x1.6c8af4p-28, +0x1.d8ef74p-55, +0x1.c4f792p-80 , # /* 235: log(1116/1024.) */
	 +0x1.4e0110p-4, +0x1.146b5cp-29, +0x1.73f7ccp-54, -0x1.d28db8p-79 , # /* 236: log(1111/1024.) */
	 +0x1.3b8758p-4, +0x1.8b1b70p-28, -0x1.20aca4p-52, -0x1.651894p-76 , # /* 237: log(1106/1024.) */
	 +0x1.28f834p-4, +0x1.43b6a4p-30, -0x1.452af8p-55, +0x1.976892p-80 , # /* 238: log(1101/1024.) */
	 +0x1.1a0fbcp-4, -0x1.e4075cp-28, +0x1.1fe618p-52, +0x1.9d6dc2p-77 , # /* 239: log(1097/1024.) */
	 +0x1.075984p-4, -0x1.4ce370p-29, -0x1.d9fc98p-53, +0x1.4ccf12p-77 , # /* 240: log(1092/1024.) */
	 +0x1.f0a30cp-5, +0x1.162a68p-37, -0x1.e83368p-61, -0x1.d222a6p-86 , # /* 241: log(1088/1024.) */
	 +0x1.cae730p-5, -0x1.1a8f7cp-31, -0x1.5f9014p-55, +0x1.2720c0p-79 , # /* 242: log(1083/1024.) */
	 +0x1.ac9724p-5, -0x1.e8ee08p-29, +0x1.a7de04p-54, -0x1.9bba74p-78 , # /* 243: log(1079/1024.) */
	 +0x1.868a84p-5, -0x1.ef8128p-30, +0x1.dc5eccp-54, -0x1.58d250p-79 , # /* 244: log(1074/1024.) */
	 +0x1.67f950p-5, -0x1.ed684cp-30, -0x1.f060c0p-55, -0x1.b1294cp-80 , # /* 245: log(1070/1024.) */
	 +0x1.494accp-5, +0x1.a6c890p-32, -0x1.c3ad48p-56, -0x1.6dc66cp-84 , # /* 246: log(1066/1024.) */
	 +0x1.22c71cp-5, -0x1.8abe2cp-32, -0x1.7e7078p-56, -0x1.ddc3dcp-86 , # /* 247: log(1061/1024.) */
	 +0x1.03d5d8p-5, +0x1.79cfbcp-31, -0x1.da7c4cp-58, +0x1.4e7582p-83 , # /* 248: log(1057/1024.) */
	 +0x1.c98d18p-6, +0x1.a01904p-31, -0x1.854164p-55, +0x1.883c36p-79 , # /* 249: log(1053/1024.) */
	 +0x1.8b31fcp-6, -0x1.356500p-30, +0x1.c3ab48p-55, +0x1.b69bdap-80 , # /* 250: log(1049/1024.) */
	 +0x1.3cea44p-6, +0x1.a352bcp-33, -0x1.8865acp-57, -0x1.48159cp-81 , # /* 251: log(1044/1024.) */
	 +0x1.fc0a8cp-7, -0x1.e07f84p-32, +0x1.e7cf6cp-58, +0x1.3a69c0p-82 , # /* 252: log(1040/1024.) */
	 +0x1.7dc474p-7, +0x1.f810a8p-31, -0x1.245b5cp-56, -0x1.a1f4f8p-80 , # /* 253: log(1036/1024.) */
	 +0x1.fe02a8p-8, -0x1.4ef988p-32, +0x1.1f86ecp-57, +0x1.20723cp-81 , # /* 254: log(1032/1024.) */
	 +0x1.ff00acp-9, -0x1.d4ef44p-33, +0x1.2821acp-63, +0x1.5a6d32p-87 , # /* 255: log(1028/1024.) */
	 0, 0, 0, 0))



## instead of  bd0_scale[i][j] in C  where i and j start with 0 !
## use         logf_mat[j+1, i+1]  in R   (i+1 and j+1 start with 1)


## /*
##  * Compute x * log (x / M) + (M - x) =
##  *      = -x * log1pmx ((M - x) / x)
##  *
##  * Deliver the result back in two parts, *yh and *yl.
##  */
ebd0.1 <- function(x, M, verbose, ...) # return  c(yh, yl)
{
    stopifnot(length(x) == 1, length(M) == 1, x >= 0, M >= 0)

    yl <- yh <- 0
    if (x == M)               return(c(yh=yh, yl=yl))
    if (x == 0) { yh <- M;    return(c(yh=yh, yl=yl)) }
    if (M == 0) { yh <- +Inf; return(c(yh=yh, yl=yl)) }

    if (is.infinite(M/x)) { yh <- M; return(c(yh=yh, yl=yl)) } # as when (x == 0)
    re <- .Call(C_R_frexp, M/x) ## FIXME: handle overflow/underflow in division 'M/x' !!!
    r <- re[["r"]]
    e <- re[["e"]]

    Sb <- 10L
    S <- 2^Sb #  = 2^10 = 1024
    N <- ncol(logf_mat)-1L # = 128

    i  <- as.integer(floor ((r - 0.5) * (2 * N) + 0.5));
    ## // now,  0 <= i <= N
    f  <- floor (S / (0.5 + i / (2.0 * N)) + 0.5);
    fg  <- .Call(C_R_ldexp, f, -(e + Sb)) # // ldexp(f, E) := f * 2^E
    if(verbose) {
	cat(sprintf("ebd0(x=%g, M=%g): M/x = (r=%.15g) * 2^(e=%d); i=%d,\n  f=%g, fg=f*2^-(e+%d)=%g\n",
                    x, M, r,e, i, f, Sb, fg))
	if (fg == Inf) {
	    cat(sprintf(" --> fg = +Inf --> return( +Inf )\n"))
            return(c(yh=fg, yl=yl))
	}
	cat("     bd0_sc[0][0..3]= ("); for(j in 1:4) cat(sprintf("%g ", logf_mat[j, 0+1L])); cat(")\n")
	cat("i -> bd0_sc[i][0..3]= ("); for(j in 1:4) cat(sprintf("%g ", logf_mat[j, i+1L])); cat(")\n")
	cat(sprintf( "  small(?)  (M*fg-x)/x = (M*fg)/x - 1 = %.16g\n", (M*fg-x)/x))
    }
    else
        if (fg == Inf) { return(c(yh=fg, yl=yl)) }

	## /* We now have (M * fg / x) close to 1.  */

	## /*
	##  * We need to compute this:
	##  * (x/M)^x * exp(M-x) =
	##  * (M/x)^-x * exp(M-x) =
	##  * (M*fg/x)^-x * (fg)^x * exp(M-x) =
	##  * (M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)
	##  *
	##  * In log terms:
	##  * log((x/M)^x * exp(M-x)) =
	##  * log((M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)) =
	##  * log((M*fg/x)^-x * exp(M*fg-x)) + x*log(fg) + (M-M*fg) =
	##  * -x*log1pmx((M*fg-x)/x) + x*log(fg) + M - M*fg =
	##  *
	##  * Note, that fg has at most 10 bits.  If M and x are suitably
	##  * "nice" -- such as being integers or half-integers -- then
	##  * we can compute M*fg as well as x * bd0_scale[.][.] without
	##  * rounding errors.
	##  */

    ADD1 <- function(d) {
        d1 <- floor (d + 0.5)
	d2 <- d - d1
        yh <<- yh+d1
        yl <<- yl+d2
    }

    if(verbose) {
        log1.. <- log1pmx((M * fg - x) / x, ...)
        d <- -x * log1..
	cat(sprintf(" 1a. before adding  -x * log1pmx(.) = -x * %g = %g\n", log1.., d))
        ADD1(d)
	cat(sprintf(" 1. after A.(-x*l..):       yl,yh = (%13g, %13g); yl+yh= %g\n",
                    yl, yh, (yl)+(yh)))
        if(fg == 1) {
            cat("___ fg = 1 ___ skipping further steps\n")
            return(c(yh=yh, yl=yl))
        }
        ##  else -- fg != 1
        cat(sprintf(" 2:  A(x*b[i,j]) and A(-x*b[0,j]), j=1:4:\n"))
        for (j in 1:4) {
            ADD1(  x * logf_mat[j, i+1L]);  # /* handles  x*log(fg*2^e) */
            cat(sprintf(" j=%d: (%13g, %13g);", j, yl, yh))
            ADD1((-x * logf_mat[j, 0+1L]) * e);  # /* handles  x*log(1/ 2^e) -- *order* important
            cat(sprintf(" (%13g, %13g); yl+yh= %g\n", yl, yh, (yl)+(yh)))
            if(!is.finite(yh)) {
                cat(" non-finite yh --> return((yh=Inf, yl=0))\n")
                return(c(yh=Inf, yl=0))
            }
        }
    } else { ## not verbose
        ADD1(-x * log1pmx ((M * fg - x) / x, ...))
        if(fg == 1) return(c(yh=yh, yl=yl))

	for (j in 1:4) {
	    ADD1(  x * logf_mat[j, i+1L])       # /* handles  x*log(fg*2^e) */
	    ADD1((-x * logf_mat[j, 0+1L]) * e)  # /* handles  x*log(1/ 2^e)  -- *order* important
            if(!is.finite(yh)) return(c(yh=Inf, yl=0))
	}
    }
    ## if(!is.finite(yh)) return(c(yh=Inf, yl=0))

    ADD1(M);
    if(verbose) cat(sprintf(" 3. after ADD1(M):            yl,yh = (%13g, %13g); yl+yh= %g\n",
                            yl, yh, (yl)+(yh)))
    ADD1(-M * fg);
    if(verbose) cat(sprintf(" 4. after ADD1(- M*fg):       yl,yh = (%13g, %13g); yl+yh= %g\n\n",
                            yl, yh, (yl)+(yh)))
    c(yh, yl) # was c(yh=yh, yl=yl)
} ## end{ ebd0.1() }

##' The vectorized version we export (and document)
ebd0 <- function(x, M, verbose = getOption("verbose"), ...) {
    ##                     n list()s of length 2 <<---  vvvvvvvvvvvvvvvv {for when length(x) > max.int}
    r <- Vectorize(ebd0.1, vectorize.args = c("x", "M"), SIMPLIFY = FALSE)(x, M=M, verbose=verbose, ...)
    ##             ^^^^^^
    n..1 <- numeric(1L)
    data.frame(yh = vapply(r, `[`, n..1, 1L),
               yl = vapply(r, `[`, n..1, 2L))
}

## #undef ADD1

## // We keep "old" bd0() : it should be faster and often sufficient
## #if 0
## // Mathematically (but not numerically) the same as  bd0() :
## double attribute_hidden
## bd0(double x, double M)
## {
## 	double yh, yl;
## 	ebd0 (x, M, &yh, &yl);
## 	return yh + yl;
## }
## #endif




## ~/R/D/r-devel/R/src/nmath/stirlerr.c     --> stirlerr()
## C Code :

##  AUTHOR
##    Catherine Loader, catherine@research.bell-labs.com.
##    October 23, 2000.
##
##  Merge in to R:
 ##	Copyright (C) 2000, The R Core Team

##  DESCRIPTION
##
##    Computes the log of the error term in Stirling's formula.
##      For n > 15, uses the series 1/12n - 1/360n^3 + ...
##      For n <=15, integers or half-integers, uses stored values.
##      For other n < 15, uses lgamma directly (don't use this to
##        write lgamma!)
##
## Merge in to R:
## Copyright (C) 2000, The R Core Team
## R has lgammafn, and lgamma is not part of ISO C
##


##  stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
##              = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
##              = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
##
##  see also lgammacor() in ./lgammacor.c  which computes almost the same!

" 3.14159'26535'89793'23846'26433'83279'50288'41971'69399'37510'58209'74944'59230'78164'06286'20899'86280'34825'34211'70679
  0     5     0     5     0     5     0     5     0     5     0     5     0     5     0     5     0     5     0     5     0
  0           1           2           3           4           5           6           7           8           9          10
"

##' [D]irect formula for stirlerr(), notably adapted to mpfr-numbers
stirlerr_simpl <- function(n, version = c("R3", "lgamma1p", "MM2", "n0"), minPrec = 128L) {
    if(notNum <- !is.numeric(n)) {
        precB <- if(isM <- inherits(n, "mpfr"))
                     max(minPrec, Rmpfr::.getPrec(n))
                 else if(isZQ <- inherits(n, "bigz") || inherits(n, "bigq"))
                     max(minPrec, Rmpfr::getPrec(n))
                 else
                     minPrec
        pi <- Rmpfr::Const("pi", precB)
        M_LN_2PI      <- log(2*pi)
        M_LN_SQRT_2PI <- M_LN_2PI/2
    } else { ## numeric
        M_LN_2PI <- 1.837877066409345483 # ~ log(2*Const("pi",60)); very slightly more accurate than log(2*pi):
        ## 1 - log(2*pi)/M_LN_2PI # 1.110223e-16 = 2^-53  is "fine" (but not zero)
        ## R's src/nmath/ had log(sqrt(2*Pi)) == log(2*Pi)/2
        M_LN_SQRT_2PI <- 0.91893853320467274153 # and again 1 - (log(2*pi)/2)/M_LN_SQRT_2PI |--> 1.11022e-16
    }
    if(notNum && !isM) {
        n <- if(isZQ)
                 Rmpfr::mpfr(n, precB)
             else ## the object-author "should" provide a method:
                 as(n, "mpfr")
    }
    ## direct formula (suffering from cancellation for largish n)
    switch(match.arg(version),
           "R3"     = lgamma(n + 1) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI, # R <= 4.3.*'s nmath/stirlerr.c
           "lgamma1p" = lgamma1p(n) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI,
           "MM2"    = { ## a few bits better for  n in [1, 5] :
                        l.n <- log(n)
                        lgamma(n) + n*(1 - l.n) + (l.n - M_LN_2PI)/2
           },
           "n0" = -(log(n) + M_LN_2PI)/2,
           stop("invalid 'version': ", match.arg(version)))
}

## Interfacing  C dpq_stirlerr() in  ../src/stirlerr.c
##                ~~~~~~~~~~~~~~
## called from  gammafn_ver() <-- R_gamma_ver()   in ../src/gamma-variants.c <<<<<<<<
##  *and* from our R level           gammaVer()
## now (2024-1-19) also from here
stirlerrC <- function(n, version = c("R3", "R4..1", "R4.4_0"))
    .Call(C_R_dpq_stirlerr, n, match(match.arg(version), eval(formals()$version)))


## R interface to C level gamma() variants   <==> ../src/gamma-variants.c
gammaVer <- function(x, version, stirlerrV = c("R3", "R4..1", "R4.4_0"),
                     traceLev = 0L) {     ##    --    ----    ----  <--> ../src/DPQpkg.h stirlerr_version_t
    ## TODO: iver <- pmatch(version, c("...", "..",)) or via match.arg()
    stopifnot(is.integer(iver <- as.integer(version)), 1 <= iver, iver <= 5,
              is.integer(traceLev <- as.integer(traceLev)), 0 <= traceLev)
    iv.stirlerr <- match(match.arg(stirlerrV), eval(formals()$stirlerrV)) # 1, 2, ..
    ## both iv.stirlerr and  iver = version needs to be integer: currently in 1..5
    .Call(C_R_gamma_ver, x, iver, traceLev, iv.stirlerr)
}


##' stirlerr() now *vectorized* in  n
stirlerr <- function(n, scheme = c("R3", "R4.4_0"),
                     cutoffs = switch(scheme
                                    , R3   = c(15, 35, 80, 500)
                                    , R4.4_0 = c(5.25, rep(6.5, 4), 7.1, 7.6, 8.25, 8.8, 9.5, 11,
                                                 14, 19,   25, 36, 81, 200, 3700, 17.4e6)
                                      )
                    , use.halves = missing(cutoffs)
                     ## had  lgamma1p = lgamma1pC # was "hard wired"  \(x) lgamma(x+1)
                     ## rather -- always the same as formals(stirlerr_simpl)$version : <==> TODO: sync!
                    , direct.ver = c("R3", "lgamma1p", "MM2", "n0")
                    , order = NA
                    , verbose = FALSE
                     )
{
    useBig <- (!is.numeric(n) &&
               (inherits(n, "mpfr") || inherits(n, "bigz") || inherits(n, "bigq")))
    if(useBig) {
        if(verbose) message(sprintf("stirlerr(n): As 'n' is \"%s\", using \"mpfr\" & stirlerrM():", class(n)))
        ##  DPQmpfr >= 0.3-1 on CRAN  where  DPQmpfr::stirlerrM() exists :
        if(requireNamespace("DPQmpfr") &&
           is.function(stirlFn <- get0("stirlerrM", asNamespace("DPQmpfr"), inherits=FALSE))) {
            ## this will warn in R CMD check before DPQmpfr is updated there!
            return(stirlFn(n)) # , precB = <n>   would be possible
        } else
            stop("Need CRAN package 'DPQmpfr' (>= 0.3-1) with its stirlerrM()")
    }
    else {
        ## Derivation: require(Rmpfr); require(gmp); n <- 2*(1:20)
        ##             Sn <- abs(BernoulliQ(n)) / (n*(n-1)) # == |B_n| / (n(n-1))  for n = 2,4,6,...
        ##             noquote(cbind(k = n/2-1, format(mpfr(Sn, 80)), format(Sn)))
        ## -------------> ../Misc/stirlerr-trms.R <-----
        ##               ~~~~~~~~~~~~~~~~~~~~~~~
        S0 <- 0.083333333333333333333       ## 1/12 */
        S1 <- 0.00277777777777777777778     ## 1/360 */
        S2 <- 0.00079365079365079365079365075 ## 1/1260
        S3 <- 0.00059523809523809523809523806 ## 1/1680
        S4 <- 0.00084175084175084175084175104 ## 1/1188
        S5 <- 0.0019175269175269175269175262  ## 691/360360
        S6 <- 0.0064102564102564102564102561  ## 1/156
        S7 <- 0.029550653594771241830065352   ## 3617/122400
        S8 <- 0.17964437236883057316493850    ## 43867/244188
        S9 <- 1.3924322169059011164274315     ## 174611/125400
        S10<- 13.402864044168391994478957     ## 77683/5796
        S11<- 156.84828462600201730636509     ## 236364091/1506960
        S12<- 2193.1033333333333333333333     ## 657931/300
        S13<- 36108.771253724989357173269     ## 3392780147/93960
        S14<- 691472.26885131306710839498     ## 1723168255201/2492028
        S15<- 15238221.539407416192283370     ## 7709321041217/505920
        S16<- 382900751.39141414141414141     ## 151628697551/396
        S17<- 10882266035.784391089015145     ## 26315271553053477373/2418179400
        S18<- 347320283765.00225225225225     ## 154210205991661/444
        S19<- 12369602142269.274454251718     ## 261082718496449122051/21106800
        oMax <- 20L # number of S<j> above

        stopifnot(length(order) == 1L,
                  is.na(order) || (order == (oi <- as.integer(order)) && 1L <= (order <- oi) && oi <= oMax))

### NB: When adding here:  **sync** with
### --  1. ../src/stirlerr.c   dpq_stirlerr()             <<<<<<<<<<<<<<<<<
### --  2. ../../DPQmpfr/R/dpq-mpfrized.R   stirlerrSer() <<<<<<<<<<<<<<<<<

        nn  <- n^2 # = n*n  (but that might integer-overflow)

        if(!is.na(order)) { # and use it -- rather than cutoffs etc
            return(switch(order, # already know it's in 1..iC
            S0/n,				# 1
            (S0-S1/nn)/n,			# 2
            (S0-(S1-S2/nn)/nn)/n,		# 3
            (S0-(S1-(S2-S3/nn)/nn)/nn)/n,	# 4
            (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n,	# 5
            (S0-(S1-(S2-(S3-(S4-S5/nn)/nn)/nn)/nn)/nn)/n,	# 6
            (S0-(S1-(S2-(S3-(S4-(S5-S6/nn)/nn)/nn)/nn)/nn)/nn)/n,	# 7
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-S7/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n,	# 8
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-S8/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n,	# 9
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-S9/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n,	# 10
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n,	# 11
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-S11/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 12
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-S12/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 13
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-S13/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 14
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-S14/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 15
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-S15/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 16
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-S16/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 17
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-S17/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 18
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-(S17-S18/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 19
            (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-(S17-(S18-S19/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n, # 20
            stop("invalid 'order': ", order)))
        }
        ## --- else use cutoffs[]  --------------------------------------------------------------------
        hasC <- !missing(cutoffs)
        scheme <- match.arg(scheme)
        stopifnot(is.numeric(cutoffs), 4 <= (nC <- length(cutoffs)), nC <= oMax, # -5+i.c >= 1
                  cutoffs[1] <= 15, # sferr_halves[] only for n <= 15
                  !is.unsorted(cutoffs))

        i.c <- nC - if(scheme == "R3") 4L else 5L # in {0,1,2,..., 6} or now {0,1,...,8}
        if(verbose) cat(sprintf("stirlerr(n, %s) :",
                                if(hasC) paste("cutoffs =", paste(formatC(cutoffs), collapse=","))
                                else paste0('scheme = "', scheme, '"')))
        r <- rep_len(NA_real_, length(n))
        if(use.halves && any(s15 <- n <= 15) &&
           any(hlf <- s15 & n+n == (n2 <- trunc(n+n))))
        {
            if(verbose) { cat(" use.halves (n <= 15):  n ="); nh <- n[hlf]
                if(length(nh) <= 10) cat(nh, "\n") else cat(head(nh), "...", tail(nh, 3),"\n") }
            r[hlf] <- sferr_halves[n2[hlf] + 1L]
        } else
            hlf <- FALSE
        S <- !hlf
        if (any(S <- S & (sml <- n <= cutoffs[1]))) {
            if(verbose) cat(" case I (n <= ",format(cutoffs[1]),"), ", sep="")
            n. <- n[S]
            if(verbose) { cat(" using direct formula for n="); str(n.) }
            r[S] <- stirlerr_simpl(n., version = direct.ver)
        }
        if (any(!sml)) { # has n > cutoffs[1]
            if(verbose) {
                cat(" case II (n > ",format(cutoffs[1]),"), ",nC, " cutoffs: (",
                    paste(cutoffs, collapse=", "),"):  n in cutoff intervals:")
                print(table(cut(n[n > cutoffs[1]], unique(c(cutoffs,Inf)))))
            }

            if(scheme != "R3" && length(i <- which(n > cutoffs[5+i.c])))        # k = 1 term:  n > c_[nC] = 1.7e7
                r[i] <- S0/n[i]
            if((scheme == "R3" && length(i <- which(n > cutoffs[4+i.c]))) ||
               (scheme != "R3" && length(i <- which(cutoffs[5+i.c] >= n & n > cutoffs[4+i.c])))) # k = 2 terms
                r[i] <- (S0-S1/nn[i])/n[i]
            if(length(i <- which(cutoffs[4+i.c] >= n & n > cutoffs[3+i.c])))	# k = 3 terms
                r[i] <- (S0-(S1-S2/nn[i])/nn[i])/n[i]
            if(length(i <- which( cutoffs[3+i.c] >= n & n > cutoffs[2+i.c])))	# k = 4 terms
                r[i] <- (S0-(S1-(S2-S3/nn[i])/nn[i])/nn[i])/n[i]
            if(length(i <- which( cutoffs[2+i.c] >= n & n > cutoffs[1+i.c])))	# k = 5 terms --(was 15 < n <= 35)
                r[i] <- (S0-(S1-(S2-(S3-S4/nn[i])/nn[i])/nn[i])/nn[i])/n[i]	#                now 26 < n <= 55
            if(i.c >= 1 && length(i <- which( cutoffs[1+i.c] >= n & n > cutoffs[0+i.c]))) # k = 6  20 < n <= 26 :
                r[i] <- (S0-(S1-(S2-(S3-(S4-S5/nn[i])/nn[i])/nn[i])/nn[i])/nn[i])/n[i]
            ## from here on: be careful to *NOT* destroy  r[hlf] values -- depends on cutoffs[i.c-1] (= 12.125) < 15
            rf <- !hlf # "rough": not half integers
            if(i.c >= 2 && length(i <- which(rf & cutoffs[0+i.c] >= n & n > cutoffs[-1+i.c]))) { # k = 7  12 < n <= 20 :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-S6/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 3 && length(i <- which(rf & cutoffs[-1+i.c] >= n & n > cutoffs[-2+i.c]))) { # k = 8   . < n <= 12 :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-S7/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 4 && length(i <- which(rf & cutoffs[-2+i.c] >= n & n > cutoffs[-3+i.c]))) { # k = 9   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-S8/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 5 && length(i <- which(rf & cutoffs[-3+i.c] >= n & n > cutoffs[-4+i.c]))) { # k = 10   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-S9/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 6 && length(i <- which(rf & cutoffs[-4+i.c] >= n & n > cutoffs[-5+i.c]))) { # k = 11   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 7 && length(i <- which(rf & cutoffs[-5+i.c] >= n & n > cutoffs[-6+i.c]))) { # k = 12   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-S11/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 8 && length(i <- which(rf & cutoffs[-6+i.c] >= n & n > cutoffs[-7+i.c]))) { # k = 13   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-S12/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i]
            }
            if(i.c >= 9 && length(i <- which(rf & cutoffs[-7+i.c] >= n & n > cutoffs[-8+i.c]))) { # k = 14   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-S13/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 14
            }
            if(i.c >= 10 && length(i <- which(rf & cutoffs[-8+i.c] >= n & n > cutoffs[-9+i.c]))) { # k = 15   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-S14/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 15
            }
            if(i.c >= 11 && length(i <- which(rf & cutoffs[-9+i.c] >= n & n > cutoffs[-10+i.c]))) { # k = 16   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-S15/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 16
            }
            if(i.c >= 12 && length(i <- which(rf & cutoffs[-10+i.c] >= n & n > cutoffs[-11+i.c]))) { # k = 17   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-S16/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 17
            }
            if(i.c >= 13 && length(i <- which(rf & cutoffs[-11+i.c] >= n & n > cutoffs[-12+i.c]))) { # k = 18   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-S17/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 18
            }
            if(i.c >= 14 && length(i <- which(rf & cutoffs[-12+i.c] >= n & n > cutoffs[-13+i.c]))) { # k = 19   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-(S17-S18/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 19
            }
            if(i.c >= 15 && length(i <- which(rf & cutoffs[-13+i.c] >= n & n > cutoffs[-14+i.c]))) { # k = 20   . < n <= . :
                n2 <- nn[i]
                r[i] <- (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-(S17-(S18-S19/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n2)/n[i] # 20
            }
        }
        ## return
        r
    }
}

## stirlerr(n) for n = 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
## const static double
sferr_halves <- c(
    Inf, #  n=0
    0.1534264097200273452913848,  ## 0.5 */
    0.0810614667953272582196702,  ## 1.0 */
    0.0548141210519176538961390,  ## 1.5 */
    0.0413406959554092940938221,  ## 2.0 */
    0.03316287351993628748511048, ## 2.5 */
    0.02767792568499833914878929, ## 3.0 */
    0.02374616365629749597132920, ## 3.5 */
    0.02079067210376509311152277, ## 4.0 */
    0.01848845053267318523077934, ## 4.5 */
    0.01664469118982119216319487, ## 5.0 */
    0.01513497322191737887351255, ## 5.5 */
    0.01387612882307074799874573, ## 6.0 */
    0.01281046524292022692424986, ## 6.5 */
    0.01189670994589177009505572, ## 7.0 */
    0.01110455975820691732662991, ## 7.5 */
    0.010411265261972096497478567, ## 8.0 */
    0.009799416126158803298389475, ## 8.5 */
    0.009255462182712732917728637, ## 9.0 */
    0.008768700134139385462952823, ## 9.5 */
    0.008330563433362871256469318, ## 10.0 */
    0.007934114564314020547248100, ## 10.5 */
    0.007573675487951840794972024, ## 11.0 */
    0.007244554301320383179543912, ## 11.5 */
    0.006942840107209529865664152, ## 12.0 */
    0.006665247032707682442354394, ## 12.5 */
    0.006408994188004207068439631, ## 13.0 */
    0.006171712263039457647532867, ## 13.5 */
    0.005951370112758847735624416, ## 14.0 */
    0.005746216513010115682023589, ## 14.5 */
    0.005554733551962801371038690  ## 15.0 */
)
