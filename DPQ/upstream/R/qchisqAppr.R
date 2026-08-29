## This uses the "initial approx" of AS 91 [as in qgamma()]
## which should be *better* than all others below:

## source("/u/maechler/R/MM/NUMERICS/dpq-functions/dpq-h.R")

## NOTA BENE:
## ---------  --> ./qgamma-fn.R for function  qchisqAppr.R()
##                  ~~~~~~~~~~~               ---------------

qchisqAppr <- function(p, df, lower.tail = TRUE, log.p = FALSE, tol = 5e-7)
{
  ## Purpose: Cheap fast approximation to qchisq()
  ## ----------------------------------------------------------------------
  ## Arguments: tol: tolerance with default = EPS2 of qgamma()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Mar 2004, 12:07

    ## Vectorized (in 'p' only)  :
    if(length(df  <- as.double(df)) > 1) stop("'df' must be length 1")
    if(length(tol <- as.double(tol)) > 1 || tol <= 0)
        stop("'tol' must be a number > 0")
    storage.mode(p) <- "double"
    n <- as.integer(length(p))
    .C(C_qchisq_appr_v,
       p, n, df, tol,
       as.logical(lower.tail), as.logical(log.p),
       q = numeric(n))$q
}

## qchisq( p,  df,    lower_tail,  log_p) ==
## qgamma( p, 0.5 * df, 2.0, lower_tail, log_p)

qgammaAppr <- function(p, shape, lower.tail = TRUE, log.p = FALSE, tol = 5e-7)
    0.5* qchisqAppr(p, df=2*shape, lower.tail=lower.tail, log.p=log.p, tol=tol)


if(FALSE) ##=== really rather use qchisqAppr.R() in  ./qgamma-fn.R
qchisq.appr.Kind <-
    function(p, nu, g = lgamma(nu/2), lower.tail = TRUE, log.p = FALSE,
             tol = 5e-7, maxit = 1000, verbose = getOption('verbose'))
{
    ## Purpose: Cheap fast "initial" approximation to qgamma()
    ## --- also to explore the different kinds / cutoffs..
    ## Return information about the KIND of qchisq() computation
    ## ----------------------------------------------------------------------
    ## Arguments: tol: tolerance with default = EPS2 of qgamma()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Mar 2004, 10:36

    if(length(nu) != 1 || length(g) != 1)
        stop("arguments must have length 1 !")
    if (nu <= 0) stop("need nu > 0")

    if(length(p) != 1) stop("argument 'p' must have length 1 !")
    if(is.na(p) || is.na(nu)) stop("'p' and 'nu' must not be NA!")
    if(log.p) stopifnot(p <= 0) else stopifnot(0 <= p && p <= 1)

    alpha <- 0.5 * nu ##/* = [pq]gamma() shape */
    c <- alpha-1
    force(g)

    Cat <- function(...) if(verbose > 0) cat(...)

    ##p. <- .DT_qIv(p, lower.tail, log.p) ##/* lower_tail prob (in any case) */
    p1 <- .DT_log(p, lower.tail, log.p)

    kind <- {
        if(nu < -1.24 * p1) "chi.small"
        else if(nu > 0.32) {
            ## 'ch' not known! if(ch > 2.2*nu + 6) "p1WH" else "WH"
            "WHchk"
        }
        else "nu.small"
    }

    switch(kind,
           "chi.small" = { ##/* for small chi-squared */

               ch <- exp((log(alpha) + p1 + g)/alpha + M_LN2)
               Cat(sprintf(" small chi-sq., ch0 = %g\n", ch))
           },
           "WHchk" = { ##/*  using Wilson and Hilferty estimate */

               x <- qnorm(p, 0, 1, lower.tail, log.p)
               p1 <- 2./(9*nu)
               ch <- nu* (x*sqrt(p1) + 1-p1)^3
               Cat(sprintf(" nu > .32: Wilson-Hilferty; x = %7g\n", x))

               kind <- (if(ch > 2.2*nu + 6) "p1WH" else "WH")
               if(kind == "p1WH") ##/* approximation for p tending to 1: */
                   ch <- -2*(.DT_Clog(p, lower.tail,log.p)
                             - c*log(0.5*ch)+ g)
           },
           "nu.small" = { ##/* small  nu : 1.24*(-log(p)) <= nu <= 0.32 */
               C7 <- 4.67
               C8 <- 6.66
               C9 <- 6.73
               C10 <- 13.32

               ch <- 0.4
               a <- .DT_Clog(p, lower.tail, log.p) + g + c*M_LN2
               Cat(sprintf(" nu=%g <= .32: a(p,nu) = %19.13g ", nu,a))

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

    list(kind= kind, it= if(kind=="nu.small") it, ch=ch)
}## qchisq.appr.Kind()


## Wilson-Hilferty  for Chisquare : (not always good !)
qchisqWH <- function(p,df, lower.tail=TRUE, log.p=FALSE) {
    ## FIXME: add 'ncp' non-central version (has nice W.-.H. trafo as well!
    p1 <- 2/(9*df)
    Z <- qnorm(p, mean = 1 - p1, sd = sqrt(p1),
               lower.tail=lower.tail, log.p=log.p)
    if(log.p) # log(df * Z^3), but Z above is already log(Z):
	log(df) + 3*Z
    else
	df*Z^3
}

## Kennedy & Gentle p.118 (according to Joe Newton's  "Comp.Statist",p.25f.)
##
## This is just a cheap version of "Phase I  Starting Approximation"
## used  in  R's  qgamma.c, i.e. AS 91 (JRSS, 1979) :
qchisqKG <- function(p, df, lower.tail=TRUE, log.p=FALSE)
{
    ## u <- .DT_qIv(p, lower.tail, log.p)
    ## lu <- log(u)
    lu <- .DT_log(p, lower.tail, log.p)
    if(df < -1.24 * lu) { ## small "p"
        exp((lu + log(df) + lgamma(df/2) + (df/2 - 1)*log(2))* (2/df))
    }
    else { ## Wilson-Hilferty:
        q0 <- qchisqWH(p, df, lower.tail=lower.tail, log.p=log.p)
        if(q0 <= 2.2*df + 6)
            q0
        else  {
            Iu <- .DT_CIv(p, lower.tail, log.p)## == 1 - u (numerically stably)
            ## the following gives NaN when u=1 (numerically),
            ## e.g., (-38, 1, lower=F, log=TRUE):
            -2 * log((Iu*gamma(df/2)) / ((q0/2)^(df/2 - 1)))
        }
    }
}

qgammaApprKG <- function(p, shape, lower.tail = TRUE, log.p = FALSE)
    0.5* qchisqKG(p, df=2*shape, lower.tail=lower.tail, log.p=log.p)
## in the case  ' small "p" ',
## this is
qgammaApprSmallP <- function(p, shape, lower.tail = TRUE, log.p = FALSE) {
    ##---- qgamma() approximation for small p -- particularly useful for small shape !

    ## u <- .DT_qIv(p, lower.tail, log.p)
    ## lu <- log(u)
    lu <- .DT_log(p, lower.tail, log.p)
    ## for shape << 1,   log(shape) + lgamma(shape)  has large cancellation
    ##  ==> use log(gamma(shape+1)) in stable form
    exp((lu + lgamma1p(shape))/shape)
}

## From which on is this good?
## If we look at  Abramowitz&Stegun gamma*(a,x) = x^-a * P(a,x)
## and its series  g*(a,x) = 1/gamma(a) * (1/a - 1/(a+1) * x + ...)
## the approximation P(a,x) = x^a * g*(a,x) ~= x^a/gamma(a+1)
## -- and hence  x = qgamma(p, a) ~= (p * gamma(a+1)) ^ (1/a)  ---
## is good as soon as  1/a >> 1/(a+1) * x
## <==>  x << (a+1)/a = (1 + 1/a)
## <==>  x < eps *(a+1)/a
## <==>  log(x) < log(eps) + log( (a+1)/a ) = log(eps) + log((a+1)/a)  ~  -36 - log(a)
##     where log(x) ~= log(p * gamma(a+1)) / a = (log(p) + lgamma1p(a))/a
## such that the above
## <==>  (log(p) + lgamma1p(a))/a < log(eps) + log((a+1)/a)
## <==>  log(p) + lgamma1p(a) < a*(-log(a)+ log(eps) + log1p(a))
## <==>  log(p) <  a*(-log(a)+ log(eps) + log1p(a)) - lgamma1p(a) =: bnd(a)
.qgammaApprBnd <- function(a, logEps = -36.043653389117156) { # = log(.Machine$double.eps)
    a*(logEps + log1p(a) - log(a)) - lgamma1p(a)
}

## NOTA BENE:  ./beta-fns.R
## =========     ~~~~~~~~~~ now has lgamma1p() "as in src/nmath/pgamma.c"
##  and that is **better**  {mainly: accurate up to a=1/2 !} than this
## for  its  lgamma1p() !

## MM FIXME:  Test that and write a vignette about it!

lgamma1p. <- function(a, cutoff.a = 1e-6, k = 3) {
    stopifnot(cutoff.a < 1/2) # as logic below uses mid := {a; cutoff.a <= a < 1/2 }
    ##  log(a * gamma(a)) == log(gamma(1+a)) == lgamma(1+a) --- notably for small a
    ##  -----------------                       -----------
    r <- a

### MM(FIXME): This expansion is *not* identical to the expansion of log(Gamma(u+1))
###  and that is actually simpler !!
###  In maple,         series(log(Gamma(1+x)), x, 9);
###  gives more good terms --> ~/maple/gamma-asympt.tex and
### then ~/maple/gamma-asympt.R  ==> lgamma1p_series() here in ./beta-fns.R

    ## Taylor-expansion:  Gamma(1+u) = 1 + u*(-gammaE + a_0* u + a_1* u^2 + a_2* u^3) + O(u^5)
    ## psi(1) = digamma(1) = -(Euler's) gamma = -Const("gamma",200)
    gammaE <- 0.57721566490153286060651209008240243104215933593992359880576723
    ## psi'(1) = trigamma(1) = pi^2/6
    ## a_0 = (psi'(1) + psi(1)^2)/2  = (pi^2/6 + gamma^2)/2 =
    ## require("Rmpfr");
    a0 <- 0.98905599532797255539539565150063470793918352072821409044319567
    ## a_1 = (psi''(1) + 3*psi(1)*psi'(1) + psi(1)^3)/6,
    ## now psi''(1) = - 2*zeta(3)  (psi2.1 <- -2*zeta(mpfr(3,200)))
    ## a_1 = (-2*zeta(3) - gammaE*(3*pi^2/6 + gammaE^2))/6 =
    a1 <- -0.90747907608088628901656016735627511492861144907256376094133062

    ## In Maple:    asympt(GAMMA(1/x), x, 5); or even better
    ##             series(log(Gamma(1+x)), x, 9);
    ## gives  even more terms --> ~/maple/gamma-asympt.tex or *.mw
    ##                            -----------------------------------
    ## a_2 = (1/160)*Pi^4+(1/3)*Zeta(3)*gamma+(1/24)*Pi^2*gamma^2+(1/24)*gamma^4
    ##      (1/160)*Pi^4+ gamma/24*(8*Zeta(3)+ Pi^2*gamma+ *gamma^4)
    a2 <- 0.97977257665513641646834022654634789065169603841314073131007745
    ## two terms u*(-g. + a_0*u) is "exact" once  a_1*u^2 < g. * eps
    ## i.e. u < sqrt( (g./a_1) * eps) = 3.47548562941137e-08
    ## but, as we found empirically, cutting off earlier *is* better
    ## sml <- a < 3.47548562941137e-08
    sml <- a < cutoff.a
    if(any(sml)) {
        a. <- a[sml]
        stopifnot(length(k) == 1, round(k) == k, 1 <= k, k <= 3)
	r[sml] <- log1p(switch(k,
			       a.*(-gammaE + a.*a0),			# k = 1
			       a.*(-gammaE + a.*(a0+a.*a1)),		# k = 2
			       a.*(-gammaE + a.*(a0+a.*(a1 + a.*a2)))	# k = 3
			       ))
    }
    if(any(ok <- !sml)) {
        if(any(lrg <- a >= 1/2)) # ==> 1+a  does not cancel
            r[lrg] <- lgamma(1+ a[lrg])
        if(any(mid <- ok & !lrg)) { ## cutoff.a <= a < 1/2
            a <- a[mid]
            r[mid] <- log(a*gamma(a))
        }
    }
    r
}
