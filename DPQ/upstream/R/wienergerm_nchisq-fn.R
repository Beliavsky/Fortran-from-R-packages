#### -*- mode: R ; delete-old-versions: never -*-

### Some docu in header of ../src/wienergerm_nchisq.c
###                        --------------------------

## Originally	~/R/MM/NUMERICS/dpq-functions/wienergerm_nchisq-fn.R
##		====================================================
## which goes back to  Jan 29, 2004

##' h(.) direct
h0 <- function(y) { y2 <- y * y; ((1-y)*log(1 - y) + y - 0.5*y2) / y2}

##' h(.) direct, version 1
h1 <- function(y) { y2 <- y * y; ((1-y)*log1p(- y) + y*(1- 0.5*y)) / y2}

##' h(.) direct (nterms = 0) or  nterm-Taylor approx
hnt <- function(y, nterms)
{
    if(nterms <= 0)# && u > 0 )
        ((1-y)*log1p(- y) + y*(1- y/2)) / (y*y)
    else ## use Taylor
        if(nterms <= 1) y/6*(1 + y/2)
    else { ## more terms: use Horner
        s <- 0
        for(k in nterms:1)
            s <- s * y + 1/((k+1)*(k+2)) ## s_k := s_{k+1} * y + c_k
        y * s
    }
}

##' h(.)  "Optimally implemented: Use 2- or 3-term Taylor when needed, direct otherwise
h2 <- function(y, eps = .Machine$double.eps, verbose = getOption("verbose"))
{
    if(verbose) cat("h(",formatC(y),")\n",sep='')
    r <- ay <- abs(y)
    ## ifelse(ay < (10/3 * eps)^ (1/2), y/6*(1 + y/2),
    ##        ifelse(ay < (5 * eps)^ (1/3), y*(1/6 + y*(1/12 + y/20)),
    ##               ifelse(y == 1, 1/2,
    ##                      ((1-y)*log1p(- y) + y*(1- y/2)) / (y*y))))
    y. <- y[c1 <- ay < (10/3 * eps)^ (1/2)] # (10/3 * eps)^ (1/2) == 2.720567e-08
    r[c1] <- y./6*(1 + y./2)
    c2 <- !c1 & ay < (5 * eps)^(1/3)        #     (5 * eps)^(1/3) == 1.035468e-05
    if(any(c2)) {
        y. <- y[c2]
        I_6 <- 1/(6+0*y.[1]) # (for Rmpfr, e.g., to have the same precision in I/6 = 1/6)
        r[c2] <- y.*(I_6 + y.*(I_6/2 + y./20))
    }
    r[c3 <- y == 1] <- 1/2
    y. <- y[c4 <- !(c1 | c2 | c3)]
    r[c4] <- ((1-y.)*log1p(- y.) + y.*(1- y./2)) / (y.*y.)
    r
}
h <- h2
##----- use the optimal one

### R(1-s) or rather 1/2* R(1-s) is what I need in z.s() for
### stable s -> 1 evaluation :
gnt <- function(u, nterms, times.u = FALSE)
{
    ## Purpose: g(u) := R(1-u) / 2
    ## ----------------------------------------------------------------------
    ## Arguments: u: numeric vector <= 1 (u > 1 ==> NaN)
    ##       nterms: number of terms in Taylor; nterms <= 0 := "direct formula"
    ##      times.u: if TRUE, compute u*g(u) instead of g(u)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 29 Jan 2004, 11:31

    if(nterms <= 0) { # && u > 0
        s <- 6*h(u)/u - 1
        return(if(times.u) s else s/u)
    }
    ## else: use Taylor
    if(nterms <= 1) { s <- (1 + 0.6*u)/2 ; return(if(times.u) u*s else s) }
    ## else more terms: use Horner
    s <- 0
    for(k in nterms:0)
        s <- s * u + 1/((k+3)*(k+4)) ## s_k := s_{k+1} * u + c_k
    return(6 * (if(times.u) u*s else s))
}

##' Final solution (for now)
g2 <- function(u, eps = .Machine$double.eps)
{
    ## Purpose: g(u) := R(1-u) / 2
    ## ----------------------------------------------------------------------
    ## Arguments: u: numeric vector <= 1 (u > 1 ==> NaN)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Jan 2004, 22:18

    au <- abs(u)
    ifelse(au < (2.5 * eps)^ (1/2), (1 + 0.6*u)/2,
           ifelse(au < (3.5 * eps)^ (1/3), (1 + u/5*(3 + 2*u))/2,
                  ifelse(au < (14/3 * eps)^(1/4),
                         (1 + u*(0.6 + u*(0.4 + u*2/7)))/2,
                         ## IMPROVE (not in 'ifelse'):
                         ## Compute the 'k' := #{terms 0:(k-1)}
                         ## else u large enough (?)
                          (6*h(u)/u - 1) / u )))
}



sW <- function(x, df, ncp)
{
    ## Purpose: s(x, df, ncp) as in Wienergerm approx.
    ## -------  but using Taylor expansion when needed: (x*ncp / df^2) << 1
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 26 Jan 2004, 17:58

    if (any(ncp < 0)) stop("negative 'ncp'")
    if (any(df <= 0)) stop("non-positive 'df'")

    mu2 <- ncp / df
    ff <- sqrt(1 + 4. * x * mu2 / df) ## /* === 1 + 2 * mu2 * s */
    e <- 2*mu2*x/df
    s <- ifelse(e < 8/7 * .Machine$double.eps ^ 0.25,
                x/df * (1 + e*(-1/2 + e*(1/2 - 5/8*e))), ## << Taylor
                (ff - 1)/ (2*mu2)) ## /* === s */
    list(ff = ff, s = s)
}

qs <- function(x,df,ncp, f.s = sW(x,df,ncp), eps1 = 1/2, sMax = 1e100)
{
    ## Purpose:
    ##     q(s) := (f - 2*h(1-s)) / s
    ##	        == 2*mu^2 + (1 - 2*h(1-2))/s
    ##	        == 2*(mu^2 - is*(log(s)*is + 1))  for  is = 1/(1-s)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Mar 2004, 10:05

    ff <- f.s$ff
    s  <- f.s$s
    if(!is.numeric(s) && length(s) == length(ff))
        stop("wrong 'f.s' argument")
    mu2 <- (ncp / df)

    ifelse(s < eps1, {is <- 1/(1-s) ; 2*(mu2 - is*(log(s)*is + 1))},
           ifelse(s > sMax, 2*mu2 + (1 + 2*h(1 - 1/s))/s,
                  (ff + 2*h(1- 1/s)) / s))
}

z0 <- function(x, df, ncp)
{
    ## Purpose: the first 'z'  in *both* Wienergerm approx.
    ## ----------------------------------------------------------------------
    ## Author: Martin Mächler, Date: 26 Jan 2004, 21:09

    f.s <- sW(x, df, ncp)
    ff <- f.s$ff
    s  <- f.s$s
    s1 <- s - 1
    qs <- qs(x,df,ncp, f.s = f.s)
    return(z = df* s1*s1 *0.5 * qs - log(qs / ff) )
}

z.f <- function(x, df, ncp)
{
  ## Purpose: compute the 'z' (before SQRT) in the 'f'irst wienergerm approx.
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date: 26 Jan 2004, 21:11
  z <- z0(x, df, ncp)
  mu2 <- ncp / df
  d2 <- 1 + 3 * mu2
  return(z + 2./9. * d2 * d2 / (df * (1 + 2*mu2)^3))
}

z.s <- function(x, df, ncp, verbose = getOption("verbose"))
{
    ## Purpose: compute the 'z' (before SQRT) in the 's'econd wienergerm approx.
    ## ----------------------------------------------------------------------
    ## Author: Martin Mächler, Date: 26 Jan 2004, 21:17

    ## the same as z0():
    f.s <- sW(x, df, ncp)
    ff <- f.s$ff
    s  <- f.s$s
    s1 <- s - 1
    qs <- qs(x,df,ncp, f.s = f.s)
    z0 <- df* s1*s1 *0.5 * qs - log(qs / ff)
    ## = z0(x, df, ncp)

    ## further:
    m2s <- (ncp / df) * s # includes Taylor when needed
    f2 <- ff * ff

    d2 <- 1 + 3 * m2s
    d.h <- -1.5 * (1 + 4 * m2s)/ f2 + 5./3.* d2*d2 /(ff*f2)

    h1s <- -h(1 - 1 / s) ## /* == h(1 - s) ~= (1-s)/6  for s ~=1

### LAST(?) FIXME:
    ## know that "t1 + t2" gives cancellation problems for s -> 1 :
    eta <- 1 - ff/qs # =(ff - 2*h1s - s * ff) / (ff - 2*h1s) # ~= 0 when s ~= 1
    gg <- eta / (s1 * s1* ff)
    t2 <- gg* (3. - (0.5 + h(eta)) * eta)
    t1 <- 2 * d2 / (s1 * f2)
    if (verbose)
        cat(sprintf(
                    "c(f= %g, s= %g, h= %g, qs=%g, z0=%g, eta=%g, g=%g, t1=%g, t2=%g)\n",
                    ff,s, h1s, qs, z0, eta, gg, t1, t2))

    return( z0 + 2*(d.h + t1 + t2) / df)
}



pchisqW.R <- function(x, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                      variant = c("s", "f"), verbose = getOption("verbose"))
{
###---- R version of 'pchisq.W' --------

    f.s <- sW(x, df, ncp)
    ##     =
    ff <- f.s$ff
    s  <- f.s$s

    if (length(x) > 1) stop("'x' must be scalar (length 1) !")
    if (x <= 0)  return(p = if(lower.tail) 0 else 1)
    variant <- match.arg(variant)

## /* start calculation */ (translated from C-

    mu2 = ncp / df

    if (s == 1) {
        warning("s == 1 --- returning p=0.5 -- we can do better!")
	return(p = 0.5)
    }
    s.1 = s - 1
    if(abs(s.1) < 1e-10) warning("s ~= 1: lost at least 10 digits of accuracy")
    h.1s = -h(1 - 1 / s)## /* == h(1 - s)
    ##h.1s = h(1 - s)
	                  ##  * MM: FIXME: depends on s which form is better */
    z = df * s.1*s.1 * (0.5 / s + mu2 - h.1s / s) - log((1 - 2 * h.1s / ff) / s)

    if (verbose)
        cat(sprintf("c(ff= %g, s= %g, h= %g, z=%g", ff,s, h.1s, z))

    if (variant == 'f') {##/* improved first order approximation */

	d.2 = 1 + 3 * mu2
	z <- z + 2./9. * d.2 * d.2 / (df * (1 + 2*mu2)^3)
    }
    else { ## /* variant = 's' : the second order approximation */

        m2s <- (ncp / df) * s # includes Taylor when needed
	f2 = ff * ff
	d.2 = 1 + 3 * m2s
	d.h = -1.5 * (1 + 4 * m2s)/ f2 + 5./3.* d.2*d.2 /(ff*f2)

	eta = ff - 2. * h.1s
	eta = (eta - s * ff) / eta

	ff = eta / (s.1 * s.1) / ff
	z <- z + 2*(d.h + 2 * d.2 / s.1 / f2 +
                    ff* (3. - (0.5 + h(eta)) * eta)) / df
    }
    if (verbose)
        cat(sprintf(" zz= %g)\n", z))

    z = sqrt(abs(z))## z < 0 is possible : e.g. x=ncp=2, df=1
    ## but the || does not help much!
    if (s < 1)
	z = -z

    pnorm(z, 0, 1, lower.tail, log.p)
}


##' Fortran/C version: ---> below for  pchisqW() -- the C only version
##'                                    =========
pchisqW. <- function(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                     Fortran = TRUE, variant = c("s", "f"))
{
    ## Purpose: Wiener germ approximation to pchisq(*, ncp=.)
    ## ----------------------------------------------------------------------
    ## Arguments: ... as for pchisq()
    ##		Fortran: logical: use Fortran or C
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 14 Jan 2004, 13:24
    variant <- match.arg(variant)
    ivariant <- match(variant, c("f", "s")) # f[irst] => 1; s[econd] => 2; int to allow > 2
    if(length(df) > 1) stop("'df' must be length 1")
    if(length(ncp) > 1) stop("'ncp' must be length 1")
    df <- as.double(df)
    ncp <- as.double(ncp)
    storage.mode(q) <- "double"
    if(Fortran)
        lapply(q, function(x)
               .Fortran(C_noncechi,
                        variant = ivariant,
                        argument = x,
                        noncentr = ncp,
                        df = df,
                        p = double(1L),
                        ifault = integer(1L))
               )
    else
        ## "nonc_chi" is not callable by .C() anymore, but C_p.. is vectorized in q
        ## we simulate the (Fortran) behavior for easier testing
	lapply(.C(C_pchisqV,
		  q, #    ^ vectorized
		  length(q),
		  noncentr = ncp,
		  df = df,
		  lower = as.logical(lower.tail),
		  log.p = as.logical(log.p),
                  variant = ivariant),
	       function(x) list(p = x))
}

##' a more usefully vectorized -- in 'q' --- version
pchisqV <- function(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                     Fortran = TRUE, variant = c("s", "f")) {
    vapply(pchisqW.(q, df=df, ncp=ncp, lower.tail=lower.tail, log.p=log.p,
                    Fortran=Fortran, variant=variant),
           `[[`, numeric(1), "p")
}


### Vectorized (in 'q' only)  (C only) version :
pchisqW <- function(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                    variant = c("s", "f"))
{
    ## Purpose: Wiener germ approximation to pchisq(*, ncp=.)
    ## ----------------------------------------------------------------------
    ## Arguments: ... as for pchisq()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 26 Jan 2004, 10:10

    variant <- match.arg(variant)
    ivariant <- match(variant, c("f","s")) # f[irst] => 1 /  s[econd] => 2
    df  <- as.double(df)
    ncp <- as.double(ncp)
    if(length(df)  > 1L) stop("'df' must be length 1")
    if(length(ncp) > 1L) stop("'ncp' must be length 1")
    storage.mode(q) <- "double"
    .C(C_pchisqV,
       q, # input *and* output
       length(q),
       noncentr = ncp,
       df = df,
       lower = as.logical(lower.tail),
       log.p = as.logical(log.p),
       variant = ivariant)[[1]]
}
