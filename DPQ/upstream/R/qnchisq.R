#### Approximations of   qnchisq() -*- delete-old-versions: never -*-
####
#### R simulation of C's qnchisq()  <<< not yet >>
####
####  =========================
####  Function definitions only
####  =========================

###-- 1 -- Normal Approximations ---------------------------------------

##- Date: Wed, 19 May 1999 16:47:23 -0400 (EDT)
##- From: Ranjan Maitra <maitra@math.umbc.edu>
##- To: S-news <s-news@wubios.wustl.edu>
##- Subject: Re: [S] quantiles of non-central chisq

##- Many thanks to Jim Stapleton for answering my request for a function to
##- calculate the quantiles of a non-central chi-squared distribution.

## MM: Where is this approx. from?

##- The function is as follows:
qchisqAppr.0 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    df + ncp + z.p * (2 * df + 4 * ncp)^0.5
}

qchisqAppr.1 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## (1+b) * q^chisq(x,df);
    ## better (?) to use q^chisq(x, df*); -> qchisqAppr.2
    (1 + ncp/(ncp +df))* df * (1 - 2/(9*df) + z.p*sqrt(2/(9*df)))^3
}

qchisqAppr.2 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## use = q^chisq(x, df*); df* = a/(1+b); a=(df+ncp); b=ncp/(df+ncp)
    a <- df + ncp
    b1 <- 1+ ncp/a # = 1+b
    df. <- a/b1
    alp <- 2/(9*df.)
    a * (1 - alp + z.p *sqrt(alp))^3
}

qchisqAppr.3 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## use = q^chisq(x, df*); df* = a/(1+b); a=(df+ncp); b=ncp/(df+ncp)
    a <- df + ncp
    df. <- a/(1+ ncp/a)
    alp <- 2/(9*df.)
    df. * (1 - alp + z.p *sqrt(alp))^3
}


## -- Cornish-Fisher (inverse Edgeworth) Expansions:
qchisqApprCF1 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## 1-term Cornish-Fisher Expansion
    mu <- df + ncp
    sig2 <- 2*(df + 2*ncp)
    gam1 <- 8*(df + 3*ncp) * sig2^(-3/2)
    mu + sqrt(sig2)* (z.p + gam1/6 * (z.p^2 - 1))
}

qchisqApprCF2 <-
    function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    z.p <- qnorm(p, lower.tail=lower.tail, log.p = log.p)
    ## 2-term Cornish-Fisher Expansion
    mu <- df + ncp
    sig2 <- 2*(df + 2*ncp)
    gam1 <- 8*(df + 3*ncp) * sig2^(-3/2)
    gam2 <- 48*(df + 4*ncp) / sig2^2
    mu + sqrt(sig2)* (z.p + gam1/6 * (z.p^2 - 1) +
                      gam2/24* z.p*(z.p^2-3) - gam1^2/36* z.p*(2*z.p^2 -5))
}



###-- 2 -- Chisq (Central) Approximations ---------------------------

qnchisqPatnaik <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Basically Abramowitz & Stegun "nu > 30" formula  (26.4.30)
    ## Also == Patnaik(1949)'s as reported in Johnson,Kotz,Bala.. p.462
    a <- df + ncp
    b1 <- 1 + ncp / a
    nu. <-  a/b1 # nu* = (nu + lam)^2 / (nu + 2*lam)
    b1 * qchisq(p, df = nu., lower.tail= lower.tail, log.p= log.p)
}

qnchisqPearson <- function(p, df, ncp = 0,
                       lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Pearson(1959) approximation to pnchisq() - using pchisq()
    ##  Very good for df > 2--5  and much better at *right tail*
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Feb 2004, 15:37

    n2 <- df + 2*ncp
    n3 <- n2 + ncp ## = df + 3*ncp
    n2n3 <- n2/n3
    cdf <- n2n3*n2*n2n3 # == n2^3 / n3^2  but with overflow/underflow protection
    -(ncp / n3 * ncp) + qchisq(p, df = cdf, lower.tail=lower.tail, log.p=log.p) / n2n3
}

qchisqCappr.2 <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Johnson,Kotz,Bala.. p.466  (29.66)
    ## is  *VERY* bad for small df or small df/ncp  (df = 0.1; df/ncp < 1};
    ## but quite strong for small ncp, or small ncp/df , even (1 , 1)
    c2 <- qchisq(p, df, lower.tail= lower.tail, log.p= log.p)
    l.n <- ncp / df # lambda / nu
    c2 * (1 + l.n *(1 + 0.5*l.n*(1 - c2 / (df+2))))
}

## Nice inversions of Normal approximations in ./pnchisq.R :
## Abdel-Aty (1954) .. Biometrika
qnchisqAbdelAty <-  function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Abdel-Aty(1954) "first approx." (Wilson-Hilferty) approximation to pnchisq()
    ## Johnson,Kotz,...: Has it *WRONGLY*  in eq. (29.61a), p.463
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 2018-08-16
    ## Explicit inversion of pnchisqAbdelAty()
    r <- df + ncp # = ν + λ  ( = k + λ  on Wikip.)
    n2 <- r + ncp # = ν + 2λ
    V <- 2*n2 / (9*r^2)
    qNp <- qnorm(p, mean = 1-V, sd = sqrt(V), lower.tail=lower.tail, log.p=log.p)
    r * qNp^3
}

qnchisqSankaran_d <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    ## Purpose: Sankaran(1959,1963) according to Johnson et al. p.463, (29.61d):
    ## Explicit inversion of pnchisqSankaran_d()
    r <- df + ncp  #              = ν + λ
    n2 <- r + ncp  # = df + 2*ncp = ν + 2λ
    n3 <- n2 + ncp # = df + 3*ncp = ν + 3λ
    q1 <- n2/r/r   # = (n + 2λ)/(n+λ)²   {with some overflow/underflow protection for large df | ncp}
    ## h1 <- - 2/3 * r*n3/n2^2            with some "protection":
    h1 <- - 2/3 * r/n2 * n3/n2 # =  h - 1
    h <- 1 + h1
    mu <- 1 + h*h1*q1*(1 + (h-2)*(1-3*h)*q1/2)
    V  <- 2*h^2 * q1 *(1 +  h1 * (1-3*h)*q1)
    qNp <- qnorm(p, mean = mu, sd = sqrt(V), lower.tail=lower.tail, log.p=log.p)
    r * qNp^(1/h)
}



###-- 3 -- Inversion Algorithms ---------------------------------------


##-- The next build on using  Newton's method newton() on pchisq() / dchisq();

## This function finds the quantile corresponding to p
## for the noncentral-chisquare distribution with
##  df degrees of freedom and noncentrality parameter delta = ncp.
qchisqN <- function(p, df, ncp = 0, qIni = qchisqAppr.0, ...)
{
    ## FIXME: implement log.p=TRUE  and also lower.tail=FALSE
    newton(x0 = qIni(p, df=df, ncp=ncp),
           G = function(x, z) pchisq(x, df = z[2], ncp = z[3]) - z[1],
           g = function(x, z) dchisq(x, df = z[2], ncp = z[3]),
           z = c(p, df, ncp),
           xMin = 0, xMax = 1, dxMax = 1, ## <== (log.p = FALSE)
           ...)
}

## >>> ../man/newton.Rd <<<
newton <- function(x0, G, g, z,
                   xMin = -Inf, xMax = Inf, warnRng = TRUE,
                   dxMax = 1000, eps = 0.0001, maxiter = 1000L,
                   warnIter = missing(maxiter) || maxiter >= 10L,
                   keepAll = NA) # FALSE, NA, TRUE
{
    ## Given the function G and its derivative g,
    ## newton uses the Newton method, beginning at x0,
    ## to find a point xp at which G is zero. G and g
    ## may each depend on a parameter z.

    ## The result is the 3-component vector of an approximation x* of xp,
    ## G(x*, z), and the number ("iter") of maxiter.
    ## x* satisfies  abs(G(x*, z)) < eps.

    ##  keepAll = TRUE  to also get  the vectors of consecutive values of
    ##  x and G(x, z);

    stopifnot(length(x0) == 1L,
              length(xM <- xMax - xMin) == 1L, xM >= 0,
              length(dxMax) == 1L, dxMax > 0,
              is.logical(keepAll), length(keepAll) == 1L)
    if(is.finite(xMin) && x0 < xMin) {
        if(warnRng) warning(sprintf("x0 < xMin(=%g) --> setting x0:=xMin", xMin))
        x0 <- xMin
    }
    if(is.finite(xMax) && x0 > xMax) {
        if(warnRng) warning(sprintf("x0 > xMax(=%g) --> setting x0:=xMax", xMax))
        x0 <- xMax
    }
    Gx <- G(x0, z)
    if((give.all <- isTRUE(keepAll))) {
        x0vec <- x0
        Gxvec <- Gx
    }
    iter <- 0L
    conv <- TRUE
    while(abs(r <- Gx/g(x0, z)) > eps) { # often ok, even if g(x0,z) == 0
        ## d := r, unless |r| is too large:
        d <- min(abs(r), dxMax) * sign(r)
        x0 <- x0 - d
        if(is.finite(xMin) && x0 < xMin) x0 <- xMin
        if(is.finite(xMax) && x0 > xMax) x0 <- xMax
        Gx <- G(x0, z)
        if(give.all) {
            Gxvec <- c(Gxvec, Gx)
            x0vec <- c(x0vec, x0)
        }
        if((iter <- iter + 1L) > maxiter) {
            if(warnIter) warning("iter >", maxiter)
            conv <- FALSE; break
        }
    }
    if(give.all)
        list(x = x0, G = Gx, it = iter, converged = conv,
             x.vec = x0vec, G.vec = Gxvec)
    else if(is.na(keepAll))
        list(x = x0, G = Gx, it = iter, converged = conv)
    else x0
} ## newton()

### From Johnson et al., p.465-466 f : Formulas by   Bol'shev and Kuzntzov (1963) [Russian]
### Useful for small  λ = ncp: error is  O(λ³)  uniformly in any finite interval of q
## --> formula (29.66), p.466
qnchisqBolKuz <- function(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
    stopifnot(df > 0)
    ## central quantile χ'²_{ν,α} :
    ch <- qchisq(p, df, ncp=0, lower.tail=lower.tail, log.p=log.p)
    lnu <- ncp/df # = λ/ν = λν⁻¹
    ## x* = .....  (29.66)
    ch * (1 + lnu*(1 + .5*lnu*(1 - ch/(df+2))))
}
