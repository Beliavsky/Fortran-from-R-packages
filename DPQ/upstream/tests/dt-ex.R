#### dt() : Density of t-dist  ----------- was Martin's ~/R/MM/NUMERICS/dpq-functions/dt-ex.R (since 2002)
#### ==-----=================

### -> ...../MM/..../pt-ex.R for the cumulative dist. density
### -> ...../MM/..../qt-ex.R for its inverse
(doExtras <- DPQ:::doExtras())
(noLdbl <- (.Machine$sizeof.longdouble <= 8)) ## TRUE when --disable-long-double
options(width = 100, nwarnings = 1e5)

###============================================================================
###========= 1. Numerics about the integration constant for nu -> Inf =========
###============================================================================

### The simple naive gamma ratio -- breaks down numerically for nu ~= 10^17
lg.ratio <- function(nu) lgamma((nu+1)/2)- lgamma(nu/2)
gratio <- function(nu) exp(lg.ratio(nu))

## Factor(nu) in dt() density :
## c(nu) = gamma((nu+1)/2) / (gamma(nu/2) * sqrt(nu))  --->  1/sqrt(2)

### and the same with  log(c(nu)) {for direct log density!}
ldt.fact.naive <- function(df) lgamma((df+1)/2)- lgamma(df/2) - log(df)/2
dt.fact.naive <- function(df) exp(ldt.fact.naive(df))

dt.fact0 <- function(df) (1 - (1 - 1/(8*df))/(4*df)) / sqrt(2)

dt.fact <- function(df, cut.val = 200) {
    ## Use asymptotic expansion  only for  df -> oo
    ifelse(df > cut.val,
           dt.fact0     (df),
           dt.fact.naive(df))
}

## on this range, they look already the same!
curve(dt.fact.naive(x), 5, 1000) ; abline(h= sqrt(1/2), col = "gray")
curve(dt.fact0     (x), 5, 1000, col = 2, add = TRUE)
## still :
curve(dt.fact (x), 100, 1e7, log = 'x')
curve(dt.fact0(x), 100, 1e7, col = 2, add = TRUE)
## look at difference
p.dtdiff <- function(df, log='x', cut = 5000)
{
    dd <- dt.fact(df, cut= cut) - dt.fact.naive(df)
    plot(df, dd, ylab = "dt.fact(*, cut) - dt.fact.naive()",
         ylim = quantile(dd, c(.001,.999)),
         type = 'l', col = 2, log = log)
    abline(v=cut, col = "green", lty = 2)
    abline(h=0, col = "gray", lty=3)
}
p.dtdiff(df = seq(5,1000, by=1/2), log= "") # all "zero"
p.dtdiff(df = 2^seq(9,16, len=10001), cut = 800)
## --> at cut: bias, then noise from about 5000 -- good pic!
## maximally about 6e-11

p.dtdiff(df = 2^seq( 9,30, len=10001))## noise up to 2e-6
p.dtdiff(df = 2^seq(20,45, len=10001))## noise up to 0.06
p.dtdiff(df = 2^seq(20,50, len=10001))## break down after 1e14


ldt.fact0 <- function(df) {
    ## formula from Maple (instead of Abramowitz & Stegun):
    t <- 1/(2*df)
    -.5* log(2) + t * (-1/2 + t*t/3)
}
ldt.fact1 <- function(df) {
    ## formula from Maple  -- 1 term more than ldt.fact0 :
    t <- 1/(2*df)
    -.5* log(2) + t * (-1/2 + t*t *(1/3 - 8/5*t*t))
}
ldt.fact <- function(df, cut.val = 200) {
    ## Use asymptotic expansion  only for  df -> oo
    ifelse(df > cut.val,
           ldt.fact1     (df),
           ldt.fact.naive(df))
}

## at the very beginning: difference
curve(ldt.fact.naive(x), .5, 10000, log = 'x', n = 2001)
abline(h= -log(2)/2, col = "gray")
curve(ldt.fact0     (x), col = 2, add = TRUE, n = 2001)
curve(ldt.fact1     (x), col = 3, add = TRUE, n = 2001)

## fact1 is `better' :
curve(abs(ldt.fact0(x) - ldt.fact.naive(x)), 2, 10000, n = 2001, log='x',
      main = "Absolute error")
curve(abs(ldt.fact1(x) - ldt.fact.naive(x)), col = 3, add = TRUE, n = 2001)
## log-zooming: -- see very small l*naive() error coming in as well
curve(abs(ldt.fact0(x) - ldt.fact.naive(x)), 2, 10000, log = 'xy', n = 2001,
      main = "Absolute error")
curve(abs(ldt.fact1(x) - ldt.fact.naive(x)), col = 3, add = TRUE, n = 2001)
## rel.err: ~ the same
curve(abs(1 - ldt.fact0(x) / ldt.fact.naive(x)), 2, 10000, log = 'xy', n = 2001,
      main = "Relative error")
curve(abs(1 - ldt.fact1(x) / ldt.fact.naive(x)), col = 3, add = TRUE, n = 2001)

## on this range, they look the same
curve(ldt.fact.naive(x), 10, 1e10, log = 'x', n = 2001)
abline(h= -log(2)/2, col = "gray")
curve(ldt.fact0     (x), col = 2, add = TRUE, n = 2001)
curve(ldt.fact1     (x), col = 3, add = TRUE, n = 2001)

## now watch breakdown of naive:
curve(ldt.fact.naive(x), 10, 1e18, log = 'x', n = 2001)## 1e14 !

curve(ldt.fact.naive(x), 10, 1e18, log = 'x', n = 2001, ylim = c(-.35,-.34))
## even from 1e9
curve(ldt.fact0     (x), col = 2, add = TRUE, n = 2001)
curve(ldt.fact1     (x), col = 3, add = TRUE, n = 2001)

### note  dt.1() checking below suggest to use ldt.fact1() down to about 100!


### Use maple to get much longer expansion:
### /u/maechler/maple/gamma-expansions.txt
### /u/maechler/maple/gamma-exp2.txt
rr <- seq(0,1.2, by = 1/128)

## c(nu) := GAMMA((nu+1)/2) / GAMMA(nu/2)  * sqrt(2 / nu)   { = above * sqrt(2) }
## r := 1/(4 * nu)
## c() = 1 - r + r*r/2 + 5/2*r^3 - 21/8*r^4 - 399/8*r^5 + 869/16*r^6
##         + 39325/16*r^7 -  334477/128*r^8

cs0 <- function(r) {
    1 - r + r*r/2 + 5/2*r^3 - 21/8*r^4 - 399/8*r^5 +
        + 869/16*r^6 + 39325/16*r^7 -  334477/128*r^8
}

cs <- function(r) {
    ## == cs0 but in Horner form ( + coefficients factored a bit)
    1 + r * (-1 + r/2 *
             (1 +  r *
              (5 + r/4 *
               (-21 + r *
                (-399 + r/2 *
                 (869 + r * (39325 -  334477/8 *r)))))))
}

all.equal(cs0(rr), cs(rr), tol = 1e-12)


## s := 1/(8 * nu) == r / 2
cs2 <- function(s) {
    ## == cs(), using s = r/2 = 1/(8 nu)
    1 + 2*s * (-1 + s *
               (1 +  s *
                (10 + s *
                 (-21 + s * 2 *
                  (-399 + s * (869 + s * (2*39325 - 334477/2 *s)))))))
}

all.equal(cs(rr), cs2(rr/2), tol = 1e-12)# TRUE

### Now using many more terms -- (thanks to maple):
### at the end we see it does not help at all: These terms "diverge"...
css0 <- function(s) {
    1 + 2*s*(-1+ s* ##                  Set of prime factors; "^" := at least ^2
             (1+  s*
              (10+ s* ##                    # 2   5
               (-21+ s* ##                  #   3   7
                (-798+ s* ##                # 2 3   7           19
                 (1738+ s* ##               # 2       11                  79
                  (157300+  s* ##           # 2^  5^  11^ 13
                   (-334477+ s* ##          #         11  13                 2339
                    (-57434806+ s* ##       # 2       11  13 17                11813
                     (119394366+ s* ##      # 2 3   7     13 17 19          677
                      (33601489740+ s* ##   # 2^3^5 7     13 17 19    29  73
                       (-68858583810+ s* ## # 2 3 5          17 19 23       509 607
                        (-28797022447980+ s*# 2^3 5          17 19 23     3607 17911
                         (58526378304180+ s*# 2^3 5          17 19 23      131301607
                          340096557365034000# 2^3 5^         17 19 23 29^ 71  127781
                          ))))))))))))))
}
css <- function(s) {
    1 + 2*s*(-1+ s* ##                    Set of prime factors; "^" := at least ^2
             (1+ s*
              (10+ s* ##                    # 2   5
               (-21+ s* ##                  #   3   7
                (-798+ s* ##                # 2 3   7           19
                 (1738+ s* ##               # 2       11                  79
                  (157300+  s* ##           # 2^  5^  11^ 13
                   (-334477+ s* 2*17* ##    #         11  13                 2339
                    (-1689259+ s* 3*19* ##  # 2       11  13 17                11813
                     (  61607 + s* 5* ##    # 2 3   7     13 17 19          677
                      (3467646 + s* 23* ##  # 2^3^5 7     13 17 19    29  73
                       (-308963 + s* 2* ##  # 2 3 5          17 19 23       509 607
                        (-64604977+ s* ##   # 2^3 5          17 19 23     3607 17911
                         (131301607+ s* ##  # 2^3 5          17 19 23      131301607
                          76299312910 ##    # 2^3 5^         17 19 23 29^ 71  127781
                          ))))))))))))))
}


all.equal(cs(rr), css(rr/2))# not at all!
all.equal(css(rr/2), css0(rr/2))## neither 8.99


p.cs <- function(nu, col = 1:k)
{
  ## Purpose: plotting
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Apr 2002, 20:42

    r <- 1/(4*nu)
    s <- r/2
    mat <- cbind(dt.fact.naive(nu)*sqrt(2)
                 ,dt.fact0    (nu)*sqrt(2)
                 ,cs0(r)
                 ,cs (r)
                 ,cs2 (s)
                 ,css0(s)
                 ,css (s)
                 )
    k <- ncol(mat)
    ## one sees an initial departure ( nu < 3)
    matplot(nu, mat, type = 'b', log = 'x', col = col,
            ylim = c(0.8, 1))# rrange(mat))##
    legend(quantile(nu,.3), mean(par("usr")[3:4]),
           c("dt.fact.naive", "dt.fact0", "cs0", "cs", "cs2", "css0", "css"),
           lty=1:k, col=col, pch = paste(1:k))
    mtext("r = 1 / (4 nu)", line = 2)
    pr <- 10^pretty(log10(r))
    axis(3, at = 1/(4*pr), labels = formatC(pr,dig=1))
    invisible(list(nu=nu,mat=mat, pr=pr))
}

(dfs <- 2^seq(0,20, by = 1/16))[1:10]
## different  at the beginning ; the longer the series the *worse* !!!!
p.cs(dfs)
## Now we see the breakdown of the naive dt.fact.0():
p.cs(nu = 2^seq(0,64, by = 1/4))

###-- Application to dt() :

dt.naive <- function (x, df, log = FALSE)
{
    n1h <- (df+1)/2
    rt <- sqrt(pi*df)
    if(log)
        lgamma(n1h) - lgamma(df/2) - log(rt)  - n1h * log1p(x^2/df)
    else ## following worse than exp(..above..) ?
        gamma(n1h)/(gamma(df/2) * rt) * (1 + x^2/df)^-n1h
}

dt.1 <- function (x, df, log = FALSE, cuts = c(log=100, 1000))
{
    ## Use asymptotic expansion for  df -> oo
    n1h <- (df+1)/2
    if(log)
        ldt.fact(df, cut= cuts[1]) -.5 *log(pi)  - n1h * log1p(x^2/df)
    else
        dt.fact(df, cut=cuts[2])/sqrt(pi) * (1 + x^2/df)^-n1h
}

curve(dt(x=9, df=x), 1, 100)
curve(dt(x=9, df=x), 1, 100, log = 'x')
curve(dt      (x=9, df=x, log = TRUE), 1, 400, log = 'x')
## no visibile difference:
curve(dt.naive(x=9, df=x, log = TRUE), 1, 400, log = 'x', add = TRUE, col=2)

curve(dt      (x=9, df=x, log = TRUE),.25,1000, log = 'x')
## no visibile difference:
curve(dt.naive(x=9, df=x, log = TRUE),.25,1000, log = 'x', add = TRUE, col=2)

### log density comparison :

### "smallish"  nu :
(dfs <- seq(1,500, by = 1/4))[1:10]
plot(dfs,
     dt      (9,df=dfs, log = TRUE) -
     dt.naive(9,df=dfs, log = TRUE),  type = 'l', col = 2)
## increasing upto 6e-13 (1000); 4e-12 (5000)
lines(dfs,
      dt  (9,df=dfs, log = TRUE) -
      dt.1(9,df=dfs, log = TRUE, cuts = 50), col = 3)# 50 too small!
lines(dfs,
      dt  (9,df=dfs, log = TRUE) -
      dt.1(9,df=dfs, log = TRUE, cuts = 80), col = 4)

### "medium"  nu :
(dfs <- seq(200, 10000, by = 1/4))[1:10]
plot(dfs,
     dt      (9,df=dfs, log = TRUE) -
     dt.naive(9,df=dfs, log = TRUE),  type = 'l', col = 2)
## increasing upto 1e-11
lines(dfs,
      dt  (9,df=dfs, log = TRUE) -
      dt.1(9,df=dfs, log = TRUE, cuts = 80), col = 4)

## larger -- using log nu
(dfs <- 2^seq(15, 25, len = 2001))[1:10]
plot(dfs,
     dt      (9,df=dfs, log = TRUE) -
     dt.naive(9,df=dfs, log = TRUE),  type = 'l', log = 'x')
## increasing upto 1e-9 (nu = 1e6), 6e-8 (nu= 4e+7)
lines(dfs,
      dt  (9,df=dfs, log = TRUE) -
      dt.1(9,df=dfs, log = TRUE, cuts = 80), col = 4)

## LARGE for break down of  dt.naive():
dfs <- 2^seq(15, 64, len = 2001)
plot(dfs, ylim = 1e-5*c(-1,1),
     dt      (9,df=dfs, log = TRUE) -
     dt.naive(9,df=dfs, log = TRUE),  type = 'l', log = 'x')
lines(dfs,
      dt  (9,df=dfs, log = TRUE) -
      dt.1(9,df=dfs, log = TRUE, cuts = 80), col = 4)

## log in y as well: -- nice pic:
require(sfsmisc) # mult.fig(), p.datum(), ..
dfs <- 2^seq(5, 64, len = 2001)
op <- mult.fig(3, main = "dt(x, df -> oo)  testing")$old.par
for(x in c(-9, 0, 90)) {
    plot(dfs, main = paste("x = ", format(x)), ylim = c(1e-16, 1e+2),
         abs(dt      (x,df=dfs, log = TRUE) -
             dt.naive(x,df=dfs, log = TRUE)),  type = 'l', log = 'xy')
    lines(dfs, ## many 0's --> missing for y-log
          abs(dt  (x,df=dfs, log = TRUE) -
              dt.1(x,df=dfs, log = TRUE, cuts = 80)), col = 4)
}
p.datum(); mtext(file.path(getwd(),"dt-ex.R"),
                 side = 1, line = 2.5, adj = 1, cex=.8)
par(op)


###============================================================================
###========= 2. t - Distributions, scaled to Variance == 1 :
###============================================================================

###--- t-distributions -- scaled to Var = 1 :
p.tdensity <-
    function(nu, nout = 501, add = FALSE, col = 2, lty = 1,
             xmax = if(add) par("usr")[2] else 6,
             xmin = if(add) par("usr")[1] else -xmax,
             main = substitute(t[n] * " - distribution, scaled to Var = 1",
                               list(n = nu)),
             ylim = c(0, max(y)), ...)
{
    if(nu <= 2) stop("Var(t_{nu}) =  oo  for nu <= 2")
    v <- if(nu < 1/.Machine$double.eps) nu / (nu - 2) else 1 ## = Var(t[nu])
    s <- sqrt(v)
    x <- seq(xmin,xmax, len = nout)
    y <- dt(x * s, df = nu) * s
    if(add)
        lines(x, y, col = col, lty = lty)
    else {
        plot(x, y, type = "l", col = col, lty = lty, ylim = ylim,
             main = main, ...)
        abline(h = 0, lty = 3, col = "gray20")
    }

    invisible(list(x=x, y=y))
}

leg.dens <- function(x = u[1] + (u[2]-u[1])/32,
                     y = u[4] - (u[4]-u[3])/32, lty = 1)
{
    abline(h = 0, lty = 3, col = "gray20")
    u <- par("usr")
    legend(x, y, paste("nu =",nus), lty = lty, col = cols)
    p.datum()
    mtext("/u/maechler/R/MM/NUMERICS/dpq-functions/dt-ex.R",
          side=4, adj=1, cex=.75)
}


nus <- c(2.2, 2.5, 3:5,7,10, Inf)
pal <- palette()
pal[pal == "white"] <- "gray70"
old.pal <- palette(pal)

cols <- 1+ seq(nus)# or something better

if(!dev.interactive(orNone=TRUE)) ## evaluate manually :
    pdf("dt-ex_t-dens.pdf")

p.tdensity(nu = nus[1], col = cols[1],
           main = expression(t[nu] * " - distributions, scaled to Var = 1"))
for(j in 2:length(nus))
    p.tdensity(nu = nus[j], add = TRUE, col = cols[j])
leg.dens()

## Only close to y = 0:
p.tdensity(nu = nus[1], col = cols[1], ylim = c(0, 0.02), xmax = 12,
           main = expression(t[nu] * " - distributions, scaled to Var = 1"))
for(j in 2:length(nus))
    p.tdensity(nu = nus[j], add = TRUE, col = cols[j], lty=j)
leg.dens(lty = seq(nus))

## y ~ 0 -- use log-log-scale and x > 0
xl <- 0.1; xU <- 50
p.tdensity(nu = nus[1], col = cols[1], ylim = c(1e-8, 1),
           xmin=xl, xmax = xU, log="xy",
           main = expression(t[nu] *
               " - distributions, scaled to Var = 1; log-log scale, x > 0"))
for(j in 2:length(nus))
    p.tdensity(nu = nus[j], add = TRUE, col = cols[j], xmin=xl,xmax=xU)
leg.dens(xl, 0.01)

### Now the same graphic for the usual unscaled t's:
x <- seq(-6,6, len = 201)
matplot(x, outer(x, nus, dt), type = "l", lty = 1, col = cols,
        main = expression(t[nu] * " - distributions (unscaled)"))
leg.dens()

palette(old.pal)# restoring to original palette


## Really small  df < 1  and even df << 1 =====================
##               ======-----      =======
## notably as stirlerr(df) used lgamma(1+df) ...

## ====> for now in  ~/R/Pkgs/Rmpfr/tests/special-fun-ex.R
##                            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###------------------- Non-central -------------------------------

## source("t-nonc-fn.R") ## __ DONE:  Replaced by library(DPQ)
require(DPQ)
##
if(!require(sfsmisc))
    lseq <- function (from, to, length)
         2^seq(log2(from), log2(to), length.out = length)

tst <- c(1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 0)
dt(tst, df = 16, ncp=1)
## [1] 0.2380318 0.2398082 0.2220446 0.0000000 0.0000000 0.0000000 0.2382217

x <- lseq(1e-3, 1e-33, length= 301)

plot(x, dt(x, df=16, ncp=1  ), type = "o", cex=.5, log = "x")
plot(x, dt(x, df=16, ncp=0.1), type = "o", cex=.5, log = "x")
plot(x, dt(x, df= 3, ncp=0.1), type = "o", cex=.5, log = "x")
plot(x, dt(x, df= 3, ncp=0.1, log=TRUE), type = "o", cex=.5, log = "x")
plot(x, dt(x, df= 3, ncp=0.001), type = "o", cex=.5, log = "x") ## <- "noise" around 1e-8
## This shows that there's room for more improvement:
curve(  dt(x, df= 3, ncp=0.001), 1e-20, 1.2e-3, log = "x", n=2^10)## ditto
curve(  dt(x, df= 3, ncp=0.001), 1e-20, 1e-5, log = "x", n=2^10)## ditto
cc <- curve(  dt(x, df= 3, ncp=0.001),  1e-8, 1e-5, log = "x", n=2^10)## ditto
## --- but *NOT* a big problem:  accuracy (rel.error) still around 10^-9
plot(x, dt(x, df=.03, ncp=1), type = "o", cex=.5, log = "x")

## MM: Check the new direct formula dtR() -- notably with Rmpfr
require(Rmpfr) || { warning("need package 'Rmpfr' from here on"); quit("no") }
stopifnot(all.equal(dntJKBf(mpfr(0, 64),  5,10), ## gave NaN
                    3.66083172640611114864e-23, tolerance=1e-20))

dntM <- dntJKBf # from DPQ, should just work

system.time(dt.5.10 <- dntM(mpfr(-4:4, 256), 5, 10))# 2.158 sec on lynne[2014]; 2023: 0.568
dt.5.10 ##--> heureka!

if(doExtras && dir.exists(M_dir <- "~/R/MM/NUMERICS/dpq-functions/t-nonc_mathematica"))
  withAutoprint({
    ## compare with the Mathematica computations {documented to be "arbitrary exact"}
    str(dtM.5.10 <- read.table(file.path(M_dir, "fp-noncent_x_5_10.out"),
                               row.names = 1, col.names=c("","dtM"),
                               colClasses="character"))
    x.dtM <- rownames(dtM.5.10)
    (dtM.5.10 <- as(dtM.5.10[,"dtM"], "mpfr"))
    system.time(dt.5.10 <- dntM(mpfr(x.dtM, 256), 5, 10))# {M=1000} --> 5.1s on lynne{2021} (21.7s lynne{'14})
    all.equal(dt.5.10, dtM.5.10) # TRUE
    all.equal(dt.5.10, dtM.5.10, tol = 1e-15) # TRUE
    all.equal(dt.5.10, dtM.5.10, tol = 1e-20) # TRUE
    all.equal(dt.5.10, dtM.5.10, tol = 1e-30) # "Mean relative difference: 3.664....e-23"
  })

### --- this is from 2006 .. back to 2002

log.dnt <- function(x, df, ncp) {
    ## dt(*, ncp) { ~/R/D/r-devel/R/src/nmath/dnt.c }
    ## uses  pt(*, ncp) internally -- for x !=0
    stopifnot(length(df) == 1, length(ncp) == 1)

    r <- x
    is.sml <- abs(x) < sqrt(df * .Machine$double.eps)
    xb <- x[ib <- !is.sml]
    pt1 <- pt(xb*sqrt((df+2)/df), df=df+2, ncp=ncp)
    pt2 <- pt(xb,                 df=df, ncp=ncp)
    r[ib] <- log(df) - log(abs(xb)) + log(abs(pt1 - pt2))
    r[is.sml] <- dt(0, df = df, ncp = ncp, log = TRUE)
    r
}
stopifnot(all.equal(log.dnt(x, df=3, ncp=0.1),
                        dt (x, df=3, ncp=0.1, log = TRUE)))

plot(x, log.dnt(x, df= 3, ncp=0.001), type = "o", cex=.5, log = "x")

ldt <- log.dnt(x, df= 3, ncp=0.001)
(ry <- rrange(ldt, 2)) + 1.00088935 # -1.08e-9  2.37e-9
plot(x, ldt, type = "o", cex=.5, log = "x", ylim = ry)

x <- lseq(1e-10, 1e-4, 2001)
plot (x, log.dnt(x, df= 3, ncp= 1e-5), type = "l", cex=.5, log = "x")
## tons of pnt() "full precision was not achieved" warnings:
lines(x, log.dnt(x, df= 3, ncp= 1e-7), col = "blue")
lines(x, dt(x, df= 3, log=TRUE), col = 2)## central t : ncp = 0; has no problem

df <- 3
ncp <- 0.1

dnt.stats <- function(x, df, ncp) {
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments: as dt()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 May 2006, 12:25

    ## simple vectorization:
    n <- max(length(x), length(df), length(ncp))
    if(n > 1) {
        x   <- rep(x,   length = n)
        df  <- rep(df,  length = n)
        ncp <- rep(ncp, length = n)
    }
    ## The relevant quantity (whose abs|difference| needs to become accurate:
    cbind(pt1 = pt(x*sqrt((df+2)/df), df=df+2, ncp=ncp),
          pt2 = pt(x, df=df, ncp=ncp))
}

del.pt <- function(x, df, ncp)
    pt(x*sqrt((df+2)/df), df=df+2, ncp=ncp) -  pt(x, df=df, ncp=ncp)


cbind(x, dnt.stats(x, df = 3, ncp = 0.1))
##--> Aha: pt1 & pt2 are very close, and there's  `` full cancellation''
dt.s <- dnt.stats(x, df=3, ncp=0.1)
plot(x, abs(dt.s[,2] - dt.s[,1]))



## now have almost the original "bad dt(*, ncp) picture" from the beginning:
plot(x, abs(dt.s[,2] - dt.s[,1])/abs(x), log = "x")

plot (x, del.pt(x, df=3, ncp=1) / x, type = "l", col=2, log = "x")
## almost the same for x < 0 :
lines (x, del.pt(-x, df=3, ncp=1) / (-x), col = 5, lty=2)
lines(x, del.pt(x, df=2, ncp=1) / x, col=3)
lines(x, del.pt(x, df=2, ncp=2) / x, col=4)
lines(x, del.pt(x, df=2, ncp=10) / x, col=4)
lines(x, del.pt(x, df=.2, ncp=10) / x, col=4)

x <- lseq(1e-14, 1e-3, length=1001)
plot (x, del.pt(x, df= 3, ncp=1) / x, type = "l", col=2, log = "x", ylim = c(0,0.2))
lines(x, del.pt(x, df= 2, ncp=1) / x, col=3)
lines(x, del.pt(x, df= 2, ncp=2) / x, col=4)
lines(x, del.pt(x, df= 2, ncp=10) / x, col=4)
lines(x, del.pt(x, df=20, ncp=10) / x, col=5)

plot (x, dt(x, df= 3, ncp=1) , type = "l", col=2, log = "x", ylim = c(0, 0.8))
lines(x, dt(x, df= 2, ncp=1) , col=3)
lines(x, dt(x, df= 2, ncp=2) , col=4)
lines(x, dt(x, df= 2, ncp=10), col=4)
lines(x, dt(x, df=20, ncp=10), col=5)
lines(x, dt(x, df=1e-2, ncp=10), col=6)

mult.fig(20, marP=-c(0,0,2,1))
for(n in 1:20) {
    plot (x, dt(x, df= 3, ncp=rlnorm(1)),
          type = "l", col=2, log = "x", ylim = c(0, 0.45))
    for(i in 1:100) {
        df <- 10/runif(1, 0.1,100); ncp <- rlnorm(1)
        lines(x, dt(x, df=df, ncp=ncp))
        if(!isTRUE(ae <- all.equal(dt(1e-8, df=df, ncp=ncp),
                                   dt(0,    df=df, ncp=ncp), tol = 1e-7)))
            cat(sprintf("df=%g, ncp=%g : not equal: %s\n", df, ncp, ae))
    }
}

df <- lseq(1e-4, 100, length=1001)
plot (df, dt(1e-5, df=df, ncp= 3), type = "l", log="xy", ylim = c(1e-9, 1.5))
lines(df, dt(1e-5, df=df, ncp= 1e-5),col=2)
lines(df, dt(1e-5, df=df, ncp= 0.1), col="pink")# "same" as ncp= 10 ^ -5
lines(df, dt(1e-5, df=df, ncp= 1),   col="pink")
lines(df, dt(1e-5, df=df, ncp= 1.5), col=3)
lines(df, dt(1e-5, df=df, ncp= 5),   col=4)
lines(df, dt(1e-5, df=df, ncp= 20),  col=5)



ncp <- seq(0,  5, length=1001)
plot (ncp, dt(1e-5, df= 100, ncp= ncp), type = "l", log="y")
##-> precision warning from 'pnt' !
## or
plot (ncp, dt(1e-5, df= 100, ncp= ncp), type = "l")
lines(ncp, dt(1e-5, df=  10, ncp= ncp), col=2)
lines(ncp, dt(1e-5, df=   3, ncp= ncp), col=3)
lines(ncp, dt(1e-5, df=   1, ncp= ncp), col=3)
lines(ncp, dt(1e-5, df= 0.1, ncp= ncp), col=4)

plot (ncp, dt(-1e-5, df= 100, ncp= ncp), type = "l")
lines(ncp, dt(-1e-5, df=  10, ncp= ncp), col=2)
lines(ncp, dt(-1e-5, df=   3, ncp= ncp), col=3)
lines(ncp, dt(-1e-5, df=   1, ncp= ncp), col=3)
lines(ncp, dt( 1e-5, df= 0.1, ncp= ncp), col=4)



##==> Explore pt(x, ) for small x a bit -- is there a simple "asymptotic" ( x --> 0 ) ?
##    ---------------
## --> Rather see ./t-nonc-tst.R (algo. in R!)  and  ./pnt-prec.R
##                ~~~~~~~~~~~~~~                     ~~~~~~~~~~~~

### NOTE:  pt() is fine ---
ncp <- lseq(1e-10,10, length=201)## for log-log plot below
## or
ncp <- seq(0,50, length=201)

plot(ncp, -pt(1e-14, df=3, ncp = ncp, log = TRUE), log =  "y", type="b", cex=0.5)
## or log-log:
plot(ncp, -pt(1e-14, df=3, ncp = ncp, log = TRUE), log = "xy", type="b", cex=0.5)
## Warning message:
## full precision was not achieved in 'pnt' <<<<< almost always here!!
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lines(ncp, -pt(1e-16, df= 3, ncp = ncp, log = TRUE), col=2)
lines(ncp, -pt(1e-12, df= 3, ncp = ncp, log = TRUE), col=4)
lines(ncp, -pt(1e-12, df= 4, ncp = ncp, log = TRUE), col= "brown")
lines(ncp, -pt(1e-10, df=40, ncp = ncp, log = TRUE), col= "purple")

## 'x' *and* 'df'  do not seem to matter really (for largish ncp'
all.equal(pt(1e-12, df= 4,  ncp = ncp, log = TRUE),
          pt(1e-16, df= 10, ncp = ncp, log = TRUE))# TRUE !!
