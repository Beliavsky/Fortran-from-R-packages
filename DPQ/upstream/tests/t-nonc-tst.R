#### Examples for non-central t
#### ==========================
###
### Johnson, Kotz, ...;  2nd ed, Vol 2.  Chapter 31, around p.520

c({}## This stems from the two files
, "/u/maechler/R/MM/NUMERICS/dpq-functions/pnt-ex.R"
, "/u/maechler/R/MM/NUMERICS/dpq-functions/t-nonc-approx.R"
, "SEE ALSO ./pnt-prec.R (still current problems!)"
, "           =========="
)
## and originally, had
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/t-nonc-fn.R")
##         - - - - - - - - - - - - - - - - - - - - - - - - - -
library(DPQ) ## source("/u/maechler/R/Pkgs/DPQ/R/t-nonc-fn.R")

stopifnot(exprs = {
    require(graphics)
    require(sfsmisc)
})

options(width = 100, nwarnings = 1e5)
source(system.file(package="DPQ", "test-tools.R", mustWork=TRUE))
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()
relErrV <- sfsmisc::relErrV

if(!dev.interactive(orNone=TRUE)) pdf("t-nonc_P-1.pdf")

### Part I  -- former  pnt-ex.R
### ---------------------------

### Find value of delta after which  pt(*,*, ncp=delta) === 0 :
### MM{2014}: This is no longer true: pt(ncp -> Inf) now "goes all the way to underflow"
###           and pt(*, log.p=TRUE) "is ok"

## This would be from theory :
(dmax <- sqrt(-2*log(2* .Machine$double.xmin))) # 37.62189
.5*exp(-dmax^2/2) # 2.225074e-308

pt.pos <- function(delta,q,df) (pt(q=q,df=df,ncp=delta) > 0) - 1/2

## see 'MM{2014}' above: this no longer works
## d1 <- uniroot(pt.pos, q=10,df=2, lower=30, upper= 40, tol=1e-10)$root
## d2 <- uniroot(pt.pos, q= 1,df=2, lower=30, upper= 40, tol=1e-10)$root
## d1 == d2 # TRUE
d3 <- uniroot(pt.pos, q=.1, df=20, lower=30, upper= 40, tol=1e-10)$root
## d4 <- uniroot(pt.pos, q=100,df=20, lower=30, upper= 40, tol=1e-10)$root
## d3 != d4 # T !!
## d1 == d4 # T
## c(d1,d3)
## 38.57550 38.56226
d5 <- uniroot(pt.pos, q=.01,df= 2, lower=30, upper= 40, tol=1e-10)$root
d6 <- uniroot(pt.pos, q=.01,df=20, lower=30, upper= 40, tol=1e-10)$root
c(d3, d5, d6)
## 37.62274 37.62189 37.62189
## unique(sort(c(d1,d2,d3,d4,d5,d6)))
##[1] 38.45745 38.46036 38.56226 38.57550


###--- When  DEBUG_pnt  is compiled in:
pt(20,20,20)
pt(20,20,25)
pt( 1,20,25)
pt(30, 1,25)
pt(100, 1,25)
##                	  p        #{iter}  result
pt(30,30,30)       #  1.84694e-196   600    .4672569
pt(35,35,35)       #  4.93855e-267   786    .4694869
pt(36,36,36)       #  1.88862e-282   827    .4698806
d <- 37  ;pt(d,d,d)#  2.65703e-298   868    .4702595
d <- 38  ;pt(d,d,d)#  1.37516e-314   902    .4706245 -- no longer
d <- 38.5;pt(d,d,d)#  6.91692e-323  (799)  s<0 => error-bound < 0  (!) wrong convergence => .4702131
##---
pt(39,39,39)#   p = 0 ==> result=0 --- no longer

## FIXME:
pt   (40,40, 38.5) # R's
pntR1(40,40, 38.5, verbose=2) ## and then edit ---> gives the *.out file
##--------- FIXME --- work via sink() or capture.output -- to get  rr data.frame()
 if(FALSE) {
t.file <- "/u/maechler/R/MM/NUMERICS/dpq-functions/pnt-40-40-38.5.out"
nr <- names(rr <- read.table(t.file, header= TRUE))
summary(rr)
matplot(rr[,"it"], rr[,-1], type='l', log="y")
##--> Graphic axis "buglet"

par(mfrow=c(4,2)); f <- FF ~ it
for(i in 2:ncol(rr)) {f[[2]] <- as.name(nr[i]); plot(f, data=rr, type='l')}

par(mfrow=c(4,2)); f <- FF ~ it
for(i in 2:ncol(rr)) {f[[2]] <- as.name(nr[i]); plot(f, data=rr, type='l', log='xy')}
}
##----------- end{FIXME}---------------------------------------------

###------------------------ dcdflib ptnc() ---- computations and plots
##			    ======= ~~~~~
##==> ~/R/Pkgs/dcdflib/tests/ ---------

## How does density look in these extremes ?
##          -------
curve(dt(x, df=10,ncp=10000), 0, 85000, col=2, log="y") # log-density
curve(dt(x, df=10,ncp=10000, log=TRUE), 1e3, 1e9, col=2, ylim = c(-30,-8),log="x")
## is it *not* log concave -- or is that a numeric error in the tails?
## ditto ?
curve(dt(x, df=10,ncp= 1e7, log=TRUE), 1e5, 1e11, col=2, ylim = c(-35,-15),log="x")
curve(dt(x, df=10,ncp= 1e7), 5e6, 3e7, col=2)
integrate(function(x) dt(x, df=10,ncp= 1e7), 5e6, 3e7) # mathematically < 1, but
## 1.004937 with absolute error < 8e-08 hm... --> seems error in tails


yl <- c(-80,0)
curve(dt(x, df=10,ncp= 10, log=TRUE), 0, 600, col=2, n=5001, ylim=yl)
curve(dt(x, df=10,ncp= 15, log=TRUE), 0, 900, col=2, n=5001, ylim=yl) -> dt10.15.L

## Log-log
curve(exp(dt(x, df=10,ncp= 15, log=TRUE)), 1, 1e7, log="xy",
      col=2, n=5001, axes=FALSE, ylab = "",
      main="dt(x, df=10,ncp= 15) -- Log-Log scale") -> edtLog
eaxis(1)
op <- par(las=2)
eaxis(2, at = 10^pretty(log10(axTicks(2, log=TRUE)), 10),
      at.small =FALSE)
par(op)
mtext(R.version.string, cex = .6, adj=1)
## NB: Mathematica (e.g) ./pnt-ex.nb  -- shows __linear__ {in log-log} tail
{}

## From: Jerry Lewis <jerry.lewis@biogenidec.com>
## To: Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: R pnt.c function
## Date: Wed, 6 Nov 2013 21:36:16 +0000

## I noticed your comment in the code

## /*----------- DEBUGGING -------------
## *
## *             make CFLAGS='-DDEBUG_pnt -g'

## * -- Feb.3, 1999; M.Maechler:
##                 - For 't > ncp > 20' (or so)              the result is completely WRONG!
## */

## while chasing down accuracy problems for calculating 1-sided normal
## tolerance limits.  For example trying to reproduce values in Odeh &
## Owen's Table 1.1  by
n <- c(102,104); round(qt(.995,n-1,-qnorm(.0001)*sqrt(n))/sqrt(n),3)
## the value R calculates for n==104 is clearly wrong, since the
## k-factor fails to be decreasing in n (inaccurate for all larger n
## as well).

## MM: graphically convincing is e.g.
n <- seq(80,120, by=0.5)
summary(ncp.n <- -qnorm(.0001)*sqrt(n))
qt995 <- qt(.995, n-1, ncp=ncp.n) / sqrt(n)
summary(warnings())# 45 identical warnings:
## full precision may not have been achieved in 'pnt{final}'
plot(n, qt995, type = "l", col=2, lwd=2) ## => visual kink in {102, 103}:
abline(v = 102:103, lty=3, col="gray")

qtSc <- function(n, p = 0.995, delta = .0001)  {
    stopifnot(length(p) == 1, length(delta) == 1, n > 1, delta > 0,
              0 < p, p < 1)
    ncp.n <- -qnorm(delta) * sqrt(n)
    qt(p, n-1, ncp=ncp.n) / sqrt(n)
}

n <- seq(80,120, by=.5); plot(n, qtSc(n), type = "l", col=2, lwd=1.5)
summary(warnings()) # --> "full precision may not have been achieved in 'pnt{final}'"
## It's the *smaller* n's which give the warning:
assign("last.warning", NULL, envir=as.environment(find("last.warning"))) # (non-API)
n <- seq(80,102, by=.5); plot(n, qtSc(n), type = "l", col=2, lwd=1.5)
try(# works "locally" but nut with the whole script
stopifnot(length(n) == length(warnings())) ##
)

## These give *no* warnings [and we know they are inaccurate!]:
n <- seq(104,120, by=.5); plot(n, qtSc(n), type = "l", col=2, lwd=1.5)

## MM: the bug really is in pt(*,df, ncp), i.e., pnt.c:
n <- seq(80,120, by=0.25)
summary(ncp.n <- -qnorm(.0001)*sqrt(n))
qq <- 4.6*sqrt(n)
plot(qq, pt(qq, n-1, ncp=ncp.n), type = "l") ## - yes!
cbind(n, qq, ncp.n, pt = pt(qq, n-1, ncp=ncp.n))## it is clearly at
## ncp =
## even simpler --
curve(pt(37,37,x), 37, 38)

## When x>0 and ncp>0, a cheap fix would substitute pf(x^2,1,df,ncp^2)
## for pt(x,df,ncp) whenever pnorm(-ncp) is negligible
## (as it should be when ncp>20).
## MM: According to
## https://en.wikipedia.org/wiki/Noncentral_t-distribution#Related_distributions
## If T is noncentral t-distributed with ν degrees of freedom and
## noncentrality parameter μ and F = T^2, then F has a noncentral
## F-distribution with 1 numerator degree of freedom, ν denominator degrees of freedom,
## and noncentrality parameter μ^2.
## MM: If t >= 0,
## pf(t^2, df1=1, df2=nu, ncp=mu^2) = P(T^2 <= t^2) = P(|T| <= t) =
## ================================
##   = P(-t <= T <= t) =  pt(t, df, ncp) - pt(-t, df, ncp)
##                        ================================
## and now, indeed, if   pt(-t, df,ncp) << pt(t, df, ncp)  is negligble, we can use pf(...)
## HOWEVER, 1) this not only depends on ncp !!
##          2) Code below suggests that noncentral pf() is not quite ok for small t^2

## MM: let's see:
ptq <-  pt(qq,          df=  n-1, ncp = ncp.n)
pfq2 <- pf(qq^2, df1=1, df2= n-1, ncp = ncp.n^2)
cbind(n, qq, ncp.n, pt = ptq, "pf.q^2" = pfq2, rel.E = 1 - ptq/pfq2)
lines(qq, pfq2, col=adjustcolor("blue", 0.4), lwd=3)
##==> he is right that  pf() is good *here*

chk.t.F <- function(x, df, ncp, tolerance = 1e-9) {
    pt. <- pt(x,           df= df, ncp = ncp)
    pf. <- pf(x^2, df1=1, df2= df, ncp = ncp^2)
    re <- relErrV(pf., pt.)
    if(any(abs(re) > tolerance))
        warning("rel.error too large: ", max(abs(re)))
    invisible(re)
}
chk.t.F(x=qq, df=n-1, ncp=ncp.n) ## too large: 0.00894

## many are not so good..
## these look good:
cols <- c(palette()[3], adjustcolor(2, 0.4))
curve(pt(x, df=11.282, ncp=30), 0, 50, n=512, col=3)
curve(pf(x^2, df1=1, df2=11.282, ncp=30^2), n=512, col=cols[2], lwd=4, add=TRUE)
legend("top", c("pt()", "pf()"), col=cols, lwd=c(1,4), bty="n")
## but not if we look at the left tail in log scale:
##                                       ==========
curve(pt(x, df=11.282, ncp=30), 0, 30, n=512, col=3, log="y",yaxt="n"); eaxis(2)
curve(pf(x^2, df1=1, df2=11.282, ncp=30^2), n=512, col=adjustcolor(2, 0.4), lwd=4, add=TRUE)
legend("topleft", c("pt()", "pf()"), col=cols, lwd=c(1,4), bty="n")
## pf() and pt() very much differ when x gets small even for this
## large ncp

## The difference should *only* be "the other"  pt() part :  pt(-x, df, ncp)
## but that looks quite bad :
## Note that it should be monotone !!!
pxy <- curve(pt(x, df=11.282, ncp=30, log=TRUE), -25, 30, n=512, col=3)
       curve(pt(x, df=5, ncp=10, log=TRUE), -20, 20, n=512, col=2, lwd=2, type="o")

## If we look at the corresponding *density* function, it's not better:
curve(dt(x, df=11.282, ncp=30, log=TRUE), -25, 30, n=512, col=3)

## What about  negative ncp ? [the same phenomenon}
## all fine "regular scale":
curve(dt(x, df=11, ncp=-10), -25, 30, n=512, col=2, lwd=2)
## "all bad" on log scale:
curve(dt(x, df=11, ncp=-10, log=TRUE), -40, 30, n=512, col=2, lwd=2)




## [Jerry, from his e-mail:] ----------------------------------------------------

## That alone would improve the minimum accuracy from 2 to 6 significant
## figures when calculating tolerance factors for the cases that Odeh and
## Owen tabulated.

## The series expansion used by Lenth (basis of R's pnt.c function) can
## be implemented in a way that can accurately handle these cases, but I
## have not digested R's c code to see where the problem lies.  The
## boundary where Lenth's p0 crosses from a normalized to a denormal binary representation is
## between the cases of n==102 and n==104, so I suspect that p0 is not
## factored out of the summation, which would accumulate reduced accuracy
## terms.

## Jerry

if(!dev.interactive(orNone=TRUE)) { dev.off(); pdf("t-nonc_P-2.pdf") }

### Part II -- former  t-nonc-fn.R
### ------------------------------

### b_nu  -----------------------------------
mult.fig(2, main = "b(nu) = E[ Chi_nu ] / sqrt(nu)")$old.par -> opar
plot(b_chi, n=1024, col=2)

## Unfortunately, the above switch gives  a small kink at nu = 300
## 2015-01-02:  b_chi() improved: no kink anymore
##              ================================= !

curve(b_chi, 200,1200)
curve(b_chiAsymp, add=TRUE, col='red')

plot(b_chi,1,10, col='red'); curve(b_chiAsymp, add=TRUE, col='blue')

## New cutoff
plot(b_chi,340.9, 341.1, n=1001)# no jump/kink visible
plot(b_chi,999.9, 1000.1, n=1001)# no jump/kink visible

nu <- seq(200, 1000, length=1001)
plot(nu, b_chi(nu) - b_chiAsymp(nu), type='l')# b_chi(nu) > b_chiAsymp() ALWAYS
## with 'one.minus=TRUE'  the reverse:  b_chi(*, one.minus=TRUE) < b_chiAsymp(*, ..TRUE)
plot(nu, b_chi(nu,one.minus=TRUE) - b_chiAsymp(nu,one.minus=TRUE), type='l')

par(opar)# reset

### ----------- pntJW39() and  pntR() ----------------------

pntR(30,30,30, verbose=2)# 600 iter
pntR(30,30,ncp = 30:40)
pntR(2,   10,ncp=1e5)#> C-code directly {underflow p=0, |ncp| = |delta| too large
pntR(2,df=10,ncp=1e4)#>  (ditto)

## t --> 0 : is it problematic ?
df <- 1
df <- 10
ncp <- 1
x <- 1e-12
x <- 1e-6
(pt1 <- pntR(x * sqrt((df+2)/df), df=df+2, ncp=ncp, errmax = 1e-14))
(pt2 <- pntR(x                  , df=df  , ncp=ncp, errmax = 1e-14))
c(pt1 - pt2, (pt1 - pt2)/x)

if(FALSE) ## useful to see things
debug(pntR)
pntR(1e-8,  df=10, ncp=1, errmax = 1e-14)
pntR(1e-15, df=10, ncp=1, errmax = 1e-14)

plot(function(t)pntJW39(t,30,30))
plot(function(t)pntJW39(t,100,30),-9,5, col="red", log="y")
plot(function(t)pt     (t,100,30),-9,5, col="blue", log="y",add=TRUE)
plot(function(t)pntJW39(t,100,30),0,5, col="red", log="y")#~ lin
plot(function(t)pntJW39(t,100,30),0,15, col="red", log="y")
plot(function(t)pt(t,100,ncp=30),0,15, col="blue", log="y", add=TRUE)

plot(function(t)pntJW39(t,100,30),0,25, col="red", log="y")
plot(function(t)pntJW39(t,100,30),0,25, col="red")

plot(function(t)pntJW39(t,100,30),0,45, col="red")

plot(function(t)pntJW39(t,100,30),20,50, col="red", log="y")
plot(function(t)pt     (t,100,30),20,50, col="blue", log="y",add=TRUE)

o <- par(las=1, mar=.1+c(4,4,4,4))
plot(function(t)abs(pt(t,100,30)-pntJW39(t,100,30)),20,50,
     col="red", log="y"); eaxis(4, at=10^-(2:8))
par(o)

## xtended x-range
o <- par(las=1, mar=.1+c(4,4,4,4))
plot(function(t)abs(pt(t,100,30)-pntJW39(t,100,30)), 17, 120, n=1001,
     col="red", log="y"); eaxis(4)
par(o)


### --------------------- qt(), qtAppr(), qtNappr(), etc ---------------------

## df=1  is pretty bad...
plot(function(t)qt(t,df=1))
plot(function(t)qtAppr(t,df=1,ncp=0),col="blue",add=TRUE)

## df=2  quite a bit better..
plot(function(t)qtAppr(t,df=2,ncp=0),col="blue")
plot(function(t)qt(t,df=2),add=TRUE)

## df=4 a bit better.. still only for alpha ~ in  (.1, .9)
plot(function(t)qtAppr(t,df=4,ncp=0),col="blue")
plot(function(t)qt(t,df=4),add=TRUE)

plot(function(t)abs(1-qt(t,df=4)/qtAppr(t,df=4,ncp=0)), main="rel.Error")
plot(function(t)abs(1-qt(t,df=10)/qtAppr(t,df=10,ncp=0)), main="rel.Error")

## max error: 10e-5,  however... still NaN's sqrt()
plot(function(t)abs(1-qt(t,df=100)/qtAppr(t,df=100,ncp=0)))

plot(function(t)qtAppr(t,df=100,ncp=100))
## catastrophe : !!!!
plot(function(t)qtAppr(t,df=1,ncp=100))


plot(function(t) t - pt     (qtAppr(t,df=4,ncp=100),df=4,ncp=100))

## --> pntJW39() uses MUCH better asymptotic than Abramowitz&Stegun (in pnt).
plot(function(t) t - pntJW39(qtAppr(t,df=4,ncp=100),df=4,ncp=100),
     col='red',add=TRUE)
## Absolute Error: very small -- just proves that the two ".appr"  CORRESPOND!
plot(function(t) t - pntJW39(qtAppr(t,df=4,ncp=100),df=4,ncp=100))


plot(function(t)qtAppr(t,df=10,ncp=1e5))

## Shows that  pt(,, ncp=1e5)  uses asymptotic form alright:
plot(function(t) t - pt     (qtAppr(t,df=10,ncp=1e5),df=10,ncp=1e5))
## --> pntJW39()  MUCH better ========== fitting to qtAppr !!
## ------ not necessarily better asymptotic than Abramowitz&Stegun (in pnt).
plot(function(t) t - pntJW39(qtAppr(t,df=10,ncp=1e5),df=10,ncp=1e5),
     col='red',add=TRUE)


###--- Diverse tests, some lifted from ../R/t-nonc-fn.R ------------------------

##		[MM: the next 2 'dntJKBch(.)' were 'dntRwrong(.)' originally ]

print(f. <- .dntJKBch(1:6, df=3, ncp=5, check=TRUE), digits = 4)
## now:[1] 0.0005083 0.0260607 0.1191377 0.1761045 0.1657718 0.1305411
## was:[1] 0.0032023 0.2728377 1.5174647 2.4481320 2.4106931 1.9481189 -- wrong "but sensible"

x <- seq(-1,12, by=1/16)
if(FALSE) ## FIXME ?! -- internal logic error:
.dntJKBch(x, df=3, ncp=5, check=TRUE)
## ## Error in (function (x, df, ncp, log = FALSE, M = 1000, check = FALSE,  :
## ##   exp(lterms[ii]) and terms[ii] are not equal:
## ##   'is.NA' value mismatch: 0 in current 340 in target

## -- This is now in help page examples  >>> ../man/dnt.Rd <<<
## fx <- dt(x, df=3, ncp=5)                  ~~~~~~~~~~~~~
## re1 <- 1 - .dntJKBch(x, df=3, ncp=5) / fx # with warnings
## re2 <- 1 -  dntJKBf (x, df=3, ncp=5) / fx
## summary(warnings()) ## "In log(x * ncp * sqrt(2)/sqrt(df + x^2)) : NaNs produced"
##           all.equal(re1[!is.na(re1)], re2[!is.na(re1)], tol=0)##  Mean relative ...: 2.068..e-5
## stopifnot(all.equal(re1[!is.na(re1)], re2[!is.na(re1)], tol=1e-6))
## matplot(x, log10(abs(cbind(re1, re2))), type = "o", cex = 1/4)


print(ft <- .dntJKBch(1:6, df=3, ncp=5), digits = 4)
## [1] 0.0005083 0.0260607 0.1191377 0.1761045 0.1657718 0.1305411
## correct !! *with* the factorial !
          all.equal(ft, dt(1:6, df=3, ncp=5), tol = 0) # 1.4959e-12
stopifnot(all.equal(ft, dt(1:6, df=3, ncp=5), tol = 1e-11))


## From: Stephen Berman ..
## To: <r-help@r-project.org>
## Subject: [R] Why does qt() return Inf with certain negative ncp values?
## Date: Mon, 13 Jun 2022 23:47:29 +0200

## Can anyone explain why Inf appears in the following results?

## MM:
options(nwarnings = 5e5)

sapply(-1:-10, function(ncp) qt(1-1*(10^(-4+ncp)), 35, ncp))
 ## [1]  3.6527153  3.0627759  2.4158355  1.7380812  1.0506904  0.3700821
 ## [7]        Inf -0.9279783 -1.5341759 -2.1085213
## There were 1060 warnings (use warnings() to see them)
summary(warnings())
## 1060 identical warnings:
## In qt(1 - 1 * (10^(-4 + ncp)), 35, ncp) :
##   full precision may not have been achieved in 'pnt{final}'

nc1 <- seq(-12, -1, by=1/32)
r1 <- vapply(nc1, function(ncp) qt(1 - 10^(-4+ncp), df=35, ncp=ncp), 1.23)
## There were 8264 warnings (use warnings() to see them)
summary(warnings())
## 8264 identical warnings:
## In qt(1 - 10^(-4 + ncp), df = 35, ncp = ncp) :
##   full precision may not have been achieved in 'pnt{final}'

## MM: better than "1 - <small>" is to use 'lower.tail=FALSE':
r1. <- vapply(nc1, function(ncp) qt(10^(-4+ncp), df=35, ncp=ncp, lower.tail=FALSE), 1.23)
summary(warnings())
matplot(nc1, cbind(r1,r1.), type="l")
## Zoom in "left tail":  ==>  {lower.tail=FALSE} did not help much
i <- nc1 <= -10
matplot(nc1[i], cbind(r1,r1.)[i,], type="l")

## Zoom in around -7
(r2 <- sapply(seq(-6.9, -7.9, -0.1), function(ncp) qt(1-1*(10^(-4+ncp)), 35, ncp)))
 ## [1] -0.2268386        Inf        Inf        Inf -0.4857400 -0.5497784
 ## [7] -0.6135402 -0.6770143 -0.7401974 -0.8030853 -0.8656810
summary(warnings())
## 3086 identical warnings:
## In qt(1 - 1 * (10^(-4 + ncp)), 35, ncp) :
##   full precision may not have been achieved in 'pnt{final}'

nc2 <- seq(-60, -55, by=1/32)/8 # {1/8, 1/32: exact decimals}
r2 <- vapply(nc2, function(ncp) qt(1 - 10^(-4+ncp), df=35, ncp=ncp), 1.23)
## There were 52444 warnings (use warnings() to see them)
summary(warnings())
## 52444 identical warnings: ....
## lower.tail=FALSE
r2. <- vapply(nc2, function(ncp)
              qt(10^(-4+ncp), df=35, ncp=ncp, lower.tail=FALSE), 1.23)
all.equal(r2., r2) ## "Mean relative difference: 6.740991e-06"
matplot(nc2, cbind(r2,r2.), type="l")
## not why such a locally linear function should not work...

## In case it matters:

## > sessionInfo()
## R Under development (unstable) (2022-06-05 r82452)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Linux From Scratch r11.0-165

##--> MM: exploring
q2 <- seq(-3/4, -1/4, by=1/128)
pq2 <- pt(q2, 35, ncp=-7, lower.tail=FALSE)
                                        # no warnings here!
### ==> Providing qtU() a simple uniroot() - based inversion of pt()
qpqU  <- qtU(pq2, 35, ncp=-7, lower.tail=FALSE, tol=1e-10)
stopifnot(all.equal(q2, qpqU, tol=1e-9)) # perfect!

## and qtAppr() is really poor in comparison
plot(pq2 ~ q2, type="l", col=2, lwd=2)
qpq2  <- qt    (pq2, 35, ncp=-7, lower.tail=FALSE) # 6262 warnings
qAp2a <- qtAppr(pq2, 35, ncp=-7, lower.tail=FALSE, method="a")
qAp2b <- qtAppr(pq2, 35, ncp=-7, lower.tail=FALSE, method="b")
qAp2c <- qtAppr(pq2, 35, ncp=-7, lower.tail=FALSE, method="c")
summary(warnings())
## 6262 identical warnings: ...
## plot error -- interestingly *all* positive: .. then JUMP to Inf
plot (qpq2-q2 ~ q2, type="l", col=2, lwd=2); abline(h=0,lty=2)

## but the approximations are *very* bad, compared:
matplot(q2, cbind(qpq2, qAp2a,qAp2b,qAp2c)-q2,
        type="l", col=2:5, lwd=2)
## leave away the "b" which is particularly bad (here):
matplot(q2, cbind(qpq2, qAp2a,qAp2c)-q2,
        type="l", col=2:5, lwd=2)
## absolute values & log-scale:
matplot(q2, abs( cbind(qpq2, qAp2a,qAp2c) - q2), log="y",
        type="l", col=2:5, lwd=2)
## Idea: what about *starting* with qtAppr() and then do  1 or a few *Newton* steps?

## unfinished ...
mat2 <- cbind(q2, pq2, qpq2, Dq=qpq2-q2, relE=relErrV(q2, qpq2))
tail(mat2, 10)


###=== Large df ; R currently uses "primitive" (<==> qtNappr(*, k=0)) for df > 1e20
##     ========
##--- *CENTRAL* case  (ncp == 0) ----------

if(!dev.interactive(orNone=TRUE)) { dev.off(); pdf("t-nonc_P-3.pdf") }

## originally from ../man/qtAppr.Rd :
qts <- function(p, df, lower.tail = TRUE, log.p = FALSE, tolUa = 1e-14) {
    cbind(qt   = qt     (p, df=df, lower.tail=lower.tail, log.p=log.p)
        , qtU  = qtU    (p, df=df, lower.tail=lower.tail, log.p=log.p)
        , qtUa = qtU    (p, df=df, lower.tail=lower.tail, log.p=log.p, tol = tolUa)
        , qtN0 = qtNappr(p, df=df, lower.tail=lower.tail, log.p=log.p, k=0)
        , qtN1 = qtNappr(p, df=df, lower.tail=lower.tail, log.p=log.p, k=1)
        , qtN2 = qtNappr(p, df=df, lower.tail=lower.tail, log.p=log.p, k=2)
        , qtN3 = qtNappr(p, df=df, lower.tail=lower.tail, log.p=log.p, k=3)
        , qtN4 = qtNappr(p, df=df, lower.tail=lower.tail, log.p=log.p, k=4)
          )
}

p.qtNappr <- function(qtsR, p, ip, errMax = 1e-17, do.diff=FALSE) {
    stopifnot(is.matrix(qtsR), (n <- length(p)) == nrow(qtsR), 1 <= ip, ip <= n)
    (iN <- startsWith(prefix="qtN", colnames(qtsR))) # also "drop" qtU()'s result
    (dqN <- (qtsR[,iN] - qtsR[,1])[ip,])
    if(do.diff) ## difference
      matplot(p[ip], dqN, type="l", col=2:6,
              main = paste0("difference qtNappr(p,df) - qt(p,df), df=",df), xlab=quote(p))
    ## |difference|
    matplot(p[ip], pmax(abs(dqN), errMax), log="y", type="l", col=2:6,
            main = paste0("abs. difference |qtNappr(p,df) - qt(p,df)|, df=",df), xlab=quote(p))
    legend("bottomright", paste0("k=",0:4), col=2:6, lty=1:5, bty="n")
    ## Rel.err.
    matplot(p[ip], pmax(abs(dqN/qtsR[ip,"qt"]), errMax), log="y", type="l", col=2:6,
            main = sprintf("rel.error qtNappr(p, df=%g, k=*)",df), xlab=quote(p))
    legend("left", paste0("k=",0:4), col=2:6, lty=1:5, bty="n")
    invisible(dqN)
}

p <- (0:100)/100
ii <- 2:100 # drop p=0 & p=1  where q*(p, .) == +/- Inf

df <- 100 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
qsp1c <- qts(p, df = df)
matplot(p[ii], qsp1c[ii,], type="l") # "all on top"
p.qtNappr(qsp1c, p, ip=ii, do.diff=TRUE)

df <- 2000 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
qsp1c <- qts(p, df=df)
p.qtNappr(qsp1c, p, ip=ii)
## qtAppr is really worse even than k=2 here (df=2000):
qtA <- qtAppr(p, df=df, ncp=0)
lines(p, pmax(1e-17, abs(qtA/qsp1c[,1] - 1)),  col=2, lwd=2, lty=2)
legend(0.6, 1e-9, "qtAppr()",  bty="n",cex=.8, col=2, lwd=2, lty=2)
## qtA is between k=1 and k=2 ...
print(digits=4, cbind(p, cbind(qsp1c[,-1], qtA) - qsp1c[,1])[ii,])

df <- 1e8 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< large <<<<<<<<<<<<
## more p-values close to {0, 1}
p. <- 2^-(2*(15:5)); p. <- sort(c(p., p, 1 - rev(p.))); plot(diff(p.))
qspL <- qts(p., df=df)
ii <- which(is.finite(qspL[,"qt"]))
p.qtNappr(qspL, p=p., ip=ii)
lines(p., pmax(1e-17, abs(qspL[,"qtU"] /qspL[,"qt"] - 1)),  col=2, lwd=2, lty=2)
legend(0.6, 1e-9, "qtU(tol=1e-5)",  bty="n",cex=.8, col=2, lwd=2, lty=2)
lines(p., pmax(1e-17, abs(qspL[,"qtUa"]/qspL[,"qt"] - 1)),  col=4, lwd=2, lty=3)
legend(0.6, 1e-14, "qtU(tol=1e-14)",  bty="n",cex=.8, col=4, lwd=2, lty=3)

print(qspL, digits=6)
