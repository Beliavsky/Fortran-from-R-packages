#### NON-central [dpq]chisq()
#### ===========      =======  (Since *central [dpq]Chisq  <===> [dpq]gamma)

library(DPQ)
### originally was  ~/R/MM/NUMERICS/dpq-functions/chisq-nonc-ex.R
###                 - - - - - - - - - - - - - - - - - - - - - - -
### with _long_ history since R 0.63.x in 1999 !

###__FIXME__ Already have 3 different ./wienergerm*.R files
###  =====   Do remove things we have twice !!!

stopifnot(exprs = {
    require(graphics)
    require(sfsmisc) # eaxis(), lseq(), p.m(), mult.fig(), sessionInfoX()
})

source(system.file(package="DPQ", "test-tools.R",
                   mustWork=TRUE))# ../inst/test-tools.R
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()

unlist(.Platform)
## For package-testing "diagnostics":
sessionInfoX(c("DPQ","Rmpfr"))

(noLdbl <- (.Machine$sizeof.longdouble <= 8)) ## TRUE when --disable-long-double

## very large ncp gave "infinite" loop in R <= 3.6.1 :
## ==> need new enough "3.6.1 patched" or R{-devel} > 3.6.x
## (okR_Lrg <- (getRversion() >  "3.6.1" ||
##              getRversion() == "3.6.1" && R.version$`svn rev` >= 77145))

(doExtras <- DPQ:::doExtras() && !grepl("valgrind", R.home()))

## save directory (to read from):
(sdir <- system.file("safe", package="DPQ"))

## on "my" platform, and if doExtras,  I'm very strict:
(myPlatf <- all(Sys.info()[c("sysname", "machine", "login")] ==
                           c("Linux",  "x86_64",  "maechler")))
(beStrict <- doExtras && !noLdbl && myPlatf)
(is32 <- .Machine$sizeof.pointer == 4) ## <- should work uniformly on Linux/MacOS/Windows

if(!dev.interactive(orNone=TRUE)) pdf("chisq-nonc-1.pdf")
.O.P. <- par(no.readonly=TRUE)
showProc.time()

### Part 1 : Densities  dchisq(*, ncp)
### ----------------------------------

### densities alone :
## ===> shows Normal limit (for lambda -> Inf;  true also for nu -> Inf)
nu <- 12
nS <- length(ncSet <- if(doExtras) 10^seq(1, 9, by=1/2) else 10^(1:6))
np <- if(doExtras) 201 else 64
cpUse <- numeric(nS); names(cpUse) <- formatC(ncSet)
mult.fig(nS, marP = -0.8,
         main = paste0("non-central chisq(*, df=",nu,") and normal approx"))$old.par -> op
for(NC in ncSet) {
    m <- NC + nu
    s <- sqrt(2*(nu + 2*NC))
    x <- seq(from= m - 3*s, to= m + 3*s, length.out = np)
    cpUse[formatC(NC)] <- system.time(y <- dchisq(x, df=nu, ncp=NC))[1]
    plot(x, y, ylim=c(0,max(y)),type = "l", ylab='f(x)', main=paste("ncp =",NC))
    lines(x, dnorm(x, mean=m, sd=s), col = 'blue')
}
par(op)# resetting mult.fig()
showProc.time()

cbind(ncSet, cpUse, "c/ncp"= cpUse / ncSet)
## fails on Win 32b: "need finite 'ylim' values" :
try(plot(cpUse ~ ncSet, log = "xy", type = 'b', col = 2))
if(doExtras) try(# fails occasionally (too many zeros -> log(0) = -Inf)
  print(summary(lmll <- lm(log(cpUse) ~ log(ncSet), subset = cpUse > 0)))
)
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)
## (Intercept) -9.099690   0.100188  -90.83 1.20e-10 ***
## log(ncSet)   0.494316   0.005548   89.09 1.35e-10 ***
##
## Residual standard error: 0.08279 on 6 degrees of freedom
## Multiple R-Squared: 0.9992,	Adjusted R-squared: 0.9991  <<- !
## F-statistic:  7938 on 1 and 6 DF,  p-value: 1.347e-10

## =>   log(cpUse) ~= -9.1 + 0.494*log(ncSet)
## <==>  cpUse  proportional to  sqrt(ncp)
##       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## further experimenting shows, that only values in the center of the density
## take longer times!  ==> exactly where the normal approx (or others!) is good

###--- Now, limit for  nu = df --> Inf :
ncp <- 16
nS <- length(dfSet <- 10^(if(doExtras) 0:11 else 0:8))
cpUse <- numeric(nS); names(cpUse) <- formatC(dfSet)
oPar <- mult.fig(nS, marP = -0.75, main = "non-central chisq(ncp = 16) and normal approx")$old.par
for(DF in dfSet) {
    m <- DF + ncp
    s <- sqrt(2*(DF + 2*ncp))
    x <- seq(from= m - 3*s, to= m + 3*s, length = np)
    cpUse[formatC(DF)] <- system.time(y <- dchisq(x, df=DF, ncp=ncp))[1]
    plot(x, y, ylim=c(0,max(y)),type = "l", ylab='f(x)', main=paste("df =",DF))
    lines(x, dnorm(x,m=m,s=s), col = 'blue', lty=2)
}
par(oPar)
cbind(dfSet, cpUse, "c/df"= cpUse / dfSet)
## remains fast!
showProc.time()

## source("~/R/MM/NUMERICS/dpq-functions/dnchisq-fn.R")# dnoncentchisq() etc

## R
curve(dnoncentchisq(x, df=3, ncp=0), 0, 10)
curve(dchisq       (x, df=3),        0, 10, add=TRUE, col='purple')
## *no* visible difference -- ok
curve(dnoncentchisq(x, df=3, ncp=1), 0, 10)
curve(dchisq       (x, df=3, ncp=1), 0, 10, add=TRUE, col='purple') # ditto

x  <- seq(0, 10, length=101)
del <- c(0:4,10,40)
res <- matrix(NA, nr=length(x), nc=length(del),
              dimnames=list(NULL, paste0("ncp=",del)))
for(id in seq(along=del))
    res[,id] <- dnoncentchisq(x=x, df=3, ncp=del[id])

matplot(x, res)
title("dnoncentchisq(*, df=3, ncp = ..)  & dchisq(..)")
l.pch <- as.character(seq_along(del))
legend("topright", paste("ncp =", del), col=1:6, pch=l.pch)
res2 <- outer(x, del, function(x,del)dchisq(x=x, 3, ncp=del))
matplot(x, res2, add=TRUE) # practically no difference visible!

signif( cbind(x, abs(1 - res/res2)) , digits=4)

matplot(x, abs(1 - res/res2)[,-1], type="b", lty=1, log="y", yaxt="n"); eaxis(2)
title("Rel.Err |1 - dnoncentchisq(*, df=3, ncp = ..) / dchisq(..)|")
legend("bottomright", paste("ncp =", del[-1]), col=1:6, lty=1, pch=l.pch, bty="n")

showProc.time()

###---- March 2008 -----  "large ncp" :

n <- if(doExtras) 1e4 else 512

## From: Martin Maechler <maechler@stat.math.ethz.ch>
## To: Peter Dalgaard <p.dalgaard@biostat.ku.dk>
## Subject: Re: non-central chisq density
## Date: Thu, 27 Mar 2008 22:22:15 +0100

## [...............]

## Hi Peter,

## I've recently looked at problems in computing the non-central
## beta density for largish non-centrality parameter,
## which made me eventually considering the non-central chisq,
## and I have recalled your R News (Vol.1, nr.1; 2001)
## article on its big improvement by finding the maximal term in
## the sum and then sum outwards;
## and of course, the same idea will be applicable to the dnbeta()
## as well.

## However, for bigger ncp things, become eventually unfeasible as
## you've also noted in your article.
## Hoever, there, you've mentioned that the new method would work
## even up to ncp = 100'000^2  which I found astonishing (but not wrong)
## because I saw much potential for underflow already in the
## central term.
## Also, I saw that the current pnchisq() does not compute the
## log-density more accurately in places where the density
## underflows to zero...

    curve(dchisq(x, df=3, ncp=30000, log=TRUE), 0, 50000)

## also not a big deal, but for me one more reason to look
## at normal / central approximation formula for "large" ncp ..
## for all (or many of) the non-central distributions.

## Anyway, I started to look a bit more closely and then saw this

  d <- 1e6; curve(dchisq(x, df=3, ncp=d, log=TRUE), .98*d, 1.02*d, n=n)

## when going to smaller ncp and looking more closely something like

   curve(dchisq(x, df=3, ncp=30000, log=TRUE), 27300, 27500, n=n)

## which you may find amusing....

## all this no need for immediate action, but something to
## consider, and as said,
## I'm looking into this first for dnbeta(), but then more
## generally.

## All in all, your R News article has been again a very nice piece
## of inspiration.
##
## Martin

##  dchisqAsym (x, df, ncp, log = FALSE)  --> ../R/dnchisq-fn.R
##  ----------                                     ~~~~~~~~~~~~

curve(dchisq(x, df=3, ncp=30000, log=TRUE), 26000, 34000, n=n)
curve(dchisqAsym(x, df=3, ncp=30000, log=TRUE),
      add=TRUE, col="purple", n=n)
curve(dnorm(x, m=3+30000, sd=sqrt(2*(3 + 2*30000)), log=TRUE),
      add = TRUE, col="blue", n=n)
##==> It seems the chisqAsym() approximation is slightly better;
##  also from this :
x <- rchisq(if(doExtras) 1e6 else 1e4, df=3, ncp=30000)
(sN <- sum(dnorm(x, m=3+30000, sd=sqrt(2*(3 + 2*30000)), log=TRUE)))
(sCh<- sum(dchisqAsym(x, df=3, ncp=30000, log=TRUE))) ## larger (less negative) <-> better
all.equal(sN, sCh) # ... 2.6887e-6" [Win 32b: 2.873e-5]

## dnchisqBessel(x, df, ncp, log = FALSE) --> ../R/dnchisq-fn.R
## -------------                                 ~~~~~~~~~~~~

## From ?pl2curves()  [ == ../man/pl2curves.Rd ] :
p.dnchiB <- function(df, ncp, log=FALSE, from=0, to = 2*ncp, p.log="", n = if(doExtras) 2001 else 512, ...)
{
    pl2curves(dnchisqBessel, dchisq, df=df, ncp=ncp, log=log,
              from=from, to=to, p.log=p.log, n=n, ...)
}

## simple check
dchR  <- dchisq       (1:30, df = 3, ncp = 1:30)
dchBI <- dnchisqBessel(1:30, df = 3, ncp = 1:30)
          all.equal(dchR, dchBI, tol = 0) # "Mean rel.diff.: 2.3378e-14"
stopifnot(all.equal(dchR, dchBI, tolerance = 1e-13))

p.dnchiB(df=1.2, ncp=500,, 200, 800)
p.dnchiB(df=1.2, ncp=500, log=TRUE)# differ in tail

p.dnchiB(df=20, ncp=500,, 200, 800)
p.dnchiB(df=20, ncp=500, log=TRUE) # ok (differ for large x)
p.dnchiB(df=20, ncp=100) # looks good
p.dnchiB(df=20, ncp=100, , 0, 500, p.log="y") # looks good too (differ large x)
p.dnchiB(df=20, ncp=100, log=TRUE, 0,500)     # the same

p.dnchiB(df=20, ncp=200, log=TRUE, 0,600)     # the same
p.dnchiB(df=35, ncp=400, log=TRUE, 0,1500)    # the same
p.dnchiB(df= 3, ncp=600, log=TRUE, 0,2500)    # the same
p.dnchiB(df= 3, ncp=800, log=TRUE, 0,3500)    # for large x --> NaN in besselI

## Originally, large  ncp  gave overflow in besselI():
dnchisqBessel(8000, df=20, ncp=5000) # 1.3197e-78
       dchisq(8000, df=20, ncp=5000) # 6.8563e-79

## This is all fine now :
x <- if(doExtras) 1000:1600 else seq(1000, 1600, by = 5)
plot (x, besselI(x, 9,   TRUE), type="l")
lines(x, besselI(x, 1.2, TRUE), col=2)
lines(x, besselI(x, 1.0, TRUE), col=2)
lines(x, besselI(x, 0.1, TRUE), col=2)
lines(x, besselI(x, 1.8, TRUE), col=2)

### OTOH: Bessel asymptotic  I_a(y) ~  exp(y) / sqrt(2*pi*y)  for y >> a
lines(x, 1/sqrt(2*pi*x),  col=3, lty=3, lwd=3)
## hmm, looks like the  nu=1.2  case, but *not* the  nu=9  one ??

lines(x, besselI(x, 2.2, TRUE), col="blue")
lines(x, besselI(x, 3.2, TRUE), col="blue")
lines(x, besselI(x, 4.2, TRUE), col="blue")
lines(x, besselI(x, 5.2, TRUE), col="blue")
lines(x, besselI(x, 6.2, TRUE), col="blue")
lines(x, besselI(x, 7.2, TRUE), col="blue")
lines(x, besselI(x, 8.2, TRUE), col="blue")
##--> Need asymptotic for besselI(x, nu) with a term that depends on nu

##--> ~/R/MM/NUMERICS/Bessel/bessel-large-x.R  and better  ~/R/Pkgs/Bessel/
##    ~~~~~~~~~~~~~~~~~~~[April 2008]  ================

showProc.time()

### Part 2 :  pchisq (non-central!)
### -------------------------------

if(!dev.interactive(orNone=TRUE)) { dev.off(); pdf("chisq-nonc-2.pdf") }

## source("/u/maechler/R/MM/NUMERICS/dpq-functions/pnchisq.R")#-> pnchisq(), pnchisqV()

## In examples ../man/pnchisqAppr.Rd ---------
## ((again there at beginning))


### Note Brian's change (which completely broke df=0 case !) for R 2.3.0:

## r37287 | ripley | 2006-02-07 23:12:38 +0100 (Tue, 07 Feb 2006) | 6 lines

## improvements to [pq]nchisq
## - use direct formula which allows for lower_tail = FALSE if ncp < 80
##   (this is often a lot more accurate).
## - use starting point and lower_tail in qnchisq
## - can be slower, so make interruptible

## --- df = 0 ------------
stopifnot(pchisq(0:10, 0,1) >= exp(-1/2)) ## gave NaN from 2.3.0 to 2.6.1
## For a series of ncp = lambda :
lam <- seq(0,100, by=.25)
p00  <- pchisq(0,     df=0, ncp=lam)
p.0 <- pchisq(1e-300, df=0, ncp=lam)
stopifnot(all.equal(p00, exp(-lam/2), tol=2e-16),# '0' (when compiled alike)
          all.equal(p.0, exp(-lam/2), tol=4e-16))# was "1e-100" aka tol=0 ..

###------
### Accuracy buglet(s) :
## df -> 0 :  pnchisq() allows this, but it's pretty wierd :
## -------
## (and S-plus does it better !!)
## Theory: (df=0, ncp=0 ) is  point mass   (1)     at 0 --> 1-P == 0 everywhere
## ------  (df=0, ncp= L) has point mass exp(-L/2) at 0
plot(function(x)pchisq(x, df=0,    ncp=0, lower=FALSE),1e-1, 5000,log="x")## fine (all 0)
plot(function(x)pchisq(x, df=0,    ncp=0, lower=FALSE),2000, 5000) ## all = 0
plot(function(x)pchisq(x, df=1e-4, ncp=0, lower=FALSE),2000, 5000) ## all = 0 still
## this is ok (2014-04)
plot(function(x)pchisq(x, df=1e-4, ncp=0, lower=FALSE),1e-1, 5000, log="xy")## !?
## The R version of this:
curve(   pnchisqV     (x, df=1e-4, ncp=0, lower=FALSE), add=TRUE, col=adjustcolor(2,1/2), lwd=4)
## central chisq is ok here:
curve( pchisq(x, df=1e-4,        lower=FALSE),add = TRUE, col = "red")
curve( pchisq(x, df=1e-4,        lower=FALSE),1e-1, 5000)#, add = TRUE)
curve( pchisq(x, df=1e-4,        lower=FALSE),1e-1, 5000, log = 'xy')

##--- but the problem persists for df > 0 for small non-zero  ncp:
curve( pchisq(x, df=0.01, ncp = 0.1), 1e-1, 5000, log="x") # fine
par(new=TRUE)
curve( pchisq(x, df=0.01, ncp = 0.1, lower=FALSE,log=TRUE),
      1e-1, 5000, log="x", ylab="", yaxt="n", col=2)
axis(4, col.axis=2); mtext("log(1 - p)", 4, col=2)
## --> underflows to -Inf [because it computes log(.)

###--- this was "noncentral-ex.R" :

x <- x10 <- 10^(-300:300)#-> this is x-range is NOT plottable!
x <- x10 <- 10^(-150:150)#-> *is* plottable

system.time(pch.x10  <- pchisq(x,x ))# 0.01 in R;  4.77 in Splus 3.4
system.time(pch.x10n <- pchisq(x,x,ncp=1e-10))#-- hangs for ever [R <= 0.63.3]
## R 1.2.x : 0.57
## in S-plus 3.4:
##- Error in .C("S_ncchisq_prob",: subroutine S_ncchisq_prob:
##- 	284 Inf value(s) in argument 2
##- Dumped
##- Timing stopped at: 4.77 0.00999999 5 0 0

stopifnot(is.na(pch.x10) == is.na(pch.x10n))#> TRUE  R & Splus 3.4
stopifnot(!any(is.na(pch.x10)))

summary(pch.x10 - pch.x10n)
## Splus:
##-  Min. 1st Qu. Median Mean 3rd Qu. Max. NA's
##-     0       0      0    0       0    0  284

## R 1.2.x:
##-       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
##- -9.054e-09  5.000e-11  5.000e-11  2.426e-01  5.000e-01  5.000e-01

## Much better now:  Max =  5e-11

## However, closer inspection reveals
summary(pch.x10[pch.x10 != 1])
##-R     Min. 1st Qu. Median    Mean 3rd Qu.    Max.
##-R  0.5000  0.5000  0.5000  0.5271  0.5000  1.0000

##-S+3.4       Min. 1st Qu. Median  Mean 3rd Qu. Max.
##-S+3.4      -Inf    -Inf   -Inf   -Inf    -Inf    1
summary(pch.x10[pch.x10 != 1 & is.finite(pch.x10)])
##-    Min. 1st Qu. Median   Mean 3rd Qu. Max.
##-  0.4912  0.5009 0.9294 0.7696       1    1


###----- less extreme x values:

rele.pch <- function(x) 1 - pchisq(x,x,ncp=1e-100) / pchisq(x,x)

(rl <- rele.pch(2^(-12:10)))
stopifnot(rl == 0)## << S-plus 3.4, later R

## The uniroot no longer works, but look at
curve(1-pchisq(x,1,1), 69, 1500, n = 2001,
      main="1-pchisq(x,1,1), x in [69, 1500]", sub = R.version.string)
axis(1, at= c(69, 1500))
abline(h=0, col="gray70")
##--> lots of noise -- but around correct value = 1
## now, all fine :
pc1 <-
    curve(pchisq(x,1,1, lower.tail=FALSE), 69, 1500, n = 2001, log = "y",
      main="1-pchisq(x,1,1), x in [69, 1500]", sub = R.version.string)

## 1-P of course jumps to 0 *much* earlier:
plot(pc1, type="l", log="y", yaxt="n", ylab="",
      main="1-pchisq(x,1,1), x in [69, 1500]", sub = R.version.string)
eaxis(2) ## whereas this underflows much much much earlier (at x ~ 100)
curve(1-pchisq(x,1,1), add=TRUE, col=adjustcolor("red", 0.5), lwd=3, n = 2001)


x <- 100:1511
p <- pchisq(x,1,1, lower=FALSE)
stopifnot(0 <= p, p <= 2e-19)
lp <- log(p)
stopifnot(is.finite(lp), lp <= -43, abs(diff(lp) + 0.475) < 0.02)
showProc.time()

### try other (df,ncp) -- compare with  Wienergerm approx.
## setwd("/u/maechler/R/MM/NUMERICS/dpq-functions/")
## source("wienergerm_nchisq-fn.R")
## dyn.load("wienergerm_nchisq.so")

p.pchUp <- function(df, ncp, from, to, log.p=FALSE, n=2001) {
    c1 <-
        curve(pchisq(x, df=df, ncp=ncp, lower.tail=FALSE, log.p=log.p),
              from=from, to=to, n=n,
              main = paste0("pchisq(x, ",formatC(df),", ",formatC(ncp),
                            ", lower=F", if(log.p)", log=T", "),",
                            "  x in [",formatC(from),", ", formatC(to), "]"),
              sub = R.version.string)
    axis(1, at= c(from, to))
    abline(h=0, col="gray70", lty=3)
    c2 <- curve(pchisqW(x, df=df, ncp=ncp, lower.tail=FALSE, log.p=log.p), n=n,
                add = TRUE, col = 2)
    legend("topright", xjust = 1/2, yjust = 1.1,
           c("normal", "Wienerg.approx"), col=1:2, lty=1)
    invisible(data.frame(x = c1$x, pchisq = c1$y, pchisqW = c2$y))
}
## in all these, Wienerg.approx. looks very good: x >> df+ncp
p.pchUp( 0.1,  0.1,     51, 1500)
p.pchUp(   1,    1,     69, 1500)
p.pchUp(   1,    1,     69, 1500, log=TRUE)# small discrepancy for largish x
p.pchUp( 2.5, 0.02,     61, 1500, log=TRUE)# (nothing visible)
p.pchUp(  25,   10,    150, 1600, log=TRUE)# discrepancy largish x
summary(warnings()) ; showProc.time()
p.pchUp( 100,  500,    980, 1900)# normal:noise(and cutoff at 1); Wienerg: smooth
p.pchUp( 100,  500,    980, 1900, log=TRUE)# pchisq() breaks down
p.pchUp( 500,  100,    897, 3000)
p.pchUp( 500,  100,    897, 3000, log=TRUE)# pchisq break down
p.pchUp( 5e3,  100,   5795, 1.e4)  # zoom in..
p.pchUp( 5e3,  100,   5777, 6000)  # --> Wiener has less noise long before!
## (but it also is systematically a bit larger - correct?)
summary(warnings()) ; showProc.time()

if(doExtras) withAutoprint({ # -----------------------------------
## Now have m + 5*s cutoff, ...
cc <- p.pchUp( 5e3,  5e3,  10400, 11e3) # still pchisq() jumps to 0 at 10866.2, too early
p.m(cc, type="l", log="y", lwd=c(1,3), col=c("black", adjustcolor("red",0.5)))
## now (larger cutoff) "fine" but also too early to jump to zero:
c2 <- p.pchUp( 1e5,  2e4, 11.6e4, 12.6e4, n=1001) # see Wienergerm-singularity !!
p.m(c2, type="l", log="y", lwd=c(1,3), col=c("black", adjustcolor("red",0.5)))

## Still shows that the   m + 5*s cutoff is too early!
p.pchUp( 5e3,  5e3,  10800, 11e3)
cc <- p.pchUp( 5e3,  5e3,  8000,  20e3)
p.m(cc, type="l", log="y", lwd=c(1,3), col=c("black", adjustcolor("red",0.5)))
p.pchUp( 1e5,  2e4, 12.25e4, 12.35e4)# m + 5*s  __much__ too early here..
showProc.time() # ~ 0.5 sec
}) ## only if(doExtras) -----------------------------------

### NOTA BENE: We have the *big* problem when s ~= 1,  x <= ncp+df
### ---------  ---------------------------------------------------
### this is unsolved and the reason we have not yet ...

### ==> conclusion:  Use Wienergerm as soon as P > 1 - 1e-12,
##      ----------   but probably much earlier
##  To *find* that  P > 1 - 1e-12  we could try a cheap  qchisq() approx.
##  unfortunately, these are quite INaccurate in the tails...


## when pnchisq.c is  RE-compiled with  -DDEBUG_pnch
## these give interesting output

## Simulate this, using pnchisq()
## Ok, "now" for  ncp <= 80,  we use direct formula
## "now" := r37287 | ripley | 2006-02-07 23:12:38
##
## ---> these no longer use old algo:

##                   Case  lt     n(#it)
pnchisq(1000, 1,1, verbose=2)#   2   -496.8   666
pnchisq(1300, 1,1)#   2   -646.6   838
pnchisq(1400, 1,1)#   2   -696.6   895
pnchisq(1420, 1,1)#   2   -706.6   906
pnchisq(1422, 1,1)#   2   -707.6   907

pnchisq(1425, 1,1)#   2   -709.6 L + x large --> 1
pnchisq(1430, 1,1)#   2   -711.6 L + x large --> 1
pnchisq(1490, 1,1)#   2   -741.6 L + x large --> 1



## With the newly (2003-01-16) visible warning [no longer; 2004-02]
pchisq(1e-5, df=100, ncp=1)
## [1] 0
##- Warning message:
##- too large x (= 1e-05) or noncentrality parameter 1 for current algorithm;
##- 	result is probably invalid!

pnchisq(1e-5, df=100, ncp=1, verbose = TRUE)
##  lt= -758.8, t=exp(lt)= 0
##- _OL_: n= 1 fx.2n = 102  > 0 ==> flag
## [1] 0

## where as
pnchisq(10e-5, df=100, ncp=1, verbose = TRUE)
##  lt= -643.7, t=exp(lt)= 2.92e-280
## _OL_: n= 1 fx.2n = 102  > 0 ==> flag

## [1] 1.771154e-280


##---------------- another bug: large x  with ncp > 0

### NOTA BENE:  Fix this with "Wiener Germ Approx." formula !!!

### now (R 1.8.x at least) ok
mult.fig(3, main = "pchisq(x >= 1497,*, ncp=)  BUG (no longer!)")$old.par -> op
curve(pchisq(x, df=1, ncp=  1), from=1,to=1e4, log='x', main="ncp = 1")
curve(pchisq(x, df=1, ncp=300), from=1,to=1e4, log='x', main="ncp = 300")
curve(pchisq(x, df=1, ncp=0), from=1,to=1e4, log='x', main="ncp = 0")
par(op)

## still (2004-01 ... 2014-01 !!) true:
## now looking closer, at the upper tail {algorithm is not good on log scale!}
curve(pchisq(x, df=1, ncp=0, lower=FALSE,log=TRUE),
      from=1,to=1e4, log='x', main="ncp = 0")# -> goes down to -700 or so

## ncp > 80 is different ..
xp <- curve(pchisq(x, df=1, ncp=300, lower=FALSE,log=TRUE), xaxt="n",
            from=1, to=1e4, log='x', main="ncp = 300, log=TRUE")# only down to ~ -25
eaxis(1, sub10=2)
## .. hmm, really bad...

## .. the reason is that we compute on (lower=TRUE, log=FALSE) scale and only then transform:
## --> gives warnings! (and 'verbose' output):
curve(pchisq(x, df=1, ncp=300, lower=FALSE),
      from=100,to=2000, log='xy', main="ncp = 300,  upper tail", axes=FALSE) -> pxy
summary(warnings())
eaxis(1, sub10=3); eaxis(2)
curve(pnchisqV(x, df=1, ncp=300, errmax = 4e-16, lower=FALSE, verbose=1),# ,log=TRUE),
      add = TRUE, col=2); mtext("ncp = 300 -- pnchisqV() pure R", col=2)
showProc.time()


## also seems to hang (or take much too long?) on Winbuilder [32 bit *and* 64 bit ]
if(.Platform$OS.type == "unix" && !noLdbl) withAutoprint({
pncRC <- pnchisqRC(pxy$x, df=1, ncp=300, lower=FALSE, verbose=1)
all.equal(pxy$y, pncRC, tol = 0)# "often" TRUE, depends on exact R version, etc
stopifnot(
    all.equal(pxy$y, pncRC, tol = if(noLdbl) 5e-14 else 0)# noLdbl: seen 1.129e-14
)
summary(warnings())
showProc.time()
})# ---------------------only if(.. "unix" ....)----------------------------


## Really large 'df' and 'x' -- "case I":
## no problem anymore:
f <- c(.9,.999,.99999, 1, 1.00001,1.111, 1.1)
x <- 1e18*f
stopifnot(exprs = {
    all.equal(pchisq(x, df=1e18, ncp=1) -> p,
              c(0,0,0, 1/2, 1,1,1))
    all.equal(p, pnchisqRC(x, df=1e18, ncp=1), tol = 4e-16) # see 0
})
## case I -- underflow protection large x --> 1
tt <- 10^-(6:12)
stopifnot(!is.unsorted(xm <- 1e18*(1 + c(-tt, 0, rev(tt)))))
(pn <- pnchisqV (xm, df=1e18, ncp=1)) #-> 0...1 is correct
pp  <- pchisq   (xm, df=1e18, ncp=1)
##
if(.Platform$OS.type == "unix") withAutoprint({ #-------------------
pp. <- pnchisqRC(xm, df=1e18, ncp=1, verbose=1)
## Pnchisq_R(x, f, th, ... lower.tail=1, log.p=0, cut_ncp=80, it_simple=110,
##   errmax=1e-12, reltol=1.77636e-15, epsS=8.88178e-16, itrmax=1000000, verbose=1)
##   --> n:= max(length(.),..) = 15
## but then does *NOT* terminate in time on Winbuilder
all(pp == pp.)# >>> TRUE :  *RC is also C code, perfect
all.equal(pp, pn, tol = 0) # see 1.6e-16  (even on Win 32b)
if(doExtras && !noLdbl)  # who knows ..
    stopifnot(pp == pp.)
stopifnot(exprs = {
    all.equal(pp, pp., tol = 1e-15) # see 0
    all.equal(pp, pn,  tol = 1e-15) # see 1.6e-16
})
})## only if( .. unix .. )

## (also "problematic" with Wienergerm: s=0)
showProc.time()#-----------------


### largish f and ncp   {no problem visible here, but see below}
curve(pchisq(x, df= 10000, ncp=300),
      from=1e-3, to=20000, log='x', main="ncp = 300")
curve(pchisq(x, df= 10000, ncp=300),
      from=2000, to=20000, log='x', main="ncp = 300")

x <- seq(3000,11000, length=201)
if(FALSE)## to see the break:
x <- seq(5500,11000, length=201)
px <- pchisq(x, df=10000, ncp=300, log=TRUE)
plot(x, px, type='l', col=2, lwd=2,
     main="pchisq(*, df=10000,ncp=300, log=TRUE))")
head(px, 20)
## for the (5500,11000): -Inf -Inf ..... -Inf -650.2379 -640.3743 -630.6..
showProc.time()

pnchisq(5500, df= 10000, ncp=300, verbose=2)
##  lt= -744.4 => exp(lt) underflow protection ln(x)= 8.612503
## _OL_: n= 1  n= 1, fx.2n = 4502 > 0
## BREAK n= 1 ; bound=  0 <= errmax rel.err=  NaN <= reltol

## New: allow large ncp {this now (that we use reltol !) *TAKES TIME*}
curve(pnchisqV(x, df= 10000, ncp= 4000),
      from=12000, to= 16000, main="df = 10000, ncp = 4000")
curve(400*dchisq(x, df= 10000, ncp= 4000), add = TRUE, col = "green")
## R's pchisq() looks fine now:
curve(pchisq(x, df= 10000, ncp= 4000), add=TRUE,
      col=adjustcolor("blue",1/4), lwd=4)

curve(pnchisqV(x, df= 16e3, ncp= 16e3),
      from=30e3, to= 35e3, main="df = 16e3, ncp = 16e3")
curve(400*dchisq(x, df= 16e3, ncp= 16e3), add = TRUE,
      col = adjustcolor("green4",.5), lwd=3)
showProc.time()

if(doExtras) withAutoprint({
## current R version: -- (also relatively slow, but much faster!) *and* non-convergence warning
rr <- curve(pchisq(x, df= 10000, ncp=3e5), type = "o", cex = 1/2,  n = 49)
summary(warnings()) ## all non-convergences (but *looks* ok)
## non-convergence in 100000 iterations : -- S.L.O.W. (~ 1 min. on 2014 lynne)
rV <- curve(pnchisqV(x, df= 10000, ncp=3e5), n = 49,
            from=3.13e5, to= 3.14e5, main="ncp = 3e5 - pnchisqV()")
summary(warnings())
identical(rr$x, rV$x)
showProc.time()#-----------------
})

### NOTA BENE:  dnchisq() has a similar sum and the following  i.Max
imaxD <- function(x,df,lambda)
    pmax(0, ceiling(0.25* (-(2+df) +sqrt((df-2)^2 + 4*lambda* x))))
## A few test comparisons with  ssR4[,,] :
## ==> imaxD() is too small (has correct order of magnitude, unless for
##     p(x) > 1-1e-4

### Investigate pnchisq() algorithm  sum(v[i] * t[i]) :
###
### rr(i) : is it increasing / decreasing / maximal / .. ?
plRpois(lambda = 1000)

mult.fig(8, main = r_pois_expr, tit.wid = 6)$old.par -> op
for(la in c(5,20,50,100,200,500,1500,5000))
    plRpois(lambda = la, do.main=FALSE)
par(op)
## -> Wow!
## 1)  Always decreasing;
## 2)  r(i,lambda) < lambda / i  and ~=~ for i < lambda)

## How well approximation from  *ratio* point of view:
## Interesting i < ~ lambda  (ok, clearly improvable):
## ----------- i >= lambda  : need better approx
plotrq <- function(lambda, i = 1:(3*round(lambda))) {
    lab <- as.expression(substitute(lambda==la, list(la=lambda)))
    plot(i, r_pois(i,lambda=lambda) / (lambda/i),
         type='b', cex=.4, col=2)
    abline(v=lambda, col='gray', lty=2)
    axis(3, at = lambda, labels = lab)
}
plotrq(10)

mult.fig(4, main = " r_pois(i) /  (lambda/i)")$old.par -> op
plotrq(20)
plotrq(50)
plotrq(100)
plotrq(500)
par(op)
showProc.time()

## How well approximation from  difference point of view:
## Interesting i < ~ lambda  (ok, clearly improvable):
## ----------- i >= lambda  : need better approx
plotDr <- function(lambda, i = 1:(4*round(lambda))) {
    lab <- as.expression(substitute(lambda==la, list(la=lambda)))
    plot(i,  (lambda/i) - r_pois(i,lam=lambda),
         type='b', cex=.4, col=2)
    abline(v=lambda, col='gray', lty=2)
    axis(3, at = lambda, label = lab)
}

mult.fig(9, main = quote(plotDr(lambda)))$old.par -> op
plotDr( 4)
plotDr(10)
plotDr(20)
plotDr(50)
plotDr(100)
plotDr(200)
plotDr(500)
plotDr(1000)
plotDr(2000)## oops: problem (no longer !)
par(op)

### Now back to the original problem:
### Using ss() terms and see where they are maximal, etc.
(pR <-           pnchisq   (1.2,df=1,ncp=3, verbose=FALSE))# iter = 12, now 13
all.equal(c(pR), pnchisq_ss(1.2,df=1,ncp=3), tol=0)# 2.19e-12, now 9.61e-14,
## 6.4e-16 on Win 32b !

(pR  <- pnchisq   (1.2, df=1, ncp=30, verbose=FALSE))# iter = 12, now 16
(pss <- pnchisq_ss(1.2, df=1, ncp=30))
all.equal(pR, pss, tol = 0)
## Now you can also pass an ss() result to pnchisq_ss():
str(ss1.2 <- ss(1.2, df=1, ncp=30))
stopifnot(identical(pss,
                    pss. <- pnchisq_ss(1.2, df=1, ncp=30, ssr = ss1.2)))
## was  2.616 e-8 (thanks to 'reltol'!)
(pR <-        pnchisq   (1.2,df=1,ncp=30, verbose=FALSE, reltol=3e-16))# 19 it.
all.equal(pR, pnchisq_ss(1.2,df=1,ncp=30), tol= 2e-16)

str(sss <- ss(1.2,df=1,ncp=30))# s[1:161], max = 3
plot(sss$s, type="h", col=2)
## i: for log-ax bug (warning) {still not nice looking
range(which(i <- 1e8*sss$s > .Machine$double.xmin * sss$s[sss$max]))# 1:160
plot(sss$s[i], type="b", col=2, log = 'xy')
## Which indices are relevant for the sum?
range(which(ii <- sss$s > .Machine$double.eps * sss$s[sss$max]))
## 1:19 -- as we had 19 iterations above!
stopifnot(sum(sss$s[ii]) == sum(sss$s))

## Left tail probabilities are now much better:
(pR <- pnchisq (1.2, df=100, ncp=30, verbose=FALSE,reltol=3e-16))
## 5.384254 e-83 , 12 iter.
       pchisq  (1.2, df=100, ncp=30)
## 4.461632 e-83, which is identical to
       pnchisq (1.2, df=100, ncp=30, reltol=1)# =^= "old" C code (1 iter!)

### What about large df and x -- #{terms} ?
str(sss <- ss(100,100, 1e-3))# 1 469
pnchisq_ss(100,100,1e-3)
pchisq    (100,100,1e-3)
((Ss <- sum(sss$s)) - sum(rev(sss$s)))/Ss # -1.9286 e-16

ss2(100,100, 1e-3)
##-  i1  i2 iN1 iN2 max
##-   1 469   1  71   1
Ns <- 2^c(-200, -15, -5, -1:15, 30, 100)
names(Ns) <- paste("2",formatC(log(Ns,2)),sep="^")
tab.ss1c <- t(sapply(Ns, function(u) ss2(100,100,ncp=u, i.max=10000)))
tab.ss1c
##-> i2 is "constant": 469 (or 468);problems from ncp >= 2^12 = 4096
tab.ss10 <- t(sapply(Ns, function(u) ss2(10,10, ncp=u, i.max=10000)))
cbind(tab.ss10, tab.ss1c) ## only from ncp ~= 2^6, things change

(t1k.1c <- t(sapply(Ns, function(u) ss2(1000,100, ncp=u))))
## even with i.max = 100000,  thing "go wrong" from ncp = 2^11
str(s.. <- ss(1000,10, 2048))
s..$s[1:400] #-- sequence diverges to +Inf -- can we better re-scale?
## (yes, we can: pnchisq() does so -- leave this for now)

## Now vary x from small to large:
(t.x.1k <- t(sapply(Ns, function(x) ss2(x,df=100, ncp=100))))# probl. from 2^11


str(s <- ss(1000,100, ncp=3000))
str(s <- ss(100,100, ncp=1000))
##  $ s  : num [1:469] 1.00e+00 4.91e+02 1.18e+05 1.86e+07 2.16e+09 ...
##  $ i1 : int 1
##  $ max: int 136
s$s[s$max] #  1.4 e-130 : down scaled
ss2(100,100, ncp=1000)
##-  i1  i2 iN1 iN2 max
##-   1 468  68 216 136
ss2(100,100, ncp=2000)
##-  i1  i2 iN1 iN2 max
##-  95 326 118 296 201

## But:
all( ss(100,100,5000)$s == 0) # TRUE -- no longer
## because  lu needs much better scaling "-lambda" is too much
## "Fixed" :
table( ss(100,100, ncp=5000)$s ) ## only values in {0, Inf}, mostly Inf !

##==> give up for these high ncp for the moment!

showProc.time()

## Instead use C - code which parallels  pnchisq()'s in C:
## dyn.load("/u/maechler/R/MM/NUMERICS/dpq-functions/pnchisq-it.so")
str(pit <- pnchisqIT(3,2,4))# 1:21
stopifnot(with(pit, all.equal(sum(terms), prob)))
##  this is a bit funny: all 0 terms
stopifnot(with(pit2 <- pnchisqIT(100,100,5000),
               all.equal(sum(terms), prob)))
all(pit2$terms == 0)# TRUE
stopifnot(with(pit3 <- pnchisqIT(100,100,5),
               all.equal(sum(terms), prob, tol=1e-15)))
str(pit3)# 1:69
str(pit <- pnchisqIT(10000,10000,5))#  567
str(pit <- pnchisqIT(10004,10000,5))#  569 terms
str(pit <- pnchisqIT(10010,10000,5))#  572 terms (i0=0)
str(pit <- pnchisqIT(12000,10000,5))# 1612 terms (i0=0)
## hmm, quite interesting:
plot(pit$terms,type='l')
par(new=TRUE)
plot(pit$terms,type='l', log = 'y',yaxt='n',col=2)# looks like -x^2 !
axis(4, col.axis=2)
summary(pit$terms) # max =  0.005150251 -- the first few 100 are unneeded

str(pit <- pnchisqIT(12000,10000,5000))# 2442 terms, i0=877; max.term= 2.5e-60 !
## many unneeded terms!
str(pit <- pnchisqIT(15000,10000,5000))# 3189 terms, i0=877; max=.003287

str(pit <- pnchisqIT(20000,10000,5))# -> 1 immediately {0 terms}

## Now use ss2.() for the "term statistics":
ss2.(15000,10000, 5000)# 3189 terms, i0=877
ss2.(1,    10000, 5000)# immediate 0 -> 1 term only
ss2.(1e5,  10000, 5000)# immediate 1 -> "0 terms"

## Takes (already quite a bit) time:
rs <- sapply(14990:15010, function(x) ss2.(x,10000,5000))
t(rs)
## as expected : n.terms gives proper 'right border':
stopifnot(rs["i2",] == rs["nT", ],
          rs["i2",] == rs["iN2", ])
## swap df & ncp ===> the *double* number of terms!
x <- c(1000,10000,14000,14500, 14990,15000,15010,15500, 15600, 15670, 15675)
rs <- sapply(x, function(x) ss2.(x,5000,10000))
cbind(t(rs), prob = sapply(x, function(x) pnchisqIT(x, 5000,10000)$prob))
##   i0   nT   i1   i2  iN1  iN2 iMax          prob
##    1    1   NA   NA    1    1    1  0.000000e+00
## 2596 4298 2596 4298 3499 4298 3907 1.697298e-137
## 2596 5290 2880 5290 4358 5290 4810  2.625129e-06
## 2596 5469 2953 5469 4459 5469 4921  1.212588e-02
## 2596 5684 3031 5684 4560 5684 5044  4.838253e-01
## 2596 5689 3032 5689 4562 5689 5047  5.016652e-01
## 2596 5694 3034 5694 4564 5694 5049  5.194951e-01
## 2596 5942 3116 5942 4675 5942 5251  9.867808e-01
## 2596 5995 3134 5995 4701 5995 5301  9.960697e-01
## 2596 6031 3146 6031 4719 6031 5336  9.984810e-01
## 2596 6034 3147 6034 4721 6034 5338  9.985853e-01 << was "1.00" !

## but we cannot go too far :
str(pp <- pnchisqIT(20000, 5000,10000))
## $ prob   : num 1
## $ i0     : int 3995
## $ n.terms: int 8283
## $ terms  : num [1:8283] 0 0 0 0 0 0 0 0 0 0 ...
1-pp$pr; 1-sum(sort(pp$terms))
## -8.999912e-12
## -8.999024e-12
##  i.e.  P > 1 is certainly not okay anymore!
## E     = df + ncp              = 15000
## sigma = sqrt( 2*(df + 2*ncp)) =   223.607
5000 / sqrt( 2*(5000 + 2*10000))
## = 22.36 i.e.  20'000 is 22.3 sigma out of E[]

showProc.time()

set.seed(635)
## 1st simulation ---------------------------------------------------------------
## Collect data: --- this took about 2 hours on "nb-mm" (P III, 700 MHz)
##                   takes 4.5 sec on ada-17 (or alo
nL <- 20
nF <- 16
nX <- length(pX <- c(0.01, (1:9)/10, 0.99, 0.9999))

sfil1 <- file.path(sdir, "tests_chisq-nonc-ssR.rds")
if(!doExtras && file.exists(sfil1)) {
  ssR_l <- readRDS_(sfil1)
  str(ssR_l)
  ## dfs :  num [1:16] 15.9 20.7 21 29.5 47.8 ...
  ## lam :  num [1:20] 5.74 8.26 8.34 8.64 10.12 ...
  ## ssR :  num [1:4, 1:20, 1:16, 1:12] 8.08 1 25 2 9.24 ...
  loadList(ssR_l)

} else { ## do run the simulation always  if(doExtras) :

lam <- sort(rlnorm(nL, 3, 1))
dfs <- sort(rlnorm(nF, 4, 1))
ssR <- array(NA, dim=c(1+3, nL,nF,nX),
             dimnames = list(c("x","iN1","iN2", "iMax"),
                             lam = formatC(lam,digits=5),
                             df  = formatC(dfs,digits=5),
                             x   = formatC(pX, width=1)))
for(iL in 1:nL) {
    lm <- lam[iL]
    cat("lam=", formatC(lm),":")
    for(iF in 1:nF) {
        f <- dfs[iF]
        x <- qchisq(pX, df=f, ncp=lm)
        for(iX in 1:nX)
            ssR[, iL,iF,iX] <- c(x[iX], ss2.(x[iX], df=f, ncp=lm)[5:7])
        cat(".")
    }; cat("\n")
}
save2RDS(list_(lam, dfs, ssR), file = sfil1)

} # {run simulation}

x. <-  ssR["x"   ,,,]
iM <-  ssR["iMax",,,]
iN1 <- ssR["iN1" ,,,]
iN2 <- ssR["iN2" ,,,]
iS <- iN2 - iN1 # the "index Spread": how many terms need to be summed
## Visualize iM(x) for some (df,lambda):
Sel <- function(i) round(quantile(i, names=FALSE))

mult.fig(mfrow=c(5,5), ## << since length(quantile()) == 5
         marP = c(-1,-1,0,0))$old.par -> op
for(iL in Sel(1:nL)) {
    lm <- lam[iL]
    for(iF in Sel(1:nF)) {
        f <- dfs[iF]
        plot(x.[iL,iF,],
             iM[iL,iF,], type = 'o', xlab = "x", ylab = "iMax",
             main=paste("df=",formatC(f),", lam=",formatC(lm)))
    }
}
par(op)
##--> 1st order, qualitatively "same" function (x)

## Same plot, but using "Wienergerm's"  sW(x,df,lam) instead of x
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/wienergerm_nchisq-fn.R")
mult.fig(mfrow=c(5,5), ## << since length(quantile()) == 5
         marP = c(-1,-1,0,0), main = "iMax vs.  sW()")$old.par -> op
for(iL in Sel(1:nL)) {
    lm <- lam[iL]
    for(iF in Sel(1:nF)) {
        f <- dfs[iF]
        plot(sW(x.[iL,iF,], df=f, ncp=lm)$s,
             iM[iL,iF,], type = 'o', xlab = "sW(x,...)", ylab = "iMax",
             main=paste("df=",formatC(f),", lam=",formatC(lm)))
    }
}
par(op)
## very similar
showProc.time()

###--- visualize 'iN1' = the first index *needed* :
## Idea: use "current" algorithm (simply summing from i=1...) when ok :

fCont <- function(ix = 6, kind = c("iN1","iN2","iMax", "d.i"),
                  pch=1, cex=.5,
                  sdat = ssR,
                  lam = as.numeric(dimnames(sdat)[["lam"]]),
                  dfs = as.numeric(dimnames(sdat)[["df"]]),
                  ## what a horrible hack ..
                  pX  = as.numeric(dimnames(sdat)[["x"]])
                  )
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 18 Feb 2004, 18:31
    kind <- match.arg(kind)
    datname <- deparse(substitute(sdat))
    nx <- dim(sdat)[4]
    if(ix < 1 || ix > nx) stop("'ix' must be in 1:",nx)
    if(kind == "d.i") {
        m <- sdat["iN2" ,,, ix] - sdat["iN1" ,,, ix]
        mtxt <- paste("Spread ", datname,"[ 'iN2 - iN1' ,,, ", ix,"]", sep='')
    } else {
        m <- sdat[kind ,,, ix]
        mtxt <- paste(datname,"[",kind," ,,, ", ix,"]", sep='')
    }
    mtxt <- paste(mtxt, " (i.e., x=",
                  formatC(100*pX[ix],digits=10,wid=1),"%-perc.)",sep='')
    if(kind == "iN1") {
        filled.contour(lam, dfs, m,
                       levels=c(1,2,3,5,10,15,20,30)-.01,
                       col = c("light gray", terrain.colors(6)),
                       plot.axes={ points(expand.grid(lam,dfs),cex=cex,pch=pch)
                                   axis(1); axis(2)},
                       plot.title={ title(mtxt,
                                          xlab="ncp (lambda)",ylab="df (nu)")}
                       )
    } else {                            # automatic levels and color
        filled.contour(lam, dfs, m,
                       color.palette = terrain.colors,
                       plot.axes={ points(expand.grid(lam,dfs),cex=cex,pch=pch)
                                   axis(1); axis(2)},
                       plot.title={ title(mtxt,
                                          xlab="ncp (lambda)",ylab="df (nu)")}
                       )
    }
}

par(.O.P.)# just in case for filled.contour() to work
fCont()# iN1, 6
fCont(1)# iN1, 11
fCont(kind="d", pch='.')##  "spread"
fCont(kind="iN2", cex=.25)## practically == "spread" (i.e. iN1 ~= 0 !)
fCont(12, kind="iN2")

showProc.time()

##
## 4th simulation: ----------------------------------------------------
##
## pX: at these quantiles(), compute pnchisq()
nx <- length(pX <- c(0.01, (1:9)/10, 0.99, 0.9999, 1-1e-6, 1-1e-9))

sfil4 <- file.path(sdir, "tests_chisq-nonc-ssR4.rds")
if(!doExtras && file.exists(sfil4)) {
    ssR_l <- readRDS_(sfil4)
    str(ssR_l)
    loadList(ssR_l)
} else {
    set.seed(41)
    ## smallish (lam,df) -- use several x --- ideally all these give "iN1"=1
    ss <- c(.1, .2, .5, 1:5, 7,10,15,c(2:6,8,10,12,15,20,30,50)*10)
    nl <- length(lam4 <- ss)
    nd <- length(dfs4 <- ss)
    ## I use "round" numbers, and can pack all the info into ssR4's dimnames:
    ssR4 <- array(NA, dim=c(1+ 3, nl,nd,nx),
                  dimnames = list(c("x", "iN1","iN2", "iMax"),
                  lam = formatC(lam4),
                  df  = formatC(dfs4),
                  pr.x= ifelse(pX < 0.999, formatC(pX),
                               paste("1-", formatC(1-pX),sep=''))))
    for(il in 1:nl) {
        lm <- lam4[il]
        for(id in 1:nd) {
            f  <- dfs4[id]
            ## use more than one x per (df,lam) pair:
            x <- qchisq(pX, df=f, ncp=lm)
            for(ix in 1:nx)
                ssR4[, il,id,ix] <- c(x[ix], ss2.(x[ix], df=f, ncp=lm)[5:7])
            cat(".")
        }; cat(il,"")
    }; cat("\n")

    save2RDS(list_(lam4, dfs4, ssR4), file=sfil4)
}

## Compute the 'iMax' value that corresponds to x = E[X] = df + ncp
## from that, with the above 'general f(x)', we can hopefully
## get a good estimated   iMax(x, df, ncp) ...
str(iM.E <- iM[,,1])
for(iL in 1:nL) {
    lm <- lam[iL]
    for(iF in 1:nF) {
        f <- dfs[iF]
        E <- f + lm # = E[X]
        iM.E[iL,iF] <- approx(x.[iL,iF,],
                              iM[iL,iF,], xout = E)$y
    }
}

str(iM.E)



persp(lam, dfs, iM.E)# pretty linear..
##->
filled.contour(lam, dfs, iM.E, xlab="ncp lambda", ylab="df nu")
## indeed: very linear in "lambda / ncp", a bit less linear in "df:
dat1 <- cbind(expand.grid(lam,dfs), c(iM.E))
names(dat1) <- c("ncp", "df", "iM")
summary(lm1 <- lm(iM ~ ncp + df, data = dat1)) ## R^2 = 99.79%
with(dat1, p.res.2x(x=ncp,y=df, residuals(lm1))) # pretty structured..
image   (x=lam,y= dfs, array(residuals(lm1),dim=dim(iM.E)))

summary(lm2 <- lm(iM ~ ncp * df, data = dat1)) ## R^2 = 99.84%
summary(lm3 <- lm(iM ~ poly(ncp, df, degree=2), data = dat1)) ## R^2 = 99.96%

## Using sqrt() -- not really better (and predict/fitted is failing !? FIXME !!
summary(lm3s <- lm(sqrt(iM) ~ poly(ncp, df, degree=2), data = dat1))

with(dat1, p.res.2x(x=ncp,y=df, residuals(lm3))) # pretty structured..
image   (x=lam,y= dfs, array(residuals(lm3),dim=dim(iM.E)))

library(mgcv)
summary(gam1. <- gam(iM ~ s(ncp) + s(df), data = dat1))
## df = 1 + 4 + 4; R^2 = 0.999; s^ = 0.226
summary(gam1.2 <- gam(iM ~ ncp + s(df), data = dat1))
## df = 2 + 4    ; R^2 = 0.999; s^ = 0.280
plot(gam1.2) #pretty square

summary(gam2. <- gam(iM ~ s(ncp, df), data = dat1)) ## 100% explained,
## but equiv.deg.freedom = 1+26.3
if(FALSE)## FAILS
summary(gam2.10 <- gam(iM ~ s(ncp, df, 10), data = dat1))
## df = 1 + 9    ; R^2 = 1.000;  s^ = 0.107
with(dat1, p.res.2x(x=ncp,y=df, residuals(gam2.)))
## ok, but
plot(fitted(gam2.), fitted(lm3)) ; abline(0,1,col=3)
## suggesting the quadratic fit being quite ok.

## OTOH, I do need an explicit formula;
## a simple 2d- regression spline instead of quadratic?

showProc.time()

###---- 2nd "simulation": -- only go for one x = E[]
nSim <- if(doExtras) 5000 else 500
sfil2 <- file.path(sdir, "tests_chisq-nonc-ssR2.rds")
if(!doExtras && file.exists(sfil2)) {
    ssR2 <- readRDS_(sfil2)
} else {
    set.seed(2)
    lam <- rlnorm(nSim, 5, 2)
    dfs <- rlnorm(nSim, 6, 1.5)
    ssR2 <- rbind(rbind(lam,dfs), matrix(NA, 3, nSim))
    for(i in 1:nSim) {
        lm <- lam[i]
        f  <- dfs[i]
        x <- f + lm
        ssR2[3:5, i] <- ss2.(x, df=f, ncp=lm)[5:7]
        cat("."); if(i %% 100 == 0) cat("\n",i)
    }
    dimnames(ssR2) <- list(c("lam","df","iN1","iN2", "iMax"),NULL)
    ssR2 <- t(ssR2)
    ssR2 <- ssR2[sort.list(ssR2[,"lam"]),]
    save2RDS(ssR2, file=sfil2)
}

## 3rd simulation: --- this takes a little while (1 min ?)
sfil3 <- file.path(sdir, "tests_chisq-nonc-ssR3.rds")
if(!doExtras && file.exists(sfil3)) {
    ssR3_l <- readRDS_(sfil3)
    str(ssR3_l)
    loadList(ssR3_l)
} else {
    set.seed(31)
    ss <- c(20,50, c(1:8,10,13,18,25)*100, 3000)
    nl <- length(lam3 <- c(ss, seq(5000, 100000, length=1+ 4*19)))
    nd <- length(dfs3 <- c(ss, seq(5000, 100000, length=1+ 2*19)))
    ssR3 <- array(NA, dim=c(3, nl,nd),
                  dimnames = list(c("iN1","iN2", "iMax"),
                  formatC(lam3), formatC(dfs3)))
    for(il in 1:nl) {
        lm <- lam3[il]
        for(id in 1:nd) {
            f  <- dfs3[id]
            x <- f + lm
            ssR3[, il,id] <- ss2.(x, df=f, ncp=lm)[5:7]
        }; cat(il,"")
    }; cat("\n")
    save2RDS(list_(lam3, dfs3, ssR3), file=sfil3)
}
showProc.time()

### now change these "3" values into a data.frame as this one:
dsR2 <- as.data.frame(ssR2)

##
d3 <- dim(ssR3)
ss3 <- matrix(ssR3, d3[1], prod(d3[2:3]))
dsR3 <- cbind(expand.grid(lam = lam3, df = dfs3), t(ss3))
colnames(dsR3)[3:5] <- dimnames(ssR3)[[1]]

dsR <- rbind(dsR2, dsR3)
rownames(dsR) <- paste(1:nrow(dsR))
## visualize "design space":
iOutl <- c(648, 1841, 5000)
plot(df ~ lam, data =dsR, log = "xy")
points(dsR[iOutl,], col=2:4, cex=2)
dsR[iOutl,]

plot(df ~ lam, data =dsR, log = "")
points(dsR[iOutl,], col=2:4, cex=2)
points(dsR[4997:5000,], col="gray", cex=2, pch=3)

with(dsR, which(lam > 100000))# 4998, 4999, 5000
## -- leave these away for the regression (high leverage!)
str(dsR. <- subset(dsR, lam <= 100000))

summary(l.1 <- lm(iMax ~ lam,      data=dsR.))## R^2(adj) = 1;  s^ = 23.73
summary(l.2 <- lm(iMax ~ lam + df, data=dsR.))##                s^ = 12.49
summary(l.3 <- lm(iMax ~ lam * df, data=dsR.))##                s^ = 12.22
summary(l.4 <- lm(iMax ~ lam * df+ I(lam^2), data=dsR.))##      s^ =  8.251
summary(l.5 <- lm(iMax ~ lam * df+ I(lam^2)+I(df^2), data=dsR.))#s^=  7.812
##               Estimate Std. Error  t value Pr(>|t|)
## (Intercept)  9.437e+00  1.105e-01    85.42  < 2e-16
## lam          5.022e-01  1.034e-05 48573.09  < 2e-16
## df           9.861e-04  1.095e-05    90.03  < 2e-16
## I(lam^2)    -1.333e-08  1.219e-10  -109.39  < 2e-16
## I(df^2)     -4.202e-09  1.257e-10   -33.44  < 2e-16
## lam:df       6.433e-10  9.405e-11     6.84 8.42e-12
summary(l.6 <- update(l.5, . ~ . + log(lam)))## R^2(adj) = 1;  s^ =  6.135 (7 p)
##               Estimate Std. Error  t value Pr(>|t|)
## (Intercept) -7.477e+00  2.348e-01   -31.85   <2e-16
## lam          5.016e-01  1.179e-05 42528.90   <2e-16
## df           8.951e-04  8.681e-06   103.11   <2e-16
## I(lam^2)    -8.579e-09  1.137e-10   -75.48   <2e-16
## I(df^2)     -3.804e-09  9.881e-11   -38.50   <2e-16
## log(lam)     3.391e+00  4.373e-02    77.54   <2e-16
## lam:df       1.546e-09  7.477e-11    20.68   <2e-16
summary(l.7 <- update(l.5, . ~ . + log(lam)*log(df)))##        s^ =  5.389 (9 p)
##-                    Estimate Std. Error   t value Pr(>|t|)
##- (Intercept)       6.315e+00  5.773e-01    10.938  < 2e-16 ***
##- lam               5.014e-01  1.066e-05 47049.885  < 2e-16 ***
##- df                4.992e-04  1.204e-05    41.449  < 2e-16 ***
##- I(lam^2)         -7.278e-09  1.034e-10   -70.356  < 2e-16 ***
##- I(df^2)          -4.404e-10  1.131e-10    -3.893 9.98e-05 ***
##- log(lam)          4.249e-02  8.164e-02     0.520   0.6028
##- log(df)          -2.062e+00  8.728e-02   -23.628  < 2e-16 ***
##- lam:df           -1.775e-10  7.703e-11    -2.305   0.0212 *
##- log(lam):log(df)  5.348e-01  1.151e-02    46.476  < 2e-16 ***

drop1(l.7)# cannot drop non-sign. log(lam)
summary(l.8 <- update(l.7,  . ~ . - log(lam)))##               s^ =  5.389 (8 p)
summary(l.9 <- update(l.8,  . ~ . - lam:df))  ##               s^ =  5.391 (7 p)
summary(l.10<- update(l.9,  . ~ . - I(df^2))) ##               s^ =  5.396 (6 p)
summary(l.11<- update(l.10, . ~ . - I(lam^2)))##               s^ =  5.396 (6 p)
summary(l.12<- update(l.11, . ~ . + log(lam)))##               s^ =  5.396 (6 p)

## dsR3 instead of dsR. :
summary(l.13<- lm(iMax ~ lam*df+ log(lam)*log(df), data=dsR3))
summary(l.14<- update(l.13, . ~ . - lam:df))
summary(l.15<- update(l.14, . ~ . - 1))

with(dsR3, p.res.2x(lam, df, residuals(l.15)))
## ok; it's really the low 'lam' (and the low 'df')
if(doExtras) { ##-> try more -----------------------------------
iMaxR3 <- matrix(dsR3$iMax, length(lam3))
persp         (log10(lam3), log10(dfs3), iMaxR3/ (lam3/2))
persp         (log10(lam3), log10(dfs3), log10(iMaxR3/ (lam3/2)))
filled.contour(log10(lam3), log10(dfs3), iMaxR3/ (lam3/2),
               plot.title = {
                   contour(log10(lam3), log10(dfs3), iMaxR3/ (lam3/2),add=TRUE)
                   title(main = "iMax / (lam/2)  ['dsR3' data]",
                         xlab = "log10(lam)", ylab = "log10(df)")
                   ##points(expand.grid(log10(lam3), log10(dfs3)), pch='.')
                   with(dsR3, points(log10(lam), log10(df), pch='.'))
               })
## almost same with  log(iMax / (lam/2)):
filled.contour(log10(lam3), log10(dfs3), log10(iMaxR3/ (lam3/2)),
               plot.title = {
                   contour(log10(lam3), log10(dfs3), log10(iMaxR3/ (lam3/2)),add=TRUE)
                   title(main = "log10(iMax / (lam/2))  ['dsR3' data]",
                         xlab = "log10(lam)", ylab = "log10(df)")
                   ##points(expand.grid(log10(lam3), log10(dfs3)), pch='.')
                   with(dsR3, points(log10(lam), log10(df), pch='.'))
               })
} #-- only if(doExtras) -------------------------------------

showProc.time()

if(doExtras && require("interp")) {
    ## same with both data --> need interp ! : s/dsR3/dsR./ :
ds1 <- subset(dsR., lam >= 1)
sr.I  <- with(ds1, interp(log(lam), log(df), iMax))
sr.Iq <- with(ds1, interp(log(lam), log(df), iMax / (lam/2)))
filled.contour(sr.I, xlab="ln(lam)", ylab="ln(df)", main="iMax")
filled.contour(sr.Iq,
               plot.title = {
                   contour(sr.Iq, add=TRUE)
                   title(main = "iMax / (lam/2)  ['ds1' data]",
                         xlab = "ln(lam)", ylab = "ln(df)")
                   with(ds1, points(log(lam), log(df), pch='.'))
               })
print(summary(l.q1 <- lm(iMax / (lam/2) ~ log(lam) * log(df), data= ds1)))
TA.plot(l.q1)
plot(resid(l.q1) ~ lam, data=ds1, pch ='.', log = 'x')
print(summary(l.q2 <- update(l.q1, .~. + I(1/lam))))
print(summary(l.q3 <- update(l.q2, .~. + I(log(lam)^2) + I(1/log(lam)))))
plot(resid(l.q3) ~ lam, data=ds1, pch ='.', log = 'x')
### --> Aha!   1/lam seems the best term !!
with(dsR., p.res.2x(lam, df, residuals(l.q3), scol=2:3))
## -- maybe try  lam^(-a)  ?
showProc.time() # 0.9
} # only if(.X.)

## This is impressive (but shows "non-fit")
with(dsR., p.res.2x(log10(lam), log10(df), residuals(l.5), scol=2:3))
with(dsR., p.res.2x(lam, df, residuals(l.5)))

with(dsR., p.res.2x(lam, df, residuals(l.10), scol=2:3))
with(dsR., p.res.2x(log(lam), log(df), residuals(l.6), scol=2:3))

plot(l.5, ask=FALSE) ## 2-3 outliers:
## 5000 : maximal lambda
## 1841 : maximal df

if(doExtras) withAutoprint({ # -----------------------------------
### Yet another idea:
summary(lq2 <- lm(I(iMax/lam) ~ (lam+ log(lam) + df + log(df))^2, data=dsR.))
lq2s <- step(lq2)
summary(lq2s, corr=TRUE, symb=TRUE)
## shows the complete non-sense (large lambda values fit very badly
with(dsR., n.plot(fitted(lq2s)*lam, iMax))

if(doExtras)## GAM -- needs tons of cpu + memory:
   summary(g.5 <- gam(iMax ~ s(lam) + s(df) + s(lam,df), data=dsR.))#s^=4.489
## -> (too) many deg.freedom s
}) #--------------------------------------------

showProc.time()


if(doExtras && require("interp")) { ## visualize more: ----------
sr2I <- with(dsR., interp(log(lam), log(df), iMax))
filled.contour(sr2I, xlab="ln(lam)", ylab="ln(df)", main="iMax")
sr2I <- with(dsR., interp(log(lam), log(df), log(iMax)))
sr2I <- with(dsR., interp(log(lam), log(df), log(iN2)))
filled.contour(sr2I, xlab="ln(lam)", ylab="ln(df)", main="ln(iN2)")
##  linear for large lambda
persp(sr2I,xlab="log(lam)", ylab="log(df)",zlab="log(iMax)",ticktype="detailed")

## restrict on those where iN1 > 1
str(dsR2r <- subset(dsR2, iN1 > 1))# only 3383 instead of 5000

sr2I <- with(dsR., interp((lam), (df), (iN2)))
filled.contour(sr2I, xlab="(lam)", ylab="(df)", main="(iN2)")
## Looks *very* nicely linear
persp(sr2I,xlab="(lam)", ylab="(df)",zlab="(iMax)",ticktype="detailed")

##
sr2I <- with(dsR., interp((lam), (df), iMax/lam))
filled.contour(sr2I, xlab="(lam)", ylab="(df)", main="iMax/lam")
persp(sr2I,xlab="(lam)", ylab="(df)",zlab="iMax/lam",ticktype="detailed")

## restrict on those where iN1 > 1
str(dsR.r <- subset(dsR., iN1 > 1))

sr2rI <- with(dsR.r, interp((lam), (df), (iN2)))
persp(sr2rI,xlab="(lam)",ylab="(df)",zlab="(iN2)",ticktype="detailed")
sr2rI <- with(dsR.r, interp(log(lam), log(df), log(iN2)))
persp(sr2rI,xlab="log(lam)",ylab="log(df)",zlab="log(iN2)",ticktype="detailed")

sr2rI <- with(dsR.r, interp(log(lam), log(df), log(iMax)))
persp(sr2rI,xlab="log(lam)",ylab="log(df)",zlab="log(iMax)",ticktype="detailed")
} else {
    cat("Define   dsR2r : \n") ; str(dsR2r <- subset(dsR., iN1 > 1))
    cat("and also dsR.r : \n") ; str(dsR.r <- subset(dsR., iN1 > 1))
}
showProc.time()
summary(ll.2 <- lm(log(iMax) ~ log(lam) + log(df), data=dsR.r))
summary(ll.3 <- lm(log(iMax) ~ log(lam) * log(df), data=dsR.r))
summary(ll.4 <- lm(log(iMax) ~ log(lam) * log(df) + I(log(lam)^2), data=dsR.r))
plot(residuals(ll.2) ~ dsR.r$lam, log='x')
plot(residuals(ll.3) ~ dsR.r$lam, log='x')
plot(residuals(ll.4) ~ dsR.r$lam, log='x')
plot(dsR.r$iMax - exp(fitted(ll.4)) ~ dsR.r$lam, log='x')

if(doExtras) {
summary(gl.4 <- gam(log(iMax) ~ s(lam) + log(df), data=dsR.r))## very bad
## but this is very good:
summary(gl.4 <- gam(log(iMax) ~ s(log(lam)) + log(df), data = dsR.r))
plot(gl.4, ask=FALSE)
if(FALSE) { # fails now
summary(gl.5 <- gam(log(iMax) ~ s(log(lam),4) + log(df)*log(lam), data=dsR.r))
plot(gl.5, ask=FALSE)
}
} # only if(.X.)
##-> try
summary(ll.5 <- lm(log(iMax) ~ (log(lam) + poly(pmax(0,log(lam)-5),2))*log(df),
                   data=dsR.r))

summary(dsR.r$iMax - exp(fitted(ll.5))) # one very negative
plot(ll.5, ask=FALSE)
showProc.time()


## First try to find formula for maximal number of terms needed
summary(l.N2.1 <- lm(iN2 ~ lam*df , data=dsR2r))
summary(l.N2.2 <- lm(iN2 ~ lam*df + I(lam^2)+I(df^2), data=dsR2r),corr=TRUE)
summary(l.N2.3 <- lm(iN2 ~ lam*df + I(lam^2)+I(lam^3), data=dsR2r),corr=TRUE)
summary(l.N2.4 <- lm(iN2 ~ lam*df + I(lam^2)+I(lam^3)+I(df^2), data=dsR2r),corr=TRUE)
plot(residuals(l.N2.2) ~ dsR2r$lam)

dsR2r$lamP20k <- pmax(0, dsR2r$lam - 20000)
dsR2r$lamM20k <- pmin(0, dsR2r$lam - 20000)
summary(l.N2.1P <- lm(iN2 ~ lam+df + I(lamP20k ^2)+I(lamM20k ^2) , data=dsR2r))


## This is to save typing: all variables 'log'ged:
dLr <- as.data.frame(lapply(dsR2r, log))
summary(l.N2. <- lm(iN2 ~ lam*df + I(lam^2)+I(df^2), data=dLr),corr=TRUE)
summary(l.N2  <- lm(iN2 ~ lam*df + I(lam^2)*I(df^2), data=dLr))
## back transformed residuals:
r <- dsR2r$iN2 - exp(fitted(l.N2))
n.plot(dsR2r$lam, r, log='x'); abline(h=0,lty=3)
## extreme (negative): 3383, also 3381, 3382
n.plot(dsR2r$lam, r, log='x', ylim=500*c(-1,1)); abline(h=0,lty=3)

showProc.time()


###---- older tests  ---------------------------------------


if(dev.interactive(TRUE)) { cat("sys.function(): "); str(sys.function()) }
if(.do.ask <- dev.interactive() && is.null(sys.function()))
    par(ask=TRUE); cat(".do.ask : ", .do.ask, "\n")
mult.fig(2)$old.par -> op
## large NC -- still (2018-08) very expensive!!
## 10^(3:10)  is (still!)  much too expensive, 10^8 alone costs 31.8 sec !
NC2 <- if(doExtras) 10^(2:7) else 10^(2:6)
for(NC in NC2) {
    cat("ncp=",NC,":\n")
    curve(dchisq(x, df=1, ncp=NC), from=NC/10,to=NC*100,
          log='x', main=paste("Density ncp =",NC))
    curve(pchisq(x, df=1, ncp=NC), from=NC/10,to=NC*100,
          log='x', main=paste("CDF    ncp =",NC))
    showProc.time()
}
par(op)
if(.do.ask) par(ask=FALSE)

## NOTE: I found that the median = qchisq(1/2, *)  is "mostly"
## ----  in (m-1, m) where m = mean = nu + lambda = df + ncp

## One exception (I've carefully searched for) where
## median < m-1  <==>   m - median > 1 : (maximal ~ 1.6):
df <- .0005; curve((df+x) - qchisq(1/2, df, ncp=x), 1, 2, col=2)
df <- .005 ; curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE, col=3)
df <- .05  ; curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE, col=4)

## These are all quite close (and quite CPU-costly !!) :
df <- 1e-4; curve((df+x) - qchisq(1/2, df, ncp=x), 0.5, 40, col=2)
abline(h=1, col='gray')
df <- 1e-4; curve((df+x) - qchisq(1/2, df, ncp=x), 1.2,   2, col=2)
if(doExtras) {
df <- 2e-4; curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE,col=3)
df <- 4e-4; curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE,col=4)
df <- 8e-4; curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE,col=5)
}
df <-16e-4; curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE,col=6)
showProc.time() # 0.41 {2019-09}

df <- 1e-100; curve((df+x) - qchisq(1/2, df, ncp=x), 1.38, 1.40, col=2)
dfs <- 2^(if(doExtras) seq(-300,-2, length=21) else -2) # as they are costly
for(df in dfs) {
    curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE,col=3)
    cat(formatC(df)," ")
}; cat("\n")
showProc.time()

if(doExtras) { ## -- show irregularity more closely ---------------
df <- 1e-300; curve((df+x) - qchisq(1/2, df, ncp=x), 1.38628, 1.38630, col=2)
for(df in dfs) {
    curve((df+x) - qchisq(1/2, df, ncp=x), add=TRUE,col=3)
    cat(formatC(df)," ")
}; cat("\n")
curve((0+x) - qchisq(1/2, df=0, ncp=x), 1.386294, 1.386295, col=2)
showProc.time() # doExtras: ~ 0.6 {2019-09}
} # only if(.X.)  -------------------------------------------------

ff <- function(ncp) (0+ncp)-qchisq(1/2, df=0, ncp=ncp)
str(oo <- optimize(ff, c(1.3,1.4), maximum=TRUE, tol=1e-15),digits=16)
##  $ maximum  : num 1.386294373354218
##  $ objective: num 1.386294355703814
qchisq(1/2, df=0, ncp = 1.386294373354218)## = 1.765e-8


## This is the case  df -> 0 where the distribution has a point mass at 0 !
x <- c(0,1e-8,1e-5,1e-3,seq(0.01, 3, length = if(doExtras) 1001 else 125))
plot (x, pchisq(x, df=1e-4, ncp = 1.4), type ='l', col=2, ylim = 0:1)
lines(x, pchisq(x, df=1e-4, ncp = 1.6), col=1)
lines(x, pchisq(x, df=1e-4, ncp = 1.2), col=3)
lines(x, pchisq(x, df=1e-4, ncp = 1.1), col=4)
lines(x, pchisq(x, df=1e-4, ncp = 0.1), col=5)

plot (x, pchisq(x, df=1e-2, ncp = 1.4), type ='l', col=2, ylim = 0:1)
lines(x, pchisq(x, df=1e-2, ncp = 1.6), col=1)
lines(x, pchisq(x, df=1e-2, ncp = 1.2), col=3)
lines(x, pchisq(x, df=1e-2, ncp = 1.1), col=4)
lines(x, pchisq(x, df=1e-2, ncp = 0.1), col=5)

plot (x, pchisq(x, df=0.1, ncp = 1.4), type ='l', col=2, ylim = 0:1)
lines(x, pchisq(x, df=0.1, ncp = 1.6), col=1)
lines(x, pchisq(x, df=0.1, ncp = 1.2), col=3)
lines(x, pchisq(x, df=0.1, ncp = 1.1), col=4)
lines(x, pchisq(x, df=0.1, ncp = 0.1), col=5)

showProc.time()

## MM: from something *not* put into ~/R/D/r-devel/R/tests/d-p-q-r-tests.R
## 1) PR#14216 (r51179, 2010-02-25)
x <- 80:200; lp <- pchisq(x, 4, ncp=1, log.p=TRUE)
stopifnot(is.finite(lp), all.equal(lp[1],-2.5519291e-14), lp < 0,
	  ## underflowed to 0, in R <= 2.10.x
	  -.4635 < (dll <- diff(log(-lp))), dll < -.4415,
	  max(abs(diff(dll))) < 3.75e-4)
##
showProc.time()


###---- again {may repeating from above --- "sorry I can't check that now"} :

x <- 250; pchisq(x, 1.01, ncp = 80, log=TRUE)

## R-2.10.0 and earlier --> quite a noisy picture ! --
## note that log P > 0  <==>  P > 1    --- of course nonsense!
xy <- curve(pchisq(x, 1.01, ncp = 80, log=TRUE), 250, 600,
            n=1001, ylim = c(-1,1)*8e-14); abline(h=0, lty=3)
## still noisy, and still slightly (5e-14) above 0 !

## bigger picture: for theta = ncp < 80, it works by using the other tail
plot(p1 <- -pchisq(1:400, 1.01, ncp = 80*(1-.Machine$double.eps), log=TRUE),
     log="y", type="o", pch=".")
## However, for  ncp >= 80 -- the  P = sum_i  term_i computation
lines(p2 <- -pchisq(1:400, 1.01, ncp = 80, log=TRUE),
      col=adjustcolor(2,0.5),lwd=2)## underflow to 0 .. not good ___ FIXME ___
## here, "the other tail",  log1p(- pnchisq(....., !lower_tail) )  does not work !!
summary(1 - p1/p2)

## From: Prof Brian Ripley <ripley@stats.ox.ac.uk>
## To: Martin Maechler <maechler@stat.math.ethz.ch>
## cc: R-core@r-project.org
## Subject: Re: p[n]chisq() warnings [was "d-p-q-r tests failures"]
## Date: Tue, 24 Nov 2009 15:14:16 +0000 (GMT)

## Martin,


## I mistyped 'mendacious": the error message lies.  I think it is
## generally wrongly worded, and should be something like 'full precision
## may not have been achieved'.

## Here is why I added the warning I added:

## MM: added 'lower.tail' 'log.p' and 'x'  arguments

##__ FIXME 2019-09:  Compare with my new  pnchis1sq()  function !!
t1 <- function(p, ncp, lower.tail = FALSE, log.p = FALSE,
               x = qchisq(p, df = 1, ncp, lower.tail=lower.tail, log.p=log.p))
{

    ## X ~ chi^2(df=1, ncp = L^2)  <==> X = Z^2 where Z ~ N(L, 1)
    ## --------------------------       -------       -----------
    ## P[ X > x ] = P[ |Z| > sqrt(x) ] = P[Z > sx] + P[Z < -sx] ,  sx := sqrt(x)

    p1 <- pchisq(x, df = 1, ncp = ncp, lower.tail=lower.tail, log.p=log.p)

    sx <- sqrt(x)
    sL <- sqrt(ncp)
    p2 <-
        if(!log.p) {
            if(lower.tail)
                pnorm(sx,  sL) -
                pnorm(sx, -sL, lower.tail=FALSE)
            else ## lower.tail = FALSE
                pnorm(sx,  sL, lower.tail=FALSE) +
                pnorm(sx, -sL, lower.tail=FALSE)
        } else { ## log scale -- MM: use logspace.add() and *.sub() for the above
            if(lower.tail)
                logspace.sub(pnorm(sx,  sL, log.p=TRUE),
                             pnorm(sx, -sL, log.p=TRUE, lower.tail=FALSE))
            else ## lower.tail = FALSE
                logspace.add(pnorm(sx,  sL, log.p=TRUE, lower.tail=FALSE),
                             pnorm(sx, -sL, log.p=TRUE, lower.tail=FALSE))
        }
    c(if(!missing(p)) c(p=p), x=x, pnchisq=p1, p.true=p2, relErr=abs(p1-p2)/p2)
}

t1(1e-12, 85)
## [1] 1.000000e-12 2.642192e+02 1.003355e-12 9.943394e-13 9.066654e-03
## Warning messages:
## 1: In qchisq(p, df = 1, ncp, lower.tail = FALSE) :
##    full precision was not achieved in 'qnchisq'
## 2: In pchisq(x, df = 1, ncp = ncp, lower.tail = FALSE) :
##    full precision was not achieved in 'pnchisq'

## so the answer is out by about 1%.  And

t1(1e-14, 100)
## [1] 1.000000e-14 5.208816e+02 0.000000e+00 6.107835e-38 1.000000e+00

## has lost all precision.  [MM: still true, Aug.2019]

## This sort of thing (because we compute 1 - answer) does not happen in
## the other tail.  So unless someone can show examples of precision
## loss, I believe that the warning in that tail should not be there (and
## would need conditional wording).

## MM: As soon as you go to log scale, completely inaccurate values around 1
##     are completely unuseful, too:
t1(x = 500, ncp=80, lower.tail=TRUE, log.p=TRUE)
##             x       pnchisq        p.true        relErr
##  5.000000e+02  3.552714e-15 -2.423206e-41 -1.466121e+26


## Brian
showProc.time()

## On Tue, 24 Nov 2009, Martin Maechler wrote:

## >>>>>> Prof Brian Ripley <ripley@stats.ox.ac.uk>
## >>>>>>     on Tue, 24 Nov 2009 12:22:48 +0000 (GMT) writes:
## >
## >    > On Tue, 24 Nov 2009, Peter Dalgaard wrote:
## >    >> Prof Brian Ripley wrote:
## >    >>> I only picked up that change this morning, and am seeing the failures
## >    >>> too.  I don't see why the warning is being given (isn't the test that
## >    >>> full accuracy was achieved?), so updating the .save file does not look
## >    >>> to me to be the solution.
## >    >>
## >    >> Hmm, I get the warnings, but it doesn't seem to stop the build for me
## >    >> and make check is failing at a different spoot:
## >
## >    > At an *earlier* spot in the check.
## >
## >    [.............]
## >
## >    >> The relevant diff is
## >    >> --- src/nmath/pnchisq.c (revision 50552)
## >    >> +++ src/nmath/pnchisq.c (revision 50553)
## >    >> @@ -40,9 +40,15 @@
## >    >> if (df < 0. || ncp < 0.) ML_ERR_return_NAN;
## >    >>
## >    >> ans = pnchisq_raw(x, df, ncp, 1e-12, 8*DBL_EPSILON, 1000000,
## >    >> lower_tail);
## >    >> -    if(!lower_tail && ncp >= 80) {
## >    >> -       if(ans < 1e-10) ML_ERROR(ME_PRECISION, "pnchisq");
## >    >> -       ans = fmax2(ans, 0.0);  /* Precaution PR#7099 */
## >    >> +    if(ncp >= 80) {
## >    >> +       if(lower_tail) {
## >    >> +           if(ans >= 1-1e-10) ML_ERROR(ME_PRECISION, "pnchisq");
## >    >> +           ans = fmin2(ans, 1.0);  /* e.g., pchisq(555, 1.01, ncp = 80) */
## >    >> +       }
## >    >> +       else { /* !lower_tail */
## >    >> +           if(ans < 1e-10) ML_ERROR(ME_PRECISION, "pnchisq");
## >    >> +           ans = fmax2(ans, 0.0);  /* Precaution PR#7099 */
## >    >> +       }
## >    >> }
## >    >> return log_p ? log(ans) : ans;
## >    >> }
## >    >>
## >    >> which warns if you get too close to 1.0 and truncates to 1.0 if you
## >    >> overshoot. All the cases tested should give the result 1.0 and thus
## >    >> trigger the warning. Are you implying that this is unintentional?
## >
## >    > I don't know nor can I guess Martin's intention, but I am confident
## >    > the warning is medacious here.
## >
## > Hmm, I don't understand "medacious".
## >
## > But anyway:  The new code of  `` pmin(ans, 1) '' is indeed necessary;
## > previously,  pchisq(x, df, ncp)  *would* return values larger
## > than one, ... somewhat embarrassingly.
## >
## > If you study a bit further, you'll find that currently,
## > pnchisq() for  ncp > 80  use identical code for TRUE or FALSE
## > lower_case;  and the old code
## > had a check for  ncp >= 80 and accuracy warnings for "upper tail"
## > and P < 1e-10.
## > The logical extension is to give the same accuracy warning for
## > "lower tail" and  P > 1 - 1e-10.

## > Of course, this is all just a workaround for the fact that our
## > current algorithm(s) are not good enough currently in those
## > extreme tail cases, and indeed,
## > I've start investigating better algorithms quite a while in the
## > past.
## > The creating of package 'Rmpfr' (for multi-precision arithmetic)
## > has BTW been influenced by my desire to get tools for exploring
## > such extreme tail misbehavior of current R algorithms.
## >
## > Here an example from one of my R scripts on this :
## >
## > ## R-2.10.0 and earlier --> quite a noisy picture ! --
## > ## note that log P > 0  <==>  P > 1    --- of course nonsense!
## > curve(pchisq(x, 1.01, ncp = 80, log=TRUE), 250, 600,
## >      n=1001, ylim = c(-1,1)*5e-14)
## >
## > So, again: these warning are a *substitute* and "cheap
## > workaround"  for now, but not
## > only for the new case that I've added, but also already for the
## > case Brian had added earlier:
## >    if(!lower_tail && ncp >= 80) {
## >       if(ans < 1e-10) ML_ERROR(ME_PRECISION, "pnchisq");
## >       ans = fmax2(ans, 0.0);  /* Precaution PR#7099 */
## >    }
## >
## > Martin
## >
## >
## >    > The save file in R-devel (which also gives the warnings) was updated
## >    > in r50552.
## >
## >    >>
## >    >> -p
## >    >>
## >    >>> Brian
## >    >>>
## >    >>> On Tue, 24 Nov 2009, Kurt Hornik wrote:
## >    >>>
## >    >>>>>>>>> Kurt Hornik writes:
## >    >>>>
## >    >>>>> I can no longer build r-patched.  Most likely from
## >    >>>>> r50553 | maechler | 2009-11-23 23:50:13 +0100 (Mon, 23 Nov 2009) | 1
## >    >>>>> line
## >    >>>>
## >    >>>>> ported r50552 [pchisq(*, ncp > 80) from trunk
## >    >>>>
## >    >>>>> I now get
## >    >>>>
## >>>>>> ##-- non central Chi^2 :
## >>>>>> xB <- c(2000,1e6,1e50,Inf)
## >>>>>> for(df in c(0.1, 1, 10))
## >    >>>>> +     for(ncp in c(0, 1, 10, 100)) stopifnot(pchisq(xB, df=df,
## >    >>>>> ncp=ncp) == 1)
## >    >>>>> There were 12 warnings (use warnings() to see them)
## >    >>>>
## >    >>>>> and as the last line does not show in the .save file, make fails.
## >    >>>>
## >    >>>>> Is anyone seeing this too?
## >    >>>>
## >    >>>> This persists for both GCC 4.3 and 4.4 for me, the warnings coming from
## >    >>>>
## >    R> xB <- c(2000,1e6,1e50,Inf)
## >    R> for(df in c(0.1, 1, 10))
## >    >>>> + for(ncp in c(0, 1, 10, 100)) stopifnot(pchisq(xB, df=df, ncp=ncp) == 1)
## >    >>>>
## >    >>>> -k
## >    >>>>
## >    >>>>> Best
## >    >>>>> -k
## >    >>>>

## Reproducing an "inverse" of Table 29.2 (p.464) of  Johnson, Kotz, Balakr.(1995) Vol.2

nu. <- c(2,4,7)
lam <- c(1,4,16,25)
(pnl <- expand.grid(ncp=lam, df=nu., KEEP.OUT.ATTRS=FALSE)[,2:1])
nl <- with(pnl, df+ncp)
pars <- rbind(cbind(pnl, q = nl/2),
              cbind(pnl, q = nl  ),
              cbind(pnl, q = nl*2))
pch   <- with(pars, pchisq(q=q, df=df, ncp=ncp))
pchAA <- with(pars, pnchisqAbdelAty  (q=q, df=df, ncp=ncp))
pchSa <- with(pars, pnchisqSankaran_d(q=q, df=df, ncp=ncp))
cbind(pars, R = pch, AA = pchAA, San = pchSa)
showProc.time()

### Reproducing part of  'Table 29.2' (p.464) of  Johnson, Kotz, Balakr.(1995) Vol.2
###
### as in ../man/pnchisqAppr.Rd -- do run over *all* current pnchisq*() approximations!

pkg <- "package:DPQ"
## NB: use versions of the functions that return numeric *vector* (of correct length) :
pnchNms <- c(paste0("pchisq", c("", "V", "W", "W.R")), # + R's own, but *not* "W." !
             ls(pkg, pattern = "^pnchisq"))
## drop some :
pnchNms <- pnchNms[!grepl("Terms$", pnchNms)]
pnchNms <- pnchNms[is.na(match(pnchNms, c("pnchisqIT", paste0("pnchisqT93.", c("a", "b")))))]
pnchF <- sapply(pnchNms, get, envir = as.environment(pkg))
## shorten the longer names for nicer tables :
n.n <- nchar(pnNms <- setNames(,pnchNms))
L8 <- n.n > 8
pnNms[n.n > 10] <- sub("pnchisq", "pn",  pnNms[n.n > 10])
pnNms[n.n >  8] <- sub("pnchisq","pnch", pnNms[n.n >  8])
names(pnchF) <- pnNms <- unname(abbreviate(pnNms, 8))
str(pnchF)
op <- options(warn = 1, digits = 5, width = 110)# warn: immediate ..
## TODO --- want also "x ~ ncp" and or "df ~ ncp"
## TODO: write a *function* that computes all this *and* stores in nicely dimnamed array
qq <- c(.001, .005, .01, .05, (1:9)/10, 2^seq(0, 10, by= 0.5))
nncp <- c(0, 1/8, 1/2, 1, 2, 5, 20, 100, 200, 1000)
ddf <- c(2:4, 7, 20, 50, 100, 1000, 1e4, 1e10) # 1e300: fails for pchisqW.R() << FIXME
AR <- array(NA_real_, # [ncp,df, q]
            dim=c(length(nncp), length(ddf), length(qq), length(pnchF)),
            dimnames= list(ncp = formatC(nncp, width=1),
                           df  = formatC( ddf, width=1),
                           q   = formatC(  qq, width=1),
                           Fn  = pnNms))
CT <- AR[,,1,1] # (w/ desired dim and dimnames)

sfil5 <- file.path(sdir, "tests_chisq-nonc-ssAp.rds")
if(!doExtras && file.exists(sfil5)) {
  ssAp_l <- readRDS_(sfil5)
  str(ssAp_l)
  AR <- ssAp_l$AR ## loadList(ssAp_l)# attach it

} else { ## do run the simulation always  if(doExtras) :

  for(incp in seq_along(nncp)) {
    cat("\n~~~~~~~~~~~~~\nncp: ", ncp <- nncp[incp], "\n=======\n")
    pnF <- if(ncp == 0) pnchF[!grepl("T93", pnNms)] else pnchF # Temme('93) : ncp > 0
    for(idf in seq_along(ddf)) {
        df <- ddf[idf]
        ct <- system.time(
          r <- vapply(pnF,
                      function(F) Vectorize(F, names(formals(F))[[1]])(qq, df=df, ncp=ncp),
                      qq)
        )[["user.self"]]
        AR[incp, idf, , names(pnF)] <- r
        CT[incp, idf] <- ct
    }
  }
  showProc.time()
  cat("User times in milli-sec.:\n")
  print(CT * 1000)
  save2RDS(list_(pnchNms, pnchF, qq, nncp, ddf,   AR, CT), file=sfil5)
} ## else *do* run ..

## Rather, show absolute and also relative "errors" ..
stopifnot(dimnames(AR)[[4]][1] == "pchisq")
## Absolute "error" , i.e., delta to R's pchisq() which is AR[,,,1] :
dAR <- AR[,,,-1] - c(AR[,,,1])
## if we were perfect, using same compilers as R etc, then these deltas are all = 0 :
summary(dAR[,,,"pnchRC"])
if(beStrict) stopifnot(dAR[,,,"pnchRC"] == 0)

apply(dAR, 4, summary) # quite some NA's for some Fn's
aperm(apply(dAR, 3:4, function(x) {u <- x[!is.na(x)]; c(min=min(u), max=max(u))}), c(2,3,1L))
ftable(apply(dAR[c("1", "2", "5"), c(1,3,4),,], c(1,2,4), median))
ftable(apply(dAR[c("1", "2", "5"), c(1,3,4),,], c(1,2,4), function(x) max(abs(x[!is.na(x)]))))

## (Note: for large df=1e10,  p*() == 0 everywhere for the small 'q' we have
##  ----  FIXME:  qq: should be "realistic"  in  mu +/- 5 sd

options(warn = 0, digits = 7)# partial revert

###----------- Much testing  pnchisqRC()  notably during my experiments
(ptol <- if(noLdbl) 8e-13 else if(doExtras) 3e-16 else if(is32) 1e-14 else 1e-15)
set.seed(123)
for(df in c(.1, .2, 1, 2, 5, 10, 20, 50, 1000,
            if(doExtras) c(1e10, 1e200))) { ## BUG!  (df=1e200, ncp=1000) takes forever
    cat("\n============\ndf = ",df,"\n~~~~~~~~~\n")
    for(ncp in c(0, .1, .2, 1, 2, 5, 10, 20, 50,
                 if(df < 1e10) c(1000, 1e4) else c(100,200))) { # BUG: large ncp take forever
        cat("\nncp = ",ncp,":  qq = ")
        qch <- if(ncp+df < 1000)
                   qchisq((1:15)/16, df=df, ncp=ncp)
               else {
                   qq <- qnchisqPatnaik((1:15)/16, df=df, ncp=ncp)
                   if(qq[1] < qq[length(qq)])
                       qq
                   else { # they all coincide ==> take mu +/- (1:3) SD
                       mu <- df+ncp
                       sigma <- sqrt(2*(df + 2*ncp))
                       mu + seq(-4,4, length.out=15)*sigma
                   }
               }
        str(qq <- c(0, qch, Inf), digits=4)
        for(lower.tail in c(TRUE, FALSE)) {
            cat(sprintf("lower.tail = %-5s :  ", lower.tail))
            for(log.p in c(FALSE, TRUE)) {
                cat("log.p=", log.p, "")
                AE <- all.equal(
                    pchisq   (qq, df=df, ncp=ncp, lower.tail=lower.tail, log.p=log.p) ,
                    pnchisqRC(qq, df=df, ncp=ncp, lower.tail=lower.tail, log.p=log.p) ,
                    tol = ptol)
                if(is.character(AE)) {
                    dd <- sub(".*:", "", AE)
                    cat("pchisq() differ by", dd,"(dd/ptol = ",as.numeric(dd)/ptol," < 100 ?)\n")
                    ## fails for first df=0.1, ncp=10000 on Windows 64-bit (winbuilder 2019-10)
                    ##  ... also on 'florence' (32-bit Fedora 28, 2019-10)
                    if(myPlatf || ncp <= 1000)
                        stopifnot(as.numeric(dd) < 100 * ptol)
                    else if (   !(as.numeric(dd) < 100 * ptol))
                        cat("not stop()ing even though dd < 100 * ptol\n")
                }
            }; cat("\n")
        }
    }# for(ncp .)
    showProc.time()
}# for(df .)
summary(warnings())

### L. Emphasis on very large  df + ncp ===============================================

##===  L 1.  very large df,  ncp/df << 1 ====================

mkPnch <- function(k, df, ncp, lower.tail=TRUE, log.p=FALSE, twoExp = -53) {
    stopifnot(is.numeric(k), length(k) > 1, k == (k <- as.integer(k)),
              is.numeric(df), length(df) == 1L, length(ncp) == 1L, ncp >= 0,
              if((k. <- min(k)) >= 0) TRUE else twoExp < -log2(-k.))
    ones <- 1 + k * 2^twoExp
    qs <- ones*(df+ncp) # df+ncp = E[chi'^2]
    xtra <-
        if(df == 1) { ## use exact formula {incl Taylor for small x = q}
            cbind(pnchi1sq = pnchi1sq(qs,  ncp, lower.tail=lower.tail, log.p=log.p))
        } else if(df == 3) {
            cbind(pnchi3sq = pnchi3sq(qs,  ncp, lower.tail=lower.tail, log.p=log.p))
        } else {
            array(NA_real_, c(length(qs), 0L))
        }
    cbind(pchisq   = pchisq          (qs,df,ncp, lower.tail=lower.tail, log.p=log.p),
          xtra,
          pcAbdelA = pnchisqAbdelAty (qs,df,ncp, lower.tail=lower.tail, log.p=log.p),
          pcBolKuz = pnchisqBolKuz   (qs,df,ncp, lower.tail=lower.tail, log.p=log.p),
          pcPatnaik= pnchisqPatnaik  (qs,df,ncp, lower.tail=lower.tail, log.p=log.p),
          pcPearson= pnchisqPearson  (qs,df,ncp, lower.tail=lower.tail, log.p=log.p),
          pcSanka_d=pnchisqSankaran_d(qs,df,ncp, lower.tail=lower.tail, log.p=log.p))
}

if(doExtras) { ## really slow because pchisq() is slow!

df <- 1e30; ncp <- 99
ks <- c(-40, -20, -15, -10, -6:6, 10, 15, 20, 40)
twoExp <- -25
##        ===
system.time(suppressWarnings(
    Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
)) ## 6.7 sec

print(Pn., digits=3) # ">>" annotated: shows bug !
##      pchisq pcAbdelA pcBolKuz pcPatnaik pcPearson pcSanka_d
##    0.00e+00      0.0      0.0       0.0       0.0       0.0
##    0.00e+00      0.0      0.0       0.0       0.0       0.0
##    0.00e+00      0.0      0.0       0.0       0.0       0.0
##    0.00e+00      0.0      0.0       0.0       0.0       0.0
##    0.00e+00      0.0      0.0       0.0       0.0       0.0
## >> 1.00e+00      0.0      0.0       0.0       0.0       0.0
## >> 1.00e+00      0.0      0.0       0.0       0.0       0.0
## >> 1.00e+00      0.0      0.0       0.0       0.0       0.0
## >> 1.00e+00      0.0      0.0       0.0       0.0       0.0
## >> 1.00e+00      0.0      0.0       0.0       0.0       0.0
##    4.17e-09      0.5      0.5       0.5       0.5       0.5
##    1.00e+00      1.0      1.0       1.0       1.0       1.0
##    1.00e+00      1.0      1.0       1.0       1.0       1.0
##    1.00e+00      1.0      1.0       1.0       1.0       1.0
##    1.00e+00      1.0      1.0       1.0       1.0       1.0
##    1.00e+00      1.0      1.0       1.0       1.0       1.0
##    .....         ....
tExp <- substitute(`pchisq*`(list(q == mu(1 + k * 2^TWOEXP), df==DF, ncp==NCP)),
                   list(TWOEXP = twoExp,  DF=df, NCP=ncp))
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)", main = tExp)

kk <- seq(min(ks), max(ks), length.out=401)
qs <- (1 + kk * 2^twoExp)*df
fq <- dchisq(qs, df, ncp)
par(new=TRUE)
plot(kk, fq, type="l", col=adjustcolor(2, 1/3), lwd=3, lty=3, axes=FALSE, ann=FALSE)

showProc.time()

}## only if(doExtras)
## BUG (FIXME) e.g. here:
 pchisq  (0.99999989*(df+ncp), df, ncp) ## --> Warning ... : not converged in 1000'000 iter
pnchisqRC(0.99999989*(df+ncp), df, ncp,
          verbose = 1) # The same with more output! ERROR on Winbuilder 64bit
## both give '1', but really should give 0
showProc.time()


### Much less extreme df ==> pchisq() *is* fast too
df <- 1e9
ncp <- 99
ks <- c(-40, -20, -15, -10, -6:6, 10, 15, 20, 40)
ks <- if(doExtras) -200:200 else seq(-200, 200, by=5)# more here for plot
twoExp <- -18 # well chosen for this "range" and behavior of pchisq()
##        ===
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
showProc.time()

tit <- substitute(`pchisq*`(list(q == mu(1 + k * 2^TWOEXP), df==DF, ncp==NCP)),
                  list(TWOEXP = twoExp,  DF=df, NCP=ncp))
matplot(ks, Pn., type = "l", xlab = quote(k), ylab = "pchisq*(q, ..)", main = tit)

cat("'Error' (difference to pchisq(*)):\n")
dP <- Pn.[,-1] - Pn.[,1]
print(cbind(ks, q=(1+ks*2^twoExp)*df, pchisq=Pn.[,1], dP), digits = 4)
matplot(ks, dP, type = "l", xlab = quote(k), main = paste("Difference", tit," - pchisq(..)"))
abline(h=0, lty=3)
## the difference to all 5 approx. is almost *IDENTICAL*
## ==> are the approximations all more accurate than pchisq() here ?
options(op) # revert
summary(warnings())

## Look at "smoothness" via first differences:
matplot(ks[-1], diff(Pn.), type = "l", xlab = quote(k))
abline(h=0, lty=3)

nk <- length(ks)
matplot(ks[-c(1,nk)], diff(Pn., differences= 2), type = "l", xlab = quote(k))
abline(h=0, lty=3)

matplot(ks[-c(1:2,nk)], diff(Pn., differences= 3), type = "l", xlab = quote(k))
abline(h=0, lty=3) ##---> start seeing noise

## Here, we see a LOT of noise : only in the first curve == pchisq() !
matplot(ks[-c(1:2,nk-1L,nk)], diff(Pn., differences= 4), type = "l", xlab = quote(k))
abline(h=0, lty=3)
## And zooming in to "zero" : log |.| scale
matplot(ks[-c(1:2,nk-1L,nk)], abs(diff(Pn., differences= 4)), log = "y", type = "l", xlab = quote(k))

showProc.time()

##===  L 2.  large ncp,  df/ncp << 1 ====================

pchiTit <- function(twoE, df, ncp, fN = "pchisq*", xtr = "", ncN = "ncp") {
    substitute(FUN(list(q == mu(1 + k* 2^TWO_E)*XTR, mu == {nu+lambda == DF+NCP}, df == DF, NCN == NCP)),
               list(FUN = fN, TWO_E = twoE, XTR=xtr, DF=df, NCN=ncN, NCP=ncp))
    ## sprintf("%s(q = μ(1 + k* 2^%g)%s, µ = ν+λ = df+ncp; df=%g, %s=%g)",
    ##         fN, twoE, xtr, df, ncN, ncp)
}

pchiTit.1 <- function(twoE, df, ncp)
    pchiTit(twoE, df, ncp, fN = "pchi*", xtr = " - pchi.1")
pchiTit.n.d <- function(twoE, df, ncp) pchiTit(twoE, df, ncp=ncp/df, ncN="ncp/df")

ks <- c(-40, -20, -15, -10, -6:6, 10, 15, 20, 40)

ncp <- 1e20; df <- 99
twoExp <- -35
##        ===
system.time(suppressWarnings(
    Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
)) ## ~ 0.3

print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
##  pchisq pcAbdelA pcBolKuz pcPatnaik pcPearson pcSanka_d
##       0 2.93e-09        1  2.93e-09  2.93e-09  2.93e-09
##       0 1.80e-03        1  1.80e-03  1.80e-03  1.80e-03
##       0 1.45e-02        1  1.45e-02  1.45e-02  1.45e-02
##       0 7.28e-02        1  7.28e-02  7.28e-02  7.28e-02
##       0 1.91e-01        1  1.91e-01  1.91e-01  1.91e-01
##       0 2.33e-01        1  2.33e-01  2.33e-01  2.33e-01
##       0 2.80e-01        1  2.80e-01  2.80e-01  2.80e-01
##       0 3.31e-01        1  3.31e-01  3.31e-01  3.31e-01
##       0 3.86e-01        1  3.86e-01  3.86e-01  3.86e-01
##       0 4.42e-01        1  4.42e-01  4.42e-01  4.42e-01
##       0 5.00e-01        1  5.00e-01  5.00e-01  5.00e-01
##       0 5.58e-01        1  5.58e-01  5.58e-01  5.58e-01
##       0 6.14e-01        1  6.14e-01  6.14e-01  6.14e-01
##       0 6.69e-01        1  6.69e-01  6.69e-01  6.69e-01
##       0 7.20e-01        1  7.20e-01  7.20e-01  7.20e-01
##       0 7.67e-01        1  7.67e-01  7.67e-01  7.67e-01
##       0 8.09e-01        1  8.09e-01  8.09e-01  8.09e-01
##       0 9.27e-01        1  9.27e-01  9.27e-01  9.27e-01
##       0 9.85e-01        1  9.85e-01  9.85e-01  9.85e-01
##       0 9.98e-01        1  9.98e-01  9.98e-01  9.98e-01
##       1 1.00e+00        1  1.00e+00  1.00e+00  1.00e+00

matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)", main= pchiTit(twoExp,df,ncp))

if(doExtras) {
## less extreme, same phenomenom:
ncp <- 1e9; df <- 99 ; twoExp <- -17
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)", main= pchiTit(twoExp,df,ncp))

## less extreme, same phenomenom: still
ncp <- 1e7; df <- 99 ; twoExp <- -14
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)", main= pchiTit(twoExp,df,ncp))

## even less extreme, same phenomenom: still
ncp <- 4e6; df <- 99 ; twoExp <- -13
##     ---
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # pchisq = Pearson = Sanka_d  ~~ AbdelA, Patnaik
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)", main= pchiTit(twoExp,df,ncp))

} # only if(doEx..)



showProc.time()

## Here pchisq() seems "perfect"
ncp <- 1e6; df <- 99 ; twoExp <- -12
##     ---
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # pchisq = Pearson = Sanka_d  ~~ AbdelA, Patnaik
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)", main= pchiTit(twoExp,df,ncp))

kk <- seq(min(ks), max(ks), length.out=401)
qs <- (1 + kk * 2^twoExp)*(df+ncp)
fq <- dchisq(qs, df, ncp)
par(new=TRUE)
plot(kk, fq, type="l", col=adjustcolor(2, 1/3), lwd=3, lty=3, axes=FALSE, ann=FALSE)
showProc.time()

##=== df=1 and df=3 ==== here we have exact formula ! ==================

## 10'000 : "too small" for asymptotic approx:
ncp <- 10000; df <- 1 ; twoExp <- -7
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.1(twoExp,df,ncp))

## now absolute *difference* to true pchi1sq() :
print(Pn.[,-c(2,4)]-Pn.[,2], digits=3)

matplot(ks, Pn.[,-c(2,4)]-Pn.[,2], type = "b", xlab = quote(k), ylab = "pchisq*(q, ..) - pchi1sq()",
        main = pchiTit.1(twoExp,df,ncp))
legend("topright", colnames(Pn.)[-c(2,4)], lty=1:5, col=1:5, bty="n")

j.dr <- 2:5 # drop
matplot(ks, Pn.[,-j.dr]-Pn.[,2], type = "b", xlab = quote(k), ylab = "pchisq*(q, ..) - pchi1sq()",
        main = pchiTit.1(twoExp,df,ncp))
legend("topright", colnames(Pn.)[-j.dr], lty=1:5, col=1:5, bty="n")

j.d2 <- 2:6 # drop
matplot(ks, Pn.[,-j.d2]-Pn.[,2], type = "b", xlab = quote(k), ylab = "pchisq*(q, ..) - pchi1sq()",
        main = pchiTit.1(twoExp,df,ncp))
legend("topright", colnames(Pn.)[-j.d2], lty=1:5, col=1:5, bty="n")

showProc.time()# --

## 1e6 : ???
ncp <- 1e6; df <- 1 ; twoExp <- -13
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
showProc.time()
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.1(twoExp,df,ncp))
## now absolute *difference* to true pchi1sq() :
print(Pn.[,-c(2,4)]-Pn.[,2], digits=3)

matplot(ks, Pn.[,-c(2,4)]-Pn.[,2], type = "b", xlab = quote(k), ylab = "pchisq*(q, ..) - pchi1sq()",
        main = pchiTit.1(twoExp,df,ncp))
legend("topright", colnames(Pn.)[-c(2,4)], lty=1:5, col=1:5, bty="n")

j.dr <- 2:5 # drop
matplot(ks, Pn.[,-j.dr]-Pn.[,2], type = "b", xlab = quote(k), ylab = "pchisq*(q, ..) - pchi1sq()",
        main = pchiTit.1(twoExp,df,ncp))
legend("topright", colnames(Pn.)[-j.dr], lty=1:5, col=1:5, bty="n")

## Convincing: here, Sanka_d  is *better* than R's pchisq() !
j.d2 <- 2:6 # drop
matplot(ks, Pn.[,-j.d2]-Pn.[,2], type = "b", xlab = quote(k), ylab = "pchisq*(q, ..) - pchi1sq()",
        main = pchiTit.1(twoExp,df,ncp))
abline(h=0, lty=3, col=adjustcolor(1, 1/4))
legend("topright", colnames(Pn.)[-j.d2], lty=1:5, col=1:5, bty="n")
showProc.time()

##     vvvv
ncp <- 1e9; df <- 3 ; twoExp <- -17
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit(twoExp,df,ncp))
showProc.time()


##===  L 3.  BOTH large ncp, large df ====================


df <- 1e9; ncp <- 1 * df ; twoExp <- -17
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.n.d(twoExp,df,ncp))

if(doExtras) { ## because it's a bit costly here :
df <- 1e8; ncp <- 2 * df ; twoExp <- -16
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.n.d(twoExp,df,ncp))

df <- 1e7; ncp <- 1/2 * df ; twoExp <- -14
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.n.d(twoExp,df,ncp))

## pnchisq still "broken"; others good
df <- 1e6; ncp <- 10 * df ; twoExp <- -14
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.n.d(twoExp,df,ncp))

## pchisq *NON*-monotone !!
df <- 4e6; ncp <- .5 * df ; twoExp <- -14
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.n.d(twoExp,df,ncp))

} # only if(.X.)
showProc.time()

## pchisq ok
df <- 2e6; ncp <- .5 * df ; twoExp <- -14
Pn. <- mkPnch(ks, df=df, ncp=ncp, twoExp=twoExp)
print(Pn., digits=3) # ">>" pcPolKuz is *NOT* for this; R's pchisq is full wrong; other "coincide"
matplot(ks, Pn., type = "b", xlab = quote(k), ylab = "pchisq*(q, ..)",
        main = pchiTit.n.d(twoExp,df,ncp))


showProc.time()

### Part 3 :  qchisq (non-central!)
### -------------------------------

if(!dev.interactive(orNone=TRUE)) { dev.off(); pdf("chisq-nonc-3.pdf") }

### Bug 875 {see also ~/R/r-devel/R/tests/d-p-q-r-tests.R
(q49.7 <- qchisq(0.025, 31, ncp=1, lower.tail=FALSE))## now ok: 49.7766
pb     <- pchisq(q49.7, 31, ncp=1, lower.tail=FALSE)
all.equal(pb, 0.025, tol=0) # 2.058e-13 [Lnx 64b]; 2.0609e-13 [Win 32b]
stopifnot(all.equal(pb, 0.025, tol= 1e-12))

##  Ensuing things I tried :
x <- seq(0, 20, len = 101)
plot(x, pnc <- pchisq(x, 5, ncp = 1.1), type = 'l', col = 'red')
xx <-        qchisq(pnc, 5, ncp = 1.1)
stopifnot(all.equal(x, xx))#TRUE
all.equal(x, xx, tol = 0) # 1.9012e-13, later 1.835842e-14 (Linux)

plot(x, pncR <- pchisq(x, 5, ncp = 1.1, lower = FALSE), type = 'l', col = 'red')
(pnc + pncR) - 1
stopifnot(all.equal(pnc + pncR, rep(1, length(pnc))))
xx0 <- qchisq(pncR, 5, ncp = 1.1, lower = FALSE)
all.equal( x, xx0, tol = 0) # 1.877586e-13; then 1.8364e-14
all.equal(xx, xx0, tol = 0) # 5.942721e-13; then 6.2907e-13, 6.2172e-14
stopifnot(all.equal(x, xx0))

plot(x, LpncR <- pchisq(x, 5, ncp = 1.1, lower = FALSE, log = TRUE),
     type = 'l', col = 'red')
Lxx0 <- qchisq(LpncR, 5, ncp = 1.1, lower = FALSE, log = TRUE)
all.equal(x, Lxx0, tol = 0)# 1.8775..e-13; 1.8364e-14
all.equal(log(pncR),    LpncR, tol = 0)# 0,             now 2.246e-16
all.equal(log(1 - pnc), LpncR, tol = 0)# 4.661586e-16;  now 2.246e-16
all.equal(log1p(- pnc), LpncR, tol = 0)# 4.626185e-16; now TRUE

showProc.time()
## source("/u/maechler/R/MM/NUMERICS/dpq-functions/qnchisq.R")#-> qnchisq.appr*()

## The values from Johnson et al (1995), Table 29.2, p.464
p.  <- c(0.95, 0.05)
nu. <- c(2,4,7)
lam <- c(1,4,16,25)
str(pars <- expand.grid(ncp=lam, df=nu., p= p., KEEP.OUT.ATTRS=FALSE)[,3:1])
## 'data.frame':  24 obs. of  3 variables:
##  $ p  : num  0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 ...
##  $ df : num  2 2 2 2 4 4 4 4 7 7 ...
##  $ ncp: num  1 4 16 25 1 4 16 25 1 4 ...

qch <- with(pars, qchisq(p=p, df=df, ncp=ncp))
p.q <- with(pars, pchisq(qch, df=df, ncp=ncp))

cbind(pars, qch, p.q, relE = signif(1 - p.q / pars$p, 4)) ## very accurate
##       p df ncp        qch  p.q       relE
## 1  0.95  2   1  8.6422039 0.95  3.331e-16
## 2  0.95  2   4 14.6402116 0.95  3.442e-15
## 3  0.95  2  16 33.0542126 0.95 -8.882e-15
## 4  0.95  2  25 45.3082281 0.95  2.887e-15
## 5  0.95  4   1 11.7072278 0.95  5.662e-15
## 6  0.95  4   4 17.3093229 0.95  5.995e-15
## 7  0.95  4  16 35.4270110 0.95 -1.199e-14
## 8  0.95  4  25 47.6127674 0.95  8.771e-15
## 9  0.95  7   1 16.0039003 0.95 -5.551e-15
## 10 0.95  7   4 21.2280338 0.95  6.661e-15
## 11 0.95  7  16 38.9700904 0.95  1.110e-14
## 12 0.95  7  25 51.0605938 0.95 -1.488e-14
## 13 0.05  2   1  0.1683911 0.05 -2.398e-14
## 14 0.05  2   4  0.6455990 0.05  3.020e-14
## 15 0.05  2  16  6.3216416 0.05  5.673e-14
## 16 0.05  2  25 12.0802051 0.05 -1.477e-13
## 17 0.05  4   1  0.9087447 0.05  5.385e-14
## 18 0.05  4   4  1.7650116 0.05  6.106e-14
## 19 0.05  4  16  7.8843284 0.05 -7.105e-15
## 20 0.05  4  25 13.7329249 0.05  8.271e-14
## 21 0.05  7   1  2.4937057 0.05 -6.128e-14
## 22 0.05  7   4  3.6642526 0.05 -2.420e-14
## 23 0.05  7  16 10.2573190 0.05  1.066e-14
## 24 0.05  7  25 16.2267524 0.05  3.775e-14
all.equal(pars$p, p.q, tol=0)# Lnx 64b: 9.2987e-15; Win 32b: 9.25e-15
stopifnot(all.equal(pars$p, p.q, tol=1e-14))
showProc.time()

## now works fine :
str(n.s3 <- newton(1, G= function(x,...) x^2 -3 , g = function(x,...) 2*x,
                   eps = 8e-16))
with(n.s3, stopifnot(converged, all.equal(x^2, 3, tol = 1e-15)))


### New comparison -- particularly for right tail:
## upper tail "1 - p"

p.qappr <- function(p, df, ncp, main = NULL,
                    kind = c("raw", "diff", "abs.Err", "rel.Err"),
                    nF = NULL, do.title= is.null(main), do.legend = TRUE,
                    ylim.range = 0.4, ...)
{
    ## Purpose: Plot comparison of different  qchisq() approximations
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Dr. Martin Maechler, Date: 27 Feb 2004, 18:19
    kind <- match.arg(kind)
    d.arg <- (l.d <- length(df)) > 1
    n.arg <- (l.n <- length(ncp)) > 1
    p.arg <- (l.p <- length(p)) > 1
    if((p.arg && (d.arg || n.arg)) || (d.arg && n.arg))
        stop("only one of the three argument should have length > 1")
    if(!(d.arg || n.arg || p.arg)) p.arg <- TRUE
    n <- max(l.d, l.n, l.p)
    Fns <- c("qchisq", "qnchisqPearson",
             "qchisqApprCF1", "qchisqApprCF2",
             "qnchisqPatnaik", "qchisqCappr.2",
             "qchisqAppr.0", "qchisqAppr.1", "qchisqAppr.2", "qchisqAppr.3"
             )
    if(is.null(nF))
        nF <- length(Fns)
    else if(is.numeric(nF) & nF > 1) Fns <- Fns[1:nF]
    else stop("invalid 'nF' argument")
    qmat <- matrix(NA, n, length(Fns), dimnames=list(NULL,Fns))
    for(i.f in 1:nF) {
        fn <- Fns[i.f]
        F <- get(fn)
        qmat[,i.f] <- do.call(F, list(p=p, df=df, ncp=ncp, lower.tail=FALSE))
    }
    cols <- 1:nF
    lwds <- c(2, rep(1,nF-1))
    ltys <- rep(1:3, length = nF)
    if(kind != "raw") {
        cols <- cols[-1]
        lwds <- lwds[-1]
        ltys <- ltys[-1]
        Fns <- Fns[-1]
        ## Inf - Inf = 0  in the following :
        "%-%" <- function(x,y)
            ifelse(is.infinite(x) & is.infinite(y) & x==y, 0, x-y)
        qm <- qmat[,-1, drop=FALSE] %-% qmat[,1]
        if(kind != "diff") {
            qm <- abs(qm)
            if(kind == "rel.Err") qm <- qm / abs(qmat[,1])
        }
        yl <- rrange(qm, r = ylim.range)
    } else {
        qm <- qmat
        yl <- range(qmat[,"qchisq"], rrange(qmat, r = ylim.range))
    }
    if(do.title && is.null(main)) main <- deparse(match.call())
    matplot(if(p.arg) p else if(d.arg) df else ncp, qm, type = 'l',
            xlab = if(p.arg) "1 - p" else if(d.arg) "df" else "ncp",
            ylim = yl, main = main, col = cols, lwd= lwds, lty= ltys, ...)
    if(do.title)
        mtext("different approximations to  qchisq()", line = 0.4, cex = 1.25)
    ## drop "qn?chisq" from names, since have it above:
    if(do.legend) {
        Fns <- sub("^qn?chisq.","*",Fns)
        pu <- par("usr")
        legend(par("xaxp")[2], par("yaxp")[2],
               Fns, xjust = 1.02, yjust = 1.02, ncol = 3,
               col = cols, lwd=lwds, lty= ltys)
    }
    invisible(qmat)
} ## end{ p.qappr() }

pU <- seq(.5, 1, length= 201)
pU <- seq( 0, 1, length= 501)[-c(1,501)]
## (I've lost the original 'pU' I had used ...)

mystats <- function(x) c(M=mean(x), quantile(x))
sum.qappr <- function(r) {
    m <- t(apply(abs(r[,-1] - r[,1]), 2,mystats))
    m[order(m[,"50%"]),]
}
op <- options(digits = 6, width = 110)# warn: immediate ..
showProc.time()

sum.qappr(p.qappr (pU, df= 1, ncp= 1))

sum.qappr(r <- p.qappr (pU, df=10, ncp= 10))
## just different pictures:
               p.qappr (pU, df=10, ncp= 10, kind = "diff", ylim.r = 1)
               p.qappr (pU, df=10, ncp= 10, kind = "abs", ylim.r = 0.01)
               p.qappr (pU, df=10, ncp= 10, kind = "rel", log = 'y')
showProc.time()

sum.qappr(p.qappr (pU, df= 1, ncp= 10))
sum.qappr(p.qappr (pU, df= 1, ncp= 10, kind="rel"))
showProc.time()

if(doExtras) ## this takes CPU !
sum.qappr(p.qappr (pU, df= 10, ncp= 1e4, kind="rel"))
showProc.time() # 2.9 sec

##--> CF2, Pea, CF1 and Patn  are  the four best ones overall
##    ---  ---  ---     ----

### Now look at upper tail only: ----- even a more clear picture
summary(pU <- 2^-seq(7,40, length=200))
sum.qappr(r <- p.qappr (pU, df= 1, ncp= 10))
               p.qappr (pU, df= 1, ncp= 10, kind="rel")
sum.qappr(r <- p.qappr (pU, df= 1, ncp= 100))
## very small ncp:
sum.qappr(r <- p.qappr (pU, df= 1, ncp= .01))
               p.qappr (pU, df= 1, ncp= .01, kind="dif", log="x", nF =6)
               p.qappr (pU, df= 1, ncp= .01, kind="rel", log="xy",nF =6)
# here, CF2 is "off" and the top is  "Cappr.2", "Pea", "Patn" ("CF1", "appr.3")

sum.qappr(r <- p.qappr (pU, df= 10, ncp= .01))
                                        # shows noise in qchisq() itself !?
               p.qappr (pU, df= 10, ncp= .01, kind="rel", log="xy")
# "CF2" +-ok;  top is   "Cappr.2", "Pea", "Patn" ("CF2", "appr.3")
sum.qappr(r <- p.qappr (pU, df= 100, ncp= .01))
                                        # shows noise in qchisq() itself !!!
               p.qappr (pU, df= 100, ncp= .01, kind="rel", log="xy")
showProc.time()

## even smaller ncp:
sum.qappr(r <- p.qappr (pU, df= 100, ncp= .001))
                                        # shows noise in qchisq() itself !!!
               p.qappr (pU, df= 100, ncp= .001, kind="rel", log="xy")

sum.qappr(r <- p.qappr (pU, df= 1, ncp= .1))
               p.qappr (pU, df= 1, ncp= .1, kind="rel", log="y")
summary(warnings())
showProc.time()

sum.qappr(r <- p.qappr (pU, df= 20, ncp= 200))
               p.qappr (pU, df= 20, ncp= 200, kind='rel', log='y')

sum.qappr(r <- p.qappr (pU, df= .1, ncp= 500))
## drop the appr.<n> : they are bad
               p.qappr (pU, df= .1, ncp= 500, kind='dif', nF = 6)
               p.qappr (pU, df= .1, ncp= 500, kind='rel', nF = 6, log='y')
sum.qappr(r <- p.qappr (pU, df= .1, ncp= 500, kind='dif', nF = 6))
               p.qappr (pU, df= .1, ncp= 500, kind='rel', log='y', nF = 6)
# order: CF2, CF1, Pea, Patn
showProc.time()


### Very large ncp, df --- Sankaran_d and Pearson had failed !

op <- options(warn = 1)# immediate ..
pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
for(DF in 10^c(50, 100,150,200, 250, 300))
  stopifnot(exprs = {
    qnchisqPearson   (pp, df=DF, ncp=100) == DF
    qnchisqSankaran_d(pp, df=DF, ncp=100) == DF
    qnchisqPatnaik   (pp, df=DF, ncp=100) == DF
})

qtol <- if(doExtras) 3e-16 else 1e-15
## Both large df & large ncp
for(NCP in 10^c(50, 100,150,200, 250, 300))
  stopifnot(exprs = {
    abs(1 - qnchisqPearson   (pp, df=2*NCP, ncp=NCP) / (3*NCP)) < qtol
    abs(1 - qnchisqSankaran_d(pp, df=2*NCP, ncp=NCP) / (3*NCP)) < qtol
    abs(1 - qnchisqPatnaik   (pp, df=2*NCP, ncp=NCP) / (3*NCP)) < qtol

    abs(1 - qnchisqPearson   (pp, df=NCP, ncp=NCP) / (2*NCP)) < qtol
    abs(1 - qnchisqSankaran_d(pp, df=NCP, ncp=NCP) / (2*NCP)) < qtol
    abs(1 - qnchisqPatnaik   (pp, df=NCP, ncp=NCP) / (2*NCP)) < qtol

    abs(1 - qnchisqPearson   (pp, df=NCP/2, ncp=NCP) / (1.5*NCP)) < qtol
    abs(1 - qnchisqSankaran_d(pp, df=NCP/2, ncp=NCP) / (1.5*NCP)) < qtol
    abs(1 - qnchisqPatnaik   (pp, df=NCP/2, ncp=NCP) / (1.5*NCP)) < qtol
})
showProc.time()
options(op) # revert


DF <- 1e200
if(FALSE)## BUG (2019-08-31):
  system.time(
      qch <- qchisq(pp, df=DF, ncp=100)
  )## gives these warnings (immediately "warn = 1"),  and then takes "forever" !!
## Warning in qchisq(pp, df = DF, ncp = 100) :
##   pnchisq(x=1e+200, ..): not converged in 10000 iter.
## Warning in qchisq(pp, df = DF, ncp = 100) :
##   pnchisq(x=1e+200, ..): not converged in 10000 iter.

## "forever":  Timing stopped at: 3871 0.184 3878  > 1 hour

showProc.time()
