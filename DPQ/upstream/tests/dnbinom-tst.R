#### Testing  1) dbinom_raw(), dnbinomR() and dnbinom.mu()
####          2) log1pmx(), logcf() etc
require(DPQ)

source(system.file(package="DPQ", "test-tools.R",
                   mustWork=TRUE))# ../inst/test-tools.R
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()
options(warnPartialMatchArgs = FALSE)

(doExtras <- DPQ:::doExtras() && !grepl("valgrind", R.home()))

do.pdf <- TRUE
do.pdf <- !dev.interactive(orNone = TRUE)

### 1. Testing  dbinom_raw(), dnbinomR() and dnbinom.mu() >>> ../R/dbinom-nbinom.R <<<
### ----------   ../man/dbinom_raw.Rd & ../man/dnbinomR.Rd

eaxis     <- sfsmisc::eaxis
relErrV   <- sfsmisc::relErrV
## "FIXME:" use relErrV() already here

###  dbinom() vs  dbinom.raw() :

for(n in 1:20) {
    cat("n=",n," ")
    for(x in 0:n)
        cat(".")
        for(p in c(0, .1, .5, .8, 1)) {
            stopifnot(all.equal(dbinom_raw(x, n, p, q=1-p, log=FALSE),
                                dbinom    (x, n, p,        log=FALSE)),
                      all.equal(dbinom_raw(x, n, p, q=1-p, log =TRUE),
                                dbinom    (x, n, p,        log =TRUE)))
    }
    cat("\n")
}
showProc.time()

###  dnbinom*() :
stopifnot(exprs = {
    dnbinomR(0, 1, 1) == 1
})

## Both of these gave errors in DPQ <= 0.6-0  -- fixed w/ ifelse():
## Error in if (x < size) log1p(-x/(size + x)) ..... the condition has length > 1
x <- c(1:9, 10^c(1:2,50*(1:6)))
(dnbX <- dnbinomR  (x, size = 1.5e156, prob = .999, log = TRUE, eps = 0.001))
(dnbM <- dnbinom.mu(x, size = 1.5e156, mu = 1e150,  log = TRUE, eps = 0.001))


### exploring 'eps' == "true" tests must be done with  Rmpfr !!

### 2. Testing  log1pmx(), logcf() etc
### ----------

do.pdf
if(do.pdf) { pdf.options(width = 9, height = 6.5)# for all {9/6.5 = 1.38 ; 4/3 < 1.38 < sqrt(2) [A4]
    pdf("dnbinom-logcf.pdf")
}

### 2a:  logcf()
##  ==   =======
x <- c((-20:3)/4, (25:31)/32) # close (but not too close) to upper bound 1

(lC   <- logcf (x, i=2, d=3, eps=1e-9))
 lCt  <- logcf (x, i=2, d=3, eps=1e-9, trace=TRUE) ; stopifnot(identical(lCt, lC))
(lR   <- logcfR(x, i=2, d=3, eps=1e-9))
          all.equal(lC, lR, tol = 0) # x86_64 F40: 6.54e-16
stopifnot(all.equal(lC, lR, tol = 4e-15)) # 1.08295e-15 (Apple clang 14.0.3)
lRt  <- logcfR(x, i=2, d=3, eps=1e-9, trace=TRUE) ; stopifnot(identical(lRt, lR))
lRt2 <- logcfR(x, i=2, d=3, eps=1e-9, trace= 2)   ; stopifnot(identical(lRt2,lR))

lR.  <- logcfR(x, i=2, d=3, eps=1e-9)
lR.t <- logcfR(x, i=2, d=3, eps=1e-9, trace=TRUE) ; stopifnot(identical(lR.t, lR.))

all.equal(lC, lR., tol = 0) # no longer TRUE, but really small
all.equal(lR, lR., tol = 0) # TRUE !!   "    "
stopifnot(all.equal(lC, lR., tol = 1e-14))
## (even though they used eps=1e-9 .. i.e., are not *so* accurate)
lR.14 <- logcf(x, i=2, d=3, eps=1e-14)# default when used in R, incl. from log1pmx():
lR.18 <- logcf(x, i=2, d=3, eps=1e-18)

showProc.time()


##require(Rmpfr) may be not, see if NS loading (via "::") is sufficient:
requireNamespace("Rmpfr") || quit("no")
##                -----      ----------
asNumeric <- Rmpfr::asNumeric
all.equal <- Rmpfr::all.equal
mpfr      <- Rmpfr::mpfr
getPrec   <- Rmpfr::getPrec
.getPrec  <- Rmpfr::.getPrec

xM <- mpfr(x, 512)
(ct.24 <- system.time(lRM24 <- logcfR   (xM, i=2, d=3, eps=1e-24, trace=TRUE))) # it=76, 0.73 sec
if(doExtras)  withAutoprint({ # -----------------------------------
  ct24 <- system.time(lR24 <- logcfR_vec(xM, i=2, d=3, eps=1e-24, trace=TRUE)); ct24 #    5.7 sec
  all.equal(lRM24, lR24, tol=0) #  TRUE
  identical(lRM24, lR24)        #  TRUE !! (not sure if on all platforms!)
  stopifnot(all.equal(lRM24, lR24, tol = 5e-16))
})

## ===> use logcfR() {the internally vectorized version from now on)

SS <- function(ch, digits = 4) sub(paste0("([0-9]{1,",digits,"})[0-9]*e"), "\\1e", ch)
## double prec <--> MPFR:    vvvv (same eps)
lRM9 <- logcfR(xM, 2,3, eps=1e-9)
## show:
SS(all.equal(Rmpfr::roundMpfr(lRM9, 64), lR, tol=0))# .. 5.1138e-16
stopifnot(all.equal(lRM9,  lR  , tol=1e-15))
SS(all.equal(lRM24, lRM9, tol=0))# .. 3.701..e-10
stopifnot(all.equal(lRM24, lRM9, tol=1e-9))
showProc.time()
## now see if small eps makes a relevant difference:
relE14 <- asNumeric(relErrV(lRM24, lR.14))
relE18 <- asNumeric(relErrV(lRM24, lR.18))
summary(cbind(relE14, relE18))
##     relE14               relE18
## Min.   :-8.766e-15   Min.   :-9.169e-17
## 1st Qu.:-9.971e-16   1st Qu.:-3.435e-17
## Median : 1.268e-15   Median :-7.028e-20
## Mean   : 1.345e-15   Mean   :-3.644e-18
## 3rd Qu.: 3.843e-15   3rd Qu.: 3.514e-17
## Max.   : 1.090e-14   Max.   : 7.480e-17
matplot(x, pmax(abs(cbind(relE14, relE18)), 2^-55), log = "y", type = "l", ylim = c(3e-17, 1e-14),
        panel.first = abline(h= 2^-(54:51), lty=3, lwd = c(1,1,3,1), col="lightgray"))
title("|relative Errors| (wrt mpfr precBits = 512) of  logcf(x, 2,3, eps = *)")
legend("left", paste("eps = ", format(c(1e-14, 1e-18))), col=1:2, lty=1:2)

## zoom into [0, 1) -- which is the domain of  y = (x/(x+2)^2  for log1pmx(x) :
N <- if(doExtras) 1024 else 128
x <- (0:(N-1))/N # close (but not too close) to upper bound 1
xM <- mpfr(x, 512)
system.time(lgRM <- logcfR(xM, i=2, d=3, eps=1e-24, trace=doExtras))# it=152; 1.5 s {doEx.. it=425, 5.6 s}
relE14 <- asNumeric(relErrV(lgRM, logcf(x, i=2, d=3, eps=1e-14)))
relE18 <- asNumeric(relErrV(lgRM, logcf(x, i=2, d=3, eps=1e-18)))

summary(cbind(relE14, relE18))
##     relE14               relE18
## Min.   :-5.935e-14   Min.   :-1.087e-16
## 1st Qu.:-2.565e-15   1st Qu.:-3.831e-17
## Median :-6.919e-16   Median :-1.261e-18
## Mean   :-2.274e-15   Mean   : 9.056e-19
## 3rd Qu.:-1.144e-16   3rd Qu.: 4.275e-17
## Max.   : 9.903e-17   Max.   : 1.155e-16
matplot(x, pmax(abs(cbind(relE14, relE18)), 2^-55), log = "y", type = "l", ylim = c(3e-17, 4e-14),
        ylab = "", main = "logcf(x, 2,3, eps = *):  |relative Errors| wrt mpfr precBits = 512",
        panel.first = abline(h= 2^-(54:51), lty=3, lwd = c(1,1,3,1), col="lightgray"), yaxt="n")
eaxis(2)
legend("topleft", paste("eps = ", format(c(1e-14, 1e-18))), col=1:2, lty=1:2, bty = "n")
## Wow !!! ===> definitely should use a smaller eps  or  tol_logcf !!

if(do.pdf) { dev.off(); pdf("dnbinom-log1pmx.pdf") }

### 2b:  log1pmx(x) --  calls  logcf(y, 2,3), for y = t^2 = (x/(x+2))^2  in  [0,1)  for x > -1
##  ==   =========

## The first part is from original ../man/log1pmx.Rd

e <- if(doExtras) 2^-12 else 2^-8; by.p <- 1/(if(doExtras) 256 else 64)
xd <- sort(c(seq(-1+e, 0+100*e, by=e), seq(by.p, 5, by=by.p))) # length 676 or 5476 if do.X.
plot(xd, log1pmx(xd), type="l", col=2, lwd=2, main = "log1pmx(x)")
abline(h=0, v=-1:0, lty=3)

xM <- mpfr(xd, 512)
## for MPFR numbers, really better than tol_logcf = eps = 1e-14 (default):
if(doExtras) {
  print( system.time(
        lg1pM <- log1pmx(xM, tol_logcf = 1e-25, eps2 = 1e-4)
  )) # 2.3 sec [new logcfR()]
  lines(xd, asNumeric(lg1pM), col=adjustcolor(4, 1/4), lwd=3)
}
lg1pdR <- log1pmx(xd, trace=TRUE, logCF=logcfR)
lines(xd, lg1pdR, col=adjustcolor(5, 1/4), lwd=5)
##
lg1pM. <- log1p(xM) - xM # good enough if(precBits are high enough)
if(doExtras) {
  print(summary(relEM <- asNumeric(relErrV(lg1pM, lg1pM.))))
  ##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
  ## -7.910e-28  0.000e+00  5.600e-33  2.282e-28  2.732e-29  1.104e-26
  stopifnot(abs(relEM) < 1e-25)
  ## the error of the log1pmx() "algorithm" even for "perfect" mpfr-accuracy:
  xM2k <- mpfr(xd, 2048) # even more bits
  lg1pM2k <- log1p(xM2k) - xM2k
  reE00 <- asNumeric( relErrV(lg1pM2k, lg1pM.) )
  print(signif(rrre <- range(reE00)), 2) # [-1.5e-151, 4.8e-151] ==> 512 bits is sufficient
  stopifnot(abs(rrre) < 1e-150)
  plot(reE00 ~ xd, type="l") # see nothing but close to x ~= 0
  plot(abs(reE00) ~ xd, type="l", log="y")
  ## zoom in close to x ~= 0:
  plot(reE00 ~ xd, type="l", subset = abs(xd) <= 1/64)
  plot(abs(reE00) ~ xd, type="l", log="y", subset = abs(xd) <= 1)
  rE.log1pm <- asNumeric(relErrV(lg1pM2k, lg1pM))
  print(summary(  rE.log1pm ) )  # -1.1e-26 7.9e-28 {"matching" tol_logcf = 1e-25 above}
}
showProc.time()
re <- asNumeric(relErrV(lg1pM., log1pmx(xd)))

## MM: From around here, "move" to vignette  <<<<<<<<<<<<<<<<<<<<<<< ../vignettes/log1pmx-etc.Rnw <<<

plot(xd, re, type="b", cex=1/4, xlab = quote(x),
     main = paste("relErrV(log1p(xM) - xM, log1pmx(x)),  xM <- mpfr(x,",min(.getPrec(xM)),")"))
abline(h = (-2:2)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))
xl <- -0.84; xr <- -0.4 ##  the "relevant x-range":
colR <- adjustcolor("tomato", 1/2)
abline(v= c(xl,xr), col=colR, lty=2, lwd=2)
axis(1,at=c(xl,xr), labels=NULL, col.axis="tomato", col=colR, lwd=2, cex.axis=3/4)
text(xl, adj=c(0,-1/2), par("usr")[3], "zoom in", col=colR)
showProc.time()

# only "relevant"  x-range: values inside (xl, xr)
iN <- xl < xd & xd < xr
x.iN <- xd[iN]
plot(x.iN, re[iN], type="b", cex=1/4, main="rel.Error of log1pmx(x)",
     xlab=quote(x), ylab=quote(rE(x)), sub = "(in relevant x-range)")
abline(h = (-2:2)*2^-52, lty=2, col=adjustcolor("gray30", 1/2))
abline(v = c(xl,xr), col=colR, lty=2, lwd=2)

## *absolute* relative errors from here on:
plot(x.iN, abs(re[iN]), type="b", cex=1/2, log="y",
     main = "| relErr( log1pmx(x) ) |  {via 'Rmpfr'}",
     ylab=quote(abs(rE(x))), xlab=quote(x),
     ylim = c(4e-17, max(abs(asNumeric(re)[iN]))))
abline(h = c(1:2,4)*2^-52, lty=2, col=adjustcolor("gray20", 1/2))
mL1 <- eval(formals(log1pmx)$minL1)
abline(v = mL1, lwd=3, col=adjustcolor(2, 1/2))
axis(3, at=mL1, sprintf("minL1=%7.4f",mL1), col.axis=2, col=2)

mL1 <- -0.7; re2 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(x.iN, minL1 = mL1)))
cat("at x=", x.iN[im <- which.max(abs(re2))], "rel.Err.=", format(re2[im], digits = 3),"\n")
## at x= -0.6677246 rel.Err.= -6.48e-16
lines(x.iN, abs(re2), col=adjustcolor(4, 1/2), lwd=2)
abline(v = mL1, lwd=3,  col=adjustcolor(4, 2/3), lty=3)
axis(3, line=-1, at=mL1, sprintf("mL1 = %g",mL1), col.axis=4, col=4)
R <- x.iN >= mL1; stopifnot(re2[R] == re[iN][R])

mL1 <- -0.66; re3   <- asNumeric(relErrV(lg1pM.[iN], log1pmx(x.iN, minL1 = mL1)))
              re3.2 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(x.iN, minL1 = mL1, eps2 = .02)))
       print(max(abs(re3))) == # 3.7029e-16
             max(abs(re3.2)) #  TRUE !  of course, eps2 does not matter as long as it is far from our x-range !!
cat("at x=", x.iN[im <- which.max(abs(re3))], "rel.Err.=", format(re3[im], digits = 3),"\n")
## at x= -0.5634766 rel.Err.= -3.70e-16
lines(x.iN, abs(re3), col=adjustcolor(6, 1/3), lwd=2)
abline(v = mL1, lwd=3,  col=adjustcolor(6, 2/3), lty=3)
axis(3, line=-2, at=mL1, sprintf("mL1 = %g",mL1), col.axis=6, col=6)
lines(lowess(x.iN, abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
lines(lowess(x.iN, abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
lines(lowess(x.iN, abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
R <- x.iN > mL1;  stopifnot(re3[R] == re2[R])
## MM: for vignette, stop here?
## ---> since  max(|re4|) {below} == max(|re3|)

mL1 <- -0.64; re4 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(x.iN, minL1 = mL1, tol_logcf = 3e-16)))
max(abs(re4)) == max(abs(re3)) # TRUE ! i.e., are *not* better (in minimax-sense)
cat("at x=", x.iN[im <- which.max(abs(re4))], "rel.Err.=", format(re4[im], digits = 3),"\n")
## at x= -0.5634766 rel.Err.= -3.7e-16  -- of course unchanged from re3 (-0.66)
lines(x.iN, abs(re4), col=adjustcolor(7, 1/3), lwd=2)
abline(v = mL1, lwd=3,  col=adjustcolor(7, 2/3), lty=3)
axis(3, line=-2, at=mL1, sprintf("mL1 = %g",mL1), col.axis=7, col=7)
lines(lowess(x.iN, abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
lines(lowess(x.iN, abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
lines(lowess(x.iN, abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
lines(lowess(x.iN, abs(re4),    f=1/50), col=adjustcolor(7, 1/2), lwd=6)
R <- x.iN > mL1;  table(re4[R] == re3[R]) # no longer

mL1 <- -0.6; re5 <- asNumeric(relErrV(lg1pM.[iN], log1pmx(x.iN, minL1 = mL1)))
max(abs(re5)) == max(abs(re3)) # TRUE (as still < -0.564.. !)
cat("at x=", x.iN[im <- which.max(abs(re5))], "rel.Err.=", format(re5[im], digits = 3),"\n")
lines(x.iN, abs(re5), col=adjustcolor(8, 1/3), lwd=2)
abline(v = mL1, lwd=3, col=adjustcolor(8, 2/3), lty=3)
axis(3, line=-1, at=mL1, sprintf("mL1 = %g",mL1), col.axis=8, col=8)
lines(lowess(x.iN, abs(re[iN]), f=1/50), col=adjustcolor("gray", 1/2), lwd=6)
lines(lowess(x.iN, abs(re2),    f=1/50), col=adjustcolor(4, 1/2), lwd=6)
lines(lowess(x.iN, abs(re3),    f=1/50), col=adjustcolor(6, 1/2), lwd=6)
lines(lowess(x.iN, abs(re4),    f=1/50), col=adjustcolor(7, 1/2), lwd=6)
lines(lowess(x.iN, abs(re5),    f=1/50), col=adjustcolor(8, 1/2), lwd=6)
R <- x.iN > mL1;  table(re5[R] == re4[R]) # no longer equal
## -0.6  now is clearly too large, -0.7 was better  {MM: why the
showProc.time()


##--- separate: is  eps2 = 0.01  "optimal" / "good enough" ??

k. <- if(doExtras) 3 else 0
str(xS <- unique(sort(c(seq(-1/32, 1/32, by=2^-(12+k.)),
                        seq(-1,1, by=2^-(7+k.))/2^(4+k.))))) ## 257 or 3841 values
plot(xS, log1pmx(xS), type="l", col=2, main = "log1pmx(x)")
abline(h=0, v = c(-.01, .01), lwd=3, col=adjustcolor(7, 2/3), lty=3)
## nothing visible at +- .01

xSm <- mpfr(xS, 1024)
xSM <- mpfr(xS, 8192)
lg1pM. <- log1p(xSm) - xSm ## the direct formally using many bits ==> not much cancellation
summary(abs(asNumeric(relErrV(log1p(xSM) - xSM, lg1pM.))))# <= 4.1e-304 : 1024 is enough!
reS <- asNumeric(relErrV(lg1pM., log1pmx(xS)))
summary(reS) ## Lnx 64b gcc 11.2.1, 20210728:
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
## -2.625e-16 -5.246e-17  1.973e-18  1.211e-18  5.286e-17  2.688e-16
plot(xS, abs(reS), type="b", cex=1/2, log="y",
     main = "| relErr( log1pmx(<small x>) ) |  {via 'Rmpfr'}"
   , ylim = c(4e-17, max(abs(asNumeric(reS))))
   , panel.last= abline(h=0, v = c(-.01, .01), lwd=3, col=adjustcolor(7, 2/3), lty=3))
## no problem visible at  +/- eps2  =  +/- 0.01
stopifnot(abs(reS) < 1e-15, mean(abs(reS)) < 6e-16) #  seen 2.688e-16, 6.35e-17
showProc.time()

## in spite of the above; ... what's the approximation error of accurate log1pmx() ?
##                                        vvvvv  (or 1e-18 above ?)
p.Err <- function(eps2 = .01, tol_logcf = 1e-17, xM=xSm) {
    stopifnot(eps2 > 0, tol_logcf > 0,
              inherits(xM, "mpfr"), (prec <- .getPrec(xM)) >= 64)
    lg1pM <- log1p(xM) - xM # direct formula suffering from cancellation
    lg1pM.<- log1pmx(xM, tol_logcf=tol_logcf, eps2=eps2)
    sys.call()
    reM <- asNumeric(relErrV(lg1pM., lg1pM))
    plot(asNumeric(xM), abs(reM), type="b", cex=1/2, log="y", xlab = quote(x),
         ylim = quantile(abs(reM), c(1/20,1), names=FALSE),
         sub = paste("min(getPrec(x)) =", min(prec)),
         main = sprintf("relative error of R log1pmx(eps2=%g, tol_logcf=%g)", eps2, tol_logcf))
    abline(h=0, v = eps2*c(-1,1), lwd=3, col=adjustcolor(7, 2/3), lty=3)
    invisible(reM)
}

### --- for tol_logcf = 1e-17 ---------------------------------------------
reSM    <- p.Err()     #  aha! Here,  eps2 = .01  was definitely too large!
reSM001 <- p.Err(.001) # still too large
xM2 <- mpfr(seq(-1,1, length=1+1024)/2^8, 1024)
reSM3_4 <- p.Err(3e-4,    xM=xM2)# still too large
reSM294 <- p.Err(2.9e-4,  xM=xM2)# still
reSM2854<- p.Err(2.85e-4, xM=xM2)# still (and practically identical to 2.9e-4)
## cbind(reSM294, reSM2854)[reSM294 != reSM2854 , ]
##            reSM294      reSM2854
## [1,]  2.572818e-36  2.988004e-46
## [2,] -2.565617e-36 -2.977915e-46
reSM284 <- p.Err(2.8e-4,  xM=xM2)# now too small -- bad!
reSM274 <- p.Err(2.7e-4,  xM=xM2)#   (ditto) too small
showProc.time()

### --- for tol_logcf = 1e-14 == R's default!
reSM     <- p.Err(       tol = 1e-14) ## eps2=.01 is *still* too large
reSM002  <- p.Err(.002,  tol = 1e-14) # still too large
reSM0017 <- p.Err(.0017, tol = 1e-14) # still ..
reSM00165<- p.Err(.00165,tol = 1e-14) # still..
reSM00163<- p.Err(.00163,tol = 1e-14) # still..  <<<<<<<<<<<<< SHOULD WE CHANGE R TO USE THIS?
reSM00162<- p.Err(.00162,tol = 1e-14) # already (very barely!) too small
reSM0016 <- p.Err(.0016, tol = 1e-14) # already (barely!) too small
showProc.time()
if(do.pdf) dev.off()


summary(warnings())
