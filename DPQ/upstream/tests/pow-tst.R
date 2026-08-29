require("DPQ")
require("Rmpfr")

source(system.file(package="DPQ", "extraR", "relErr-plots.R", mustWork=TRUE))
##-> drawEps.h() , mtextVersion(), myPlatform() ...
(doExtras <- DPQ:::doExtras())
(noLdbl <- (.Machine$sizeof.longdouble <= 8)) ## TRUE when --disable-long-double
options(width = 100, nwarnings = 1e5)

k1 <- 44300 + (1:800)
p1 <- (63/64)
x <- p1; stopifnot(  (1-x)-1 == -x , 1-x == 2^-6 )
px <- p1 ^ k1
isNorm <- px >= 2^-1022 # is "normal", i.e., not subnormal
k1 <- k1[isNorm]
px <- px[isNorm]
dput(px,, "digits17")
dput(px,, "hex")


x <- p2 <- 1023/1024
stopifnot(  (1-x)-1 == -x , 1-x == 2^-10 )

k2 <- 700000 + (1:1000)
dput(p2^ k2, control = "digits17") ## last 2-3 digits printed *differ*
dput(p2^ k2, control = "hex")

x <- p3 <- (1- 2^-14)
stopifnot(  (1-x)-1 == -x , 1-x == 2^-14 )
log(1e-300)/ log(p3) # 11317321

k3 <- 11321012 + (1:1000)
dput(p3^ k3, control = "digits17") ## last 2-3 digits printed *differ*
dput(p3^ k3, control = "hex")

require(sfsmisc) # for relErrV(), eaxis()
## relative error of  x ^ k  computation - comparing with truth = <x_mpfr_{precBits])
relEP <- function(x, k, precBits, # see below about  precBits needed
                  mp = mpfr(x, precBits)^k,
                  FN = function(x,k) x^k) {
    stopifnot(is.numeric(x), inherits(mp, "mpfr"))
    asNumeric(relErrV(mp, FN(x,k)))
}

(myOS <- myPlatform())
if(!dev.interactive(TRUE)) pdf(paste0("pow-tst_", myOS, ".pdf"),
                               width = 9, height=5)

stopifnot(.Machine$double.xmin == 2^-1022) # will use  2^-1022

re1. <- relEP(p1, k1, precBits = 2^12) # high precision
re1  <- relEP(p1, k1, precBits = 256)
stopifnot(re1. == re1) # ===> 256 bit suffice here
summary(re1)
summary(abs(re1))

re1pow <- relEP(p1, k1, precBits = 256, FN = pow)
re1pdi <- relEP(p1, k1, precBits = 256, FN = pow_di)
stopifnot(identical(re1pow, re1pdi))# because pow() smartly uses pow_di() here
## however
all.equal(re1, re1pdi, tolerance = 0) # 2201.618 !!
summary(re1pdi) # negatively biased !
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
## -9.179e-14 -9.106e-14 -9.070e-14 -7.740e-14 -9.025e-14  0.000e+00
plot(k1, re1pdi, type="l"); abline(v = min(k1[p1^k1 < 2^-1022]), lty=3)

re1pow <- relEP(p1, k1, precBits = 256, FN = .pow)
summary(re1pow)
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.  (Linux Fedora 38):
## -1.070e-16 -3.082e-17  0.000e+00  7.063e-19  3.507e-17  1.075e-16

## Explore more about pow_di()'s "flaw"
k1. <- c(2:64, as.integer(unique(round(lseq(66, 44980, length = 1000)))))
re1.pdi <- relEP(p1, k1., precBits = 256, FN = pow_di)
plot(k1., re1.pdi, type="l", log="x", xaxt="n",
     main = "rel.error {wrt MPFR} of (63/64)^k -- using R_pow_di()")
abline(h=0, lty=3); eaxis(1, sub10=3)

plot(k1., re1.pdi, type="o", cex=1/4, log="x", ylim = c(-50, 2)*2^-53, xaxt="n",
     main = "rel.error of (63/64)^k -- R_pow_di() -- zoomed in")
abline(h=(-2:2)*2^-53, lty=c(3,3,1), col=adjustcolor(1, 1/2)); eaxis(1, sub10=3)

all.equal(p1^k1., pow(p1, k1., FALSE), tolerance=0) # TRUE (everywhere?)

## comparison:  R's  x^y i.e. R_pow() is *much* better -- starting already at ca k >= 50 (for this x=p1):
plot(k1., abs(re1.pdi), type="l", log="xy", yaxt="n", ylim = c(1e-17, max(abs(re1.pdi))),
     xlab = quote(k), main = "|rel.error| {wrt MPFR} of"~ (frac(63,64))^k)
lines(k1., abs(relEP(p1, k1., 256)), col = adjustcolor(2, 1/2))
legend("top", paste("using", c("R_pow_di(x,k)", "x ^ k")), lwd = 1,
       col = c(palette()[[1]], adjustcolor(2, 1/2)), bty="n")
eaxis(2); drawEps.h(); abline(h = 8e-18, col=adjustcolor(4, 1/2), lwd=2)
rug(k1.[k1. <= 50], col="gray30")

plot(k1, re1, type="l", main = "rel.error {wrt MPFR} of  (63/64)^k")
mtextVersion()
abline(v = min(k1[p1^k1 < 2^-1022]), lty=3)
ll <- c(1,1,2,1,1)
abline(h=(-2:2)*2^-53, lty=3-ll, lwd=ll, col=adjustcolor(ll, 1/2))
## looks perfect (in Linux; *not* windows

## ditto (pow_di exploration) for an even more critical p
p1 <- 0.999 # not exactly representable
print(p1, digits=19)
k <- c(2:50, as.integer(lseq(55, 700000, length=1000)))
pk <- p1^k
tail(cbind(k, pk))
plot(pk ~ k, type = "l", log="xy")
re.pdi <- relEP(p1, k, precBits = 256, FN = pow_di)

plot(k, abs(re.pdi), type="l", log="xy", xaxt="n", yaxt="n", ylim = c(1e-17, max(abs(re.pdi))),
     xlab = quote(k), main = substitute(abs(rel.error) ~~ "of" ~~ P^k, list(P = format(p1))))
lines(k, abs(relEP(p1, k, 256)), col = adjustcolor(2, 1/2))
legend("top", paste("using", c("R_pow_di(x,k)", "x ^ k")), lwd = 1,
       col = c(palette()[[1]], adjustcolor(2, 1/2)), bty="n")
eaxis(1, sub=3); eaxis(2); drawEps.h(); abline(h = 8e-18, col=adjustcolor(4, 1/2), lwd=2)
rug(k[k <= 50], col="gray30")
axis(1, at=2:6)
mtextVersion()

## finaly: pow_di exploration an "easy" p:
p1 <- 0.001 # not exactly representable
k <- 2:105
pk <- p1^k
tail(cbind(k, pk)) # subnormals ..
plot(pk ~ k, type = "l", log="xy")
re.pdi <- relEP(p1, k, precBits = 256, FN = pow_di)

plot(k, abs(re.pdi), type="l", log="xy", xaxt="n", yaxt="n", ylim = c(1e-17, max(abs(re.pdi))),
     xlab = quote(k), main = substitute(abs(rel.error) ~~ "of" ~~ P^k, list(P = format(p1))))
lines(k, abs(relEP(p1, k, 256)), col = adjustcolor(2, 1/2))
legend("top", paste("using", c("R_pow_di(x,k)", "x ^ k")), lwd = 1,
       col = c(palette()[[1]], adjustcolor(2, 1/2)), bty="n")
eaxis(1, sub=3); eaxis(2); drawEps.h(); abline(h = 8e-18, col=adjustcolor(4, 1/2), lwd=2)
rug(k, col="gray30")
axis(1, at=2:6)
mtextVersion()


re2  <- relEP(p2, k2, 2^ 9)
re2. <- relEP(p2, k2, 2^14)
stopifnot(re2. == re2)
summary(re2)
summary(abs(re2))

str(rat2 <- c(MASS:::.rat(p2) $ rat))
plot(k2, re2, type="l", xlab = quote(k), main = substitute("rel.Err of " ~ (N/Z)^k,
                                                           list(N=rat2[1], Z=rat2[2])),
     ylim = range(0, re2))
if(any((px <- p2^k2) < 2^-1022)) abline(v = min(k2[px < 2^-1022]), lty=3)
abline(h=(-2:2)*2^-53, lty=3-ll, lwd=ll, col=adjustcolor(ll, 1/2))

re3 <- relEP(p3, k3, prec = 256)
summary(re3)
str(rat3 <- c(MASS:::.rat(p3, max.denominator=1e5) $ rat))

plot(k3, re3, type="l", xlab = quote(k),  ylim = range(0, re3),
     main = substitute("rel.Err of " ~ (N/Z)^k, list(N=rat3[1], Z=rat3[2])))
if(any((px <- p3^k3) < 2^-1022)) abline(v = min(k3[px < 2^-1022]), lty=3)
abline(h=(-2:2)*2^-53, lty=3-ll, lwd=ll, col=adjustcolor(ll, 1/2))

##===> all these are perfect on Linux Fedora 36 & 38;
## but increasingly (p1 --> p2 --> p3)  worse on Windows

## even more extreme
e2 <- -(10:21)
 ps <- 1 - 2 ^ e2
(pM <- 1 - mpfr(2, 256) ^ e2) # 12 'mpfr' numbers of precision  256   bits
##  0.9990234375  0.99951171875 ... 0.999999523162841796875
stopifnot( (1 - ps) -1 == -ps, ps == pM) ## <==> all ps are exactly representable 1 - 2^{-m}

print(pM, scientific = FALSE, drop0trailing = TRUE)


## The range of interesting k -- such that
## ps^k  does *not* underflow (in double prec) but *is* ~= 10^{-300}:
(ks <- as.integer(round( log(1e-300)/ log(ps) )))
## [1]    707009    1414363    2829071    5658488   11317321  22634987  45270320
## [8]  90540985  181082315  362164975  724330295 1448660934
(k0 <- as.integer(signif(ks, 4)))

cbind(psN <- setNames(seq_along(ps), paste0("1-2^",e2)))
reL  <- lapply(psN, function(i) relEP(ps[i], k0[i] + 0:999, precBits = 2^8))
## can use names(reL) to reconstruct ps exactly:
stopifnot(identical(ps, vapply(lapply(names(reL), str2lang), eval, 0.9)))
if(doExtras) { ## much higher precision -- just to show it's *not* needed
    print(system.time(
        reLL <- lapply(psN, function(i) relEP(ps[i], k0[i] + 0:999, precBits = 2^ 14))
    )) # 2.407 sec
    stopifnot(all.equal(reL, reLL, tolerance = 0)) # exactly the same
}

absreL <- lapply(reL, abs)
t(sapply(reL,    summary))
t(sapply(absreL, summary)) # very nice on Linux:
##                 Min.      1st Qu.       Median         Mean      3rd Qu.         Max.
## 1-2^-10 4.698362e-20 1.991498e-17 4.117806e-17 4.292841e-17 6.362789e-17 1.085069e-16
## 1-2^-11 2.644379e-20 2.143401e-17 4.148793e-17 4.352628e-17 6.504935e-17 1.048650e-16
## 1-2^-12 5.185148e-21 2.465148e-17 4.910034e-17 4.719418e-17 6.823981e-17 1.001215e-16
## 1-2^-13 1.732619e-19 2.110998e-17 4.188601e-17 4.175295e-17 6.264898e-17 8.759273e-17
## 1-2^-14 1.028253e-19 2.692311e-17 5.463435e-17 5.234885e-17 7.736664e-17 1.027465e-16
## 1-2^-15 2.173905e-20 1.919186e-17 3.737047e-17 3.683430e-17 5.442273e-17 7.258842e-17
## 1-2^-16 7.547357e-20 1.964997e-17 4.034049e-17 4.078038e-17 6.051778e-17 8.354223e-17
## 1-2^-17 6.468787e-20 2.113894e-17 4.282366e-17 4.184208e-17 6.196293e-17 8.264685e-17
## 1-2^-18 3.053764e-20 2.343907e-17 4.581431e-17 4.515700e-17 6.764085e-17 8.893000e-17
## 1-2^-19 3.693339e-20 1.922513e-17 4.318434e-17 4.283787e-17 6.588232e-17 8.900600e-17
## 1-2^-20 5.398257e-20 2.190955e-17 4.092451e-17 4.095826e-17 6.083059e-17 8.062014e-17
## 1-2^-21 4.708387e-20 2.429069e-17 4.717160e-17 4.758539e-17 7.079720e-17 9.742768e-17

d.absre <- as.data.frame(absreL, optional=TRUE)
stopifnot(is.numeric(m.absrel <- as.matrix(d.absre)))
if(doExtras) { # not particularly interesting ..
    matplot(0:999, m.absrel, type="l", log="y", ylim = c(1e-17, max(m.absrel)))
    drawEps.h(lty=2)
                                        #
    ## smoothed |rel.errors|
    s.absre <- d.absre
    s.absre[] <- lapply(absreL, function(y) lowess(y, f = 1/20)$y) # robust smoothing
    matlines(0:999, s.absre, lwd=3, col = adjustcolor(1:12, 1/2))
}

## different "log-scale"  k sets for different ps / ks:
k2p   <- lapply(setNames(ks, names(psN)), \(k)
                as.integer(2^seq(7, floor(log2(k)), by = 1/8)))
stopifnot( sapply(k2p, is.integer) ) # -> so pow_di() / __powi() .. could be used
reL2  <- lapply(psN, function(i) relEP(ps[i], k2p[[i]], precBits = 2^8))
absre2 <- lapply(reL2, abs)

K <- length(ks) # 12 here

plot(k2p[[K]], reL2[[K]], type = "l", log = "x", col=K,
     xlab = quote(k), ylab = "relErr{ p ^ k }",
     ylim = range((-1:1)*2^-53, unlist(reL2)), xaxt="n"); eaxis(1)
for(i in rev(seq_along(ks))[-1L])
    lines(k2p[[i]], reL2[[i]], col=i)
abline(h = (-2:2)*2^-53, lty=c(3,3,2,3,3), lwd=c(1,1,2,1,1), col=adjustcolor(c(1,1,2,1,1), 1/2))

plot(k2p[[K]], absre2[[K]], type = "l", log = "xy", col=K,
     xlab = quote(k), ylab = "|relE|", main = expression(abs(rel.Err( p ^ k ))),
     ylim = range(2^-53, pmax(8e-18, unlist(absre2))), xaxt="n", yaxt="n"); eaxis(1); eaxis(2)
(pExpr <- as.call(c(quote(list), lapply(names(k2p), str2lang))))
mtext(substitute(p == group("{", EE, "}"), list(EE = pExpr)), cex= 2/3)
for(i in rev(seq_along(ks))[-1L])
    lines(k2p[[i]], absre2[[i]], col=i)
drawEps.h(lty=3, lwd=2)
## looks perfect (in Linux; *not* windows

t(sapply(absre2, summary))
summary(are2 <- unlist(absre2))

if(.Platform$OS.type == "unix") { # && !noLdbl
    stopifnot(m.absrel < 2.23e-16, # see max(.) == 1.085e-16 < 2^-53 = 1.11022e-16
              are2     < 2.23e-16) # see max(.) == 1.032e-16
}

###---- x ^ y -- with *non-integer* y ---------------

## keep x = ps  { above (1-2^-k) }

yp <- lapply(setNames(ks, names(psN)), \(k)
             2^seq(7, floor(log2(k)), by = 1/8) - 1/4) # -1/4: sure *not* be int

reLd  <- lapply(psN, function(i) relEP(ps[i], yp[[i]], precBits = 2^8))
absred <- lapply(reLd, abs)

K <- length(ks) # 12 here

plot(yp[[K]], reLd[[K]], type = "l", log = "x", col=K,
     xlab = quote(k), ylab = "relErr{ p ^ k }",
     ylim = range((-1:1)*2^-53, unlist(reLd)), xaxt="n"); eaxis(1)
for(i in rev(seq_along(ks))[-1L])
    lines(yp[[i]], reLd[[i]], col=i)
abline(h = (-2:2)*2^-53, lty=c(3,3,2,3,3), lwd=c(1,1,2,1,1), col=adjustcolor(c(1,1,2,1,1), 1/2))

palROBG <- colorRampPalette(c("red", "darkorange2", "blue", "seagreen"), space = "Lab")
palette(adjustcolor(palROBG(12), 3/4))
## pmax(.) for 'y': to "show" the '0' also
plot(yp[[K]], pmax(1e-20, absred[[K]]), type = "l", log = "xy", col=K,
     xlab = quote(y), ylab = "|relE|", main = quote(list(abs(rel.Err( x ^ y )), y~"*not* integer")),
     ylim = range(2^-53, pmax(8e-18, unlist(absred))), xaxt="n"); eaxis(1, sub10=3)
mtext(substitute(x == group("{", EE, "}"), list(EE = pExpr)), cex= 2/3)
for(i in rev(seq_along(ks))[-1L])
    lines(yp[[i]], absred[[i]], col=i)
abline(h = c(1,2,4)*2^-53, lty=3, lwd=2, col=adjustcolor(2, 1/2))
axis(4, las=2, line=-1, at=c(1,2,4)*2^-53, labels = expression(2^-53, 2^-52, 2^-51),
     col.axis=adjustcolor(2, 1/2), col=NA, col.ticks=NA)
legend("bottomright", legend=do.call(expression, as.list(pExpr[-1])),
       col=palette(), lwd=3, cex=2/3)
palette("default")# revert
## looks perfect (in Linux; *not* Windows)

t(sapply(absred, summary))
c(summary(ared <- unlist(absred)))

if(.Platform$OS.type == "unix") { # && !noLdbl
    stopifnot(ared < 2.23e-16) # see max(.) == 1.033987e-16 < 2^-53
}

