### ----------------- Testing  expm1x() { := e^x - 1 - x , numerically stable} ----
require(DPQ)
require(sfsmisc)


##------------------------ Testing *accuracy* via  Rmpfr {R <--> GNU MPFR C-library} ----
for(pkg in c("Rmpfr"))
    if(!requireNamespace(pkg)) {
        cat("no CRAN package", sQuote(pkg), " ---> no tests here.\n")
        q("no")
    }
require("Rmpfr")

#-------------------------------------------------------

do.pdf <- TRUE
do.pdf <- !dev.interactive(orNone = TRUE)
do.pdf
if(do.pdf) endPdf <- if(interactive()) sfsmisc::end.pdf else dev.off

## relErrV <- sfsmisc::relErrV
## eaxis   <- sfsmisc::eaxis


## direct formula - may be really "bad" :
expm1x.0 <- function(x) exp(x) -1 - x
## less direct formula - improved (but still not universally ok):
expm1x.1 <- function(x) expm1(x)  - x

## a symmetric set of negative and positive
x <- unique(c(2^-seq(-3/8, 54, by = 1/8), seq(7/8, 3, by = 1/128)))
x <- x0 <- sort(c(-x, 0, x)) # negative *and* positive

## Mathematically,  expm1x() = exp(x) - 1 - x  >= 0  (and == 0 only at x=0):
em1x <- expm1x(x)
stopifnot(em1x >= 0, identical(x == 0, em1x == 0))

xxTrue <- expm1x.1(mpfr(x, 1024))
relE  <- asNumeric(relErrV(xxTrue, em1x))
relE1 <- asNumeric(relErrV(xxTrue, expm1(x)-x))

if(do.pdf) pdf("expm1x_relE-1.pdf")
plot(x,  abs(relE), log="y", type = "b", cex=1/2, ylim = c(2e-17, max(abs(relE))))
lines(x, abs(relE1), col = adjustcolor(2, 1/2), lwd=3)
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50)))
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)

plot (abs(x), abs(relE), log="xy", type = "b", cex=1/2, ylim = c(2e-17, max(abs(relE))))
lines(abs(x), abs(relE1), col = adjustcolor(2, 1/2), lwd=3)
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50)))
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)

iL <- which(abs(relE) > 2^-52)
print(cbind(x = x[iL], relE = relE[iL] , relE1 = relE1[iL]), digits = 4)
##          x       relE      relE1
##  9.686e-16 -2.255e-16 -1.591e-01
##  5.611e-12  2.321e-16  1.572e-06
##  3.250e-08 -2.773e-16 -5.176e-09
##  3.328e-05 -2.512e-16  2.289e-12
##  1.704e-02 -2.318e-16  3.110e-15
##  1.211e+00  2.221e-16  2.221e-16
stopifnot(abs(relE) < 5e-16, sum(abs(relE) > 2^-52) <= 6 + 10)# the above + "some slack"


### -------------------- Relative error of Taylor series approximations -----------
##
if(do.pdf && .Device == "pdf") { dev.off(); pdf("expm1x_relE_Taylor.pdf") }

twoP <- seq(-0.75, 54, by = 1/8)
x <- 2^-twoP
x <- sort(c(-x,x)) # negative *and* positive
e1xAll <- cbind(expm1x.0 = expm1x.0(x),
                expm1x.1 = expm1x.1(x),
                vapply(1:15, \(k) expm1xTser(x, k=k), x))
colnames(e1xAll)[-(1:2)] <- paste0("k=",1:15)

xM <- mpfr(x, 1024) # high accuracy (to push cancellation out of dbl.prec. range)
expm1xM <- expm1x.1(xM)
expm1x.relE <- asNumeric(relErrV(expm1xM, e1xAll))

pl.relEexpm1x <- function(x, relE, ind = TRUE, type = "l", las = 2, ...,
                          leg.x = "top", leg.y = NULL, leg.ncol = 4,
                          leg.cex = if(.Device == "pdf") 4/5 else 3/4)
{
    stopifnot(is.numeric(N <- ncol(relE)), N >= 2)
    legs <- c("exp(x) - 1 - x", "expm1(x) - x",
              paste0("expm1xTser(x, ", colnames(relE)[-(1:2)], ")"))
    stopifnot(N == length(legs))
    ## matplot() default has  col = 1:6, lty = 1:5
    col <- rep_len(1:6, N)[ind]
    lty <- rep_len(1:5, N)[ind]
    matplot(x, relE[, ind], type=type, las=las, xaxt = "n", ...,
            col = col, lty = lty, ylab = deparse1(substitute(relE)),
            main = expression("Accuracy (Relative Error) of" ~~~~~ e^x - 1 - x ~~ "Computations"))
    eaxis(1, sub10 = c(-2, 2))
    legend(leg.x, leg.y, legend = legs[ind], ncol=leg.ncol,
           cex = leg.cex, bty = "n", col = col, lty = lty)
}

## no x-log scale here at first:
pl.relEexpm1x(x, expm1x.relE, ylim = c(-1,1)*1000, leg.x = "topright", leg.ncol = 3) # non sense
rug(x)
pl.relEexpm1x(x, expm1x.relE, xlim = c(-1,1)*0.5, ylim = c(-1,1)*1000)
pl.relEexpm1x(x, expm1x.relE, xlim = c(-1,1)*0.1, ylim = c(-1,1)*1000) # still "non sense"

I  <- x > 0 # positive
pl.relEexpm1x(x[I], expm1x.relE[I,], log="x", ylim = c(-1,1)*100, leg.x = "topright", leg.ncol = 3)
N <- x < 0  # negative
matlines(-x[N], expm1x.relE[N,]) # not much changes in picture -- good ?!

## not showing the worst (= the first):
pl.relEexpm1x(x[I], expm1x.relE[I,], log="x", ylim = c(-1,1) /  2 , ind = -1)
                                        # see how Taylor k=1 and k=2 diverge
## zoom-in (ylim):
pl.relEexpm1x(x[I], expm1x.relE[I,], log="x", ylim = c(-1,1) /100 , ind = -1)
                                        # see how Taylor k=1,2,3,4  diverge
matlines(-x[N], expm1x.relE[N,-1], col = c(2:6,1), lty = c(2:5,1))
## with correct  colors -- shows different error *sign* of k=<odd> !

## Much more relevant: *relative* error (and log-scale) :
pl.relEexpm1x(x[I], abs(expm1x.relE[I,]), log="xy", yaxt="n", leg.cex = .9) ; eaxis(2, nintLog = 20)
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50)))
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)
if("nintLog" %in% names(formals(grid))) grid(nintLog = 20) else grid() # in the future R >= 4.5
matlines(-x[N], expm1x.relE[N,]) # onle little changes in picture -- but visible for "large x" (x ~= 1) good!
## draw some vertical lines "at" cutpoints x[k]:
abl <- function(v, lty = 3, col = adjustcolor(1, 1/2), lwd = 2, ...) {
    ## abline(v = v, lty=lty, col=col, lwd=lwd, ...)
    uy <- par("usr")[3:4]
    segments(x0 = v, y0 = 10^uy[1], y1 = 2^-52, lty=lty, col=col, lwd=lwd, ...)
    if(getOption("scipen") > -2) { op <- options(scipen = -2); on.exit(op) }
    form <- function(v) sub("e-0", "e-", format(v, digits=3))
    axis(1, at=v, labels=form(v), padj = -3, cex = 0.75)
}
abl(v = 5.00e-16)# k = 1 <==> ord=2
abl(v = 4.40e-8)
abl(v = 2.02e-5) # k = 3 <==> ord=4
abl(v = 3.95e-4)
abl(v = 3.00e-3) # k = 5 <==> ord=6
## abl(v = 0.0111)

if(do.pdf && .Device == "pdf") dev.off()


###--------------------- Older Experiments for finding cutoffs and visualization --- originally in ../R/expm1x.R

(doExtras <- DPQ:::doExtras() && !grepl("valgrind", R.home())) # TRUE, typically when interactive
if(!doExtras) quit("no")

## if(doExtras) : ----------------------------------------------

x <- unique(c(2^-seq(-3/8, 54, by = 1/8),
              seq(7/8, 3, by = 1/128)))

## x <- 2^-seq(-2, 27, by = 1/4) # at least *one* < 4.4e-8
## x <- 2^-seq(-2, 12, by = 1/2) # for debugging
x <- x0 <- sort(c(-x, 0, x)) # negative *and* positive

## cutx = c(  4.4e-8, 0.10, 0.385, 1.1, 2)  # cutoff x[k]
## k     = c(2,      9,   12,   17)
## r <- ax <- abs(x)
## ## vectorizing using findInterval():
## in.x <- findInterval(ax, c(0, cutx, Inf), all.inside = TRUE)

xxx <- expm1x(x)
stopifnot(identical(x == 0, xxx == 0))

if(do.pdf) pdf("expm1x_part-2.pdf")

plot (x, xxx, type='b', log="y")
lines(x, expm1(x)-x, col = adjustcolor(2, 1/2), lwd = 3) ## should nicely cover ..
lines(x, exp(x)-1-x, col = adjustcolor(4, 1/4), lwd = 5) ## should nicely cover ..
cuts <- c(4.4e-8, 0.10, 0.385, 1.1, 2)[-1] # *not* drawing 4.4e-8
v <- c(-rev(cuts), 0, cuts); stopifnot(!is.unsorted(v))
abline(v = v, lty = 3, col=adjustcolor("gray20", 1/2))

stopifnot(diff(xxx[x <= 0]) <= 0)
stopifnot(diff(xxx[x >= 0]) >= 0)

xxTrue <- expm1x.1(mpfr(x, 1024))
relE  <- asNumeric(relErrV(xxTrue, xxx))
relE1 <- asNumeric(relErrV(xxTrue, expm1(x)-x))

plot(x,  abs(relE), log="y", type = "b", cex=1/2, ylim = c(2e-17, max(abs(relE))))
lines(x, abs(relE1), col = adjustcolor(2, 1/2), lwd=3)
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50)))
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)
abline(v = v, lty = 3, col=adjustcolor("gray20", 1/2))

plot(abs(x), abs(relE), log="xy", type = "b", cex=1/2, ylim = c(2e-17, max(abs(relE))))
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50)))
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)
lines(abs(x), abs(relE1), col = adjustcolor(2, 1/2), lwd=3)
abline(v = v, lty = 3, col=adjustcolor("gray20", 1/2))

iL <- which(abs(relE) > 2^-52)
cbind(x = x[iL], relE = relE[iL] , relE1 = relE1[iL])
##          x       relE      relE1
##  9.686e-16 -2.255e-16 -1.591e-01
##  5.611e-12  2.321e-16  1.572e-06
##  3.250e-08 -2.773e-16 -5.176e-09
##  3.328e-05 -2.512e-16  2.289e-12
##  1.704e-02 -2.318e-16  3.110e-15
##  1.211e+00  2.221e-16  2.221e-16


expm1xBnds <- DPQ ::: expm1xBnds ## from ../R/expm1x.R
##                                       ^^^^^^^^^^^^^  see there, for the default double precision P = 52

### Using  P = 128  as for mpfr(x, 128) :
xBnd128 <- lapply(2:18, expm1xBnds, P = 128L)
dBnd128 <- do.call(rbind, lapply(xBnd128, data.frame)) # <- nice data.frame from the  uniroot() lists
dBnd128 <- cbind(ord = 2:18, within(dBnd128, { rm(init.it); x. <- exp(root) }))
if(getOption("digits") > 5) options(digits = 5) # op has 'digits' already
dBnd128
##    ord    root     f.root iter estim.prec        x.
## 1    2 -89.416  0.000e+00    1  9.442e+01 1.469e-39
## 2    3 -44.361 -7.105e-15    2  5.000e-08 5.421e-20
## 3    4 -29.208 -3.553e-15    2  5.000e-08 2.066e-13
## 4    5 -21.559  0.000e+00    1  7.844e+01 4.333e-10
## 5    6 -16.926  0.000e+00    1  8.307e+01 4.459e-08
## 6    7 -13.806  1.776e-15    2  5.000e-08 1.009e-06
## 7    8 -11.556  1.776e-15    2  5.000e-08 9.580e-06
## 8    9  -9.851  0.000e+00    1  9.015e+01 5.267e-05
## 9   10  -8.513  1.776e-15    2  5.000e-08 2.009e-04
## 10  11  -7.431 -8.882e-16    2  5.000e-08 5.925e-04
## 11  12  -6.538  0.000e+00    1  9.346e+01 1.448e-03
## 12  13  -5.786 -8.882e-16    2  5.000e-08 3.071e-03
## 13  14  -5.143  0.000e+00    1  9.486e+01 5.838e-03
## 14  15  -4.587 -8.882e-16    2  5.000e-08 1.018e-02
## 15  16  -4.101  1.776e-15    2  5.000e-08 1.655e-02
## 16  17  -3.672  0.000e+00    1  9.633e+01 2.544e-02
## 17  18  -3.289 -8.882e-16    2  5.000e-08 3.730e-02

### Using  P = 1024  as for mpfr(x, 1024) :
xBnd1024 <- lapply(2:18, expm1xBnds, P = 1024L)
dBnd1024 <- do.call(rbind, lapply(xBnd1024, data.frame))
dBnd1024 <- cbind(ord = 2:18, within(dBnd1024, { rm(init.it); x. <- exp(root) }))
dBnd1024
##    ord     root      f.root iter estim.prec          x.
## 1    2 -710.476 -1.1369e-13   12 5.0000e-08 2.7813e-309
## 2    3 -354.891  0.0000e+00    9 3.7264e+02 7.4583e-155
## 3    4 -236.228 -2.8422e-14   10 5.0000e-08 2.5555e-103
## 4    5 -176.824  0.0000e+00    8 1.8817e+02  1.6074e-77
## 5    6 -141.138  0.0000e+00    7 1.4929e+02  5.0663e-62
## 6    7 -117.316  0.0000e+00    6 1.2387e+02  1.1227e-51
## 7    8 -100.279  0.0000e+00    2 1.0533e+02  2.8153e-44
## 8    9  -87.484  0.0000e+00    1 9.2484e+01  1.0144e-38
## 9   10  -77.519  0.0000e+00    1 8.2519e+01  2.1567e-34
## 10  11  -69.537  0.0000e+00    1 7.4537e+01  6.3154e-31
## 11  12  -62.998 -7.1054e-15    2 5.0000e-08  4.3701e-28
## 12  13  -57.541 -7.1054e-15    2 5.0000e-08  1.0242e-25
## 13  14  -52.917  0.0000e+00    1 5.7917e+01  1.0432e-23
## 14  15  -48.949 -7.1054e-15    2 5.0000e-08  5.5177e-22
## 15  16  -45.505  0.0000e+00    1 5.4495e+01  1.7274e-20
## 16  17  -42.488  7.1054e-15    2 5.0000e-08  3.5302e-19
## 17  18  -39.82   0.000e+00     1  6.018e+01  5.077e-18

## Much more relevant: rel.error:

if(.Device == "pdf") dev.off()
if(do.pdf) pdf.do("expm1x_abs_relErr_L.pdf", paper = "a4r") # A4 rotated : horizontal

## zoom in larger x  {after recomputing with more x there}
## above had  2^twoP; twoP <- seq(-0.75, 54, by = 1/8)
xL <- 2^-seq(-1, 13, by = 1/64)
expm1xAllxL <- cbind(expm1x.0 = expm1x.0(xL), expm1x.1 = expm1x.1(xL),
                     vapply(1:16, \(k) expm1xTser(xL, k=k), xL))
colnames(expm1xAllxL)[-(1:2)] <- paste0("k=",1:16)

expm1x.relExL <- asNumeric(relErrV(expm1x.1(mpfr(xL, 1024)), expm1xAllxL))


pl.relEexpm1x(xL, abs(expm1x.relExL), log="xy", yaxt="n", leg.cex = .8) ; eaxis(2, nintLog = 20)
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50)))
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)
if("nintLog" %in% names(formals(grid))) grid(nintLog = 20) else grid() # in the future R >= 4.5
axis(1, at = c(0.2, 0.5, 2))
## draw some vertical lines "at" cutpoints x[k]:
## abl(v = 2.02e-5) # k = 3 <==> ord=4
abl(v = 3.95e-4)
abl(v = 3.00e-3)# k = 5 <==> ord=6
abl(v = 0.0111) # k = 6
abl(v = 0.031)
abl(v = 0.06) # k = 8
abl(v = 0.10)
abl(v = 0.185)# k = 10

## abl(v = 0.27)

## zoom around x = 1:
str(x1 <- seq(3/16, 1.5, by = 1/512))
expm1xAllx1 <- cbind(expm1x.0 = expm1x.0(x1), expm1x.1 = expm1x.1(x1),
                     vapply(10:18, \(k) expm1xTser(x1, k=k), x1))
colnames(expm1xAllx1)[-(1:2)] <- paste0("k=",10:18)

expm1x.relEx1 <- asNumeric(relErrV(expm1x.1(mpfr(x1, 1024)), expm1xAllx1))

if(.Device == "pdf") endPdf()
if(do.pdf) pdf.do("expm1x_abs_relErr_x=1.pdf", paper = "a4r") # A4 rotated : horizontal

pl.relEexpm1x(x1, abs(expm1x.relEx1), log="xy", yaxt="n", ylim = c(1e-17, 1e-9), leg.cex = .8) ; eaxis(2)
abline(h = 2^-(52:51), lty=3, col=paste("gray",c(20,50))); mtext("around x ~= 1")
axis(4, at=2^-52, labels = quote(2^{-52}), col.ticks=NA, las=2, cex.axis=3/4, hadj = 1/2)
if("nintLog" %in% names(formals(grid))) grid(nintLog = 20) else grid() # in the future R >= 4.5
axis(1, at = c(0.2, 0.5, 2))
## abl(v = 0.0111) # k = 6
## abl(v = 0.031)
## abl(v = 0.06) # k = 8
## abl(v = 0.10)
abl(v = 0.185)# k = 10
abl(v = 0.27) # k = 11
abl(v = 0.385)
abl(v = 0.49) # k = 13
abl(v = 0.63)
abl(v = 0.79) # k = 15
abl(v = 0.97) # k = 16
abl(v = 1.12) # k = 17
abl(v = 1.25) # k = 18

if(.Device == "pdf") endPdf()
