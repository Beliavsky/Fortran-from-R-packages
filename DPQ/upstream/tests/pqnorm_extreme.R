#### Extreme-tail pnorm() and qnorm() approximations
#### ===============================================
### notably  pnormAsymp()  and  qnormAsymp()
###          ~~~~~~~~~~~~       ~~~~~~~~~~~~
## partly inspired by MM's NUMERICS/dpq-functions/qnorm-asymptotic.R
## and also ../man/pnormAsymp.Rd
##          ../man/qnormAsymp.Rd

library(DPQ)
library(sfsmisc)# {we import it in DPQ}

p.ver <- function() mtext(R.version.string, cex=3/4, adj=1)
readUser <- function(i, max.i) {
    if(interactive() && i != max.i) {
        cat("[Enter] to continue: ")
        cat(readLines(stdin(), n=1), "\n")
    }
}
# "bisque4" , "coral2" ..
p.epsC <- function(col = "coral2", alpha = 3/4, lwd1 = 3) {
    lg <- par("ylog")
    mult <- if(lg) c(1/2, 1, 2) else -2:2 # multiplicities
    is1 <- abs(mult) == 1
    abline(h = .Machine$double.eps * mult,
           col = adjustcolor(col, alpha),
           lty = ifelse(is1, 2, 5),
           lwd = ifelse(is1, lwd1, 1))
}

absP <- function(re) pmax(abs(re), 2e-17) # not zero, so log-scale "shows" it


if(!dev.interactive(orNone=TRUE)) pdf("pqnorm-extreme.pdf")

x <- sort(2 ^ c(seq(3, 10, by=1/16), seq(10.125, 17, by=1/8)))
str(x)
summary(x)

pn <- pnorm(x, lower.tail=FALSE, log.p=TRUE)
## look at the asymptotic approximations:
k.s <- 0:5; nks <- paste0("k=", k.s); k. <- setNames(k.s, nks)
pnAsy <- sapply(k., function(k) pnormAsymp(x, k=k, lower.tail=FALSE, log.p=TRUE))

relEpnAsy <- apply(pnAsy, 2, relErrV, target = pn)
tit <- "relErrV(pnormAsym(x, k=k, lower.tail=FALSE, log.p=TRUE)\n  vs  pnorm(x, ..)"
leg.k <- function(x, y=NULL, ...) legend(x,y, legend=paste0("k=",0:5), lwd=2, col=1:6, lty=1:6, bty="n", ...)
matpl <- function(x, y, ...) {
    matplot(x,y, type="l", lwd=2, log="x", main=tit, xaxt="n", ...)
    eaxis(1, sub10=3); grid()
}

matpl(x, relEpnAsy, ylim = c(-8,1)*5e-6 ); leg.k("center")
matpl(x, relEpnAsy, ylim = c(-1,1)*1e-8 ); leg.k("topright")
matpl(x, relEpnAsy, ylim = c(-1,1)*1e-10); leg.k("topright")
matpl(x, relEpnAsy, ylim = c(-1,1)*1e-12); leg.k("topright")
matpl(x, relEpnAsy, ylim = c(-1,1)*1e-14); leg.k("topright"); p.epsC()
matpl(x, relEpnAsy, ylim = c(-1,1)*1e-15); leg.k("topright"); p.epsC()
abline(v = 38581.371, col=adjustcolor(4, 1/2), lwd=3) # the Rmpfr maximal value (for "independent accurate pnorm")

## They seem to be of opposite sign - mostly (and apart from +/- epsC)
## *but*  k=5 (at the beginning; (and then +/- epsC)

## ?proof?  --- use Rmpfr-version pnormAsymp()  but there,
## .......  pnorm(x, log.p=TRUE) *underflows* to -Inf as soon as x > 38581.371

## Start with quantiles so we know the  "truth according to pnorm()" :
qs <- c(2^seq(0, 35, by=1/256), Inf) # => s >= 1.84  --> no NaN in xs_5 etc
## == the "true" qnorm() according to pnorm()

lp  <- pnorm(qs, lower.tail=FALSE, log.p=TRUE)
qnp <- qnorm(lp, lower.tail=FALSE, log.p=TRUE)
s <- -lp # =  -log(1 - Phi(qs))
summary(r <- sqrt(s))
## store in 'dat' including current "inaccuracy" of qnorm():
head(dat <- data.frame(l2q=log2(qs), qs, # lp,
                       s,  # r = sqrt(s), t = log(s),
                       relE_qn = relErrV(qs, qnp)))
				##^^^^^^  will depend much on R version
stopifnot(dat[nrow(dat), "relE_qn"] == 0)

## relative Error (log x)
plot(relE_qn ~ s, dat, log="x", type="l", col = 2); p.ver()

## |rel.Error| ==> log-log scale
plot(abs(relE_qn) ~ s, dat, log="xy", type="l", col = 2,
     axes=FALSE, main = "|relative Error| of  qnorm(-s, log.p=TRUE, lower=F)  wrt pnorm()")
eaxis(1); eaxis(2); p.ver(); p.epsC()

## now look at the asymptotic approximations:
k.s <- 0:5; nks <- paste0("k=", k.s); k. <- setNames(k.s, nks)
qnAsy <- sapply(k., function(ord) qnormAsymp(lp=lp, order=ord))
stopifnot(identical(qnAsy,
         sapply(k., function(ord) qnormAsymp(p=lp, lower.tail=FALSE, log.p=TRUE, order=ord))
))

relEAsy <- apply(qnAsy, 2, relErrV, target = qs)


matplot(sqrt(s), relEAsy, type="l", log="x", xlab = quote(r == sqrt(s)),
        main = "relative Error of qnormAsymp(.., k=*)"); p.ver()
legend("top", nks, col=1:6, lty=1:5, bty="n")

mp <- c(-1,1)
xL <- list(NULL, c(2, 1000),  c(2, 1e5), c(5, 1e5), c(10, 1e8), c(10, 1e10))
yL <- list(NULL,   mp * 0.1,  mp * 1e-4, mp * 1e-7, mp * 1e-11, mp * 2e-15)
stopifnot(length(xL) == length(yL))
for(j in seq_along(xL)) {
    matplot(sqrt(s), relEAsy, type="l", log="x", xaxt="n", lwd=2,
            xlab = quote(r == sqrt(s)),
            xlim = xL[[j]], ylim = yL[[j]],
            main = "relative Error of qnormAsymp(.., k=*)")
    eaxis(1, sub10=2) # not here: p.ver()
    legend("top", nks, col=1:6, lty=1:5, lwd=2, bty="n")
    readUser(j, length(xL))
}

for(j in seq_along(xL)) {
    if(!is.null(yli <- yL[[j]])) yli <- c(2e-17, 100*yli[2])
    if(!is.null(xli <- xL[[j]])) xli[2] <- 10*xli[2]
    cat("j=",j,"; rbind(xlim,ylim) = \n"); print(rbind(xlim=xli, ylim=yli))
    cols <- adjustcolor(1:6, 3/4)
    ## abs() does not show 0 {-> -Inf);  absP() will show them
    matplot(sqrt(s), abs(relEAsy), type="l", log="xy", xaxt="n", yaxt="n", col=cols, lwd=2,
            xlab = quote(r == sqrt(s)),
            xlim = xli, ylim = yli,
            main = "relative Error of qnormAsymp(.., k=*)")
    matlines(sqrt(s), absP(relEAsy), col=adjustcolor(cols, 1/2))
    eaxis(1, sub10=2) # not here: p.ver()
    eaxis(2, sub10=c(-3,2))
    legend("top", nks, col=cols, lwd=2, lty=1:5, bty="n")
    readUser(j, length(xL))
}

## Real tests :
table(relEAsy[s > 1e17,] * 2^52)
## -0.5    0    1
##   71 9434  125

table(relEAsy[s > 1e16,] * 2^52)
## -0.5     0     1     2     3     4     5
##   71 11550   311   135    78    18    17


## drop k = 0 from now on -------------  The following *is* platform dependent .. let's see

re1 <- abs(relEAsy[s > 95e6, "k=1"])
table(re1 * 2^52)
##    0  0.5    1
## 5264   84   93

plot(abs(relEAsy[,"k=2"]) ~ s, subset = s > 1e5 & s < 1e8, log="xy", type="l", xaxt="n")
eaxis(1); p.epsC()

re2 <- abs(relEAsy[s > 2e5, "k=2"])
table(re2 * 2^52)
##    0  0.5    1
## 6353   85  141

plot(abs(relEAsy[,"k=3"]) ~ s, subset = s > 500 & s < 1e5, log="xy", type="l", xaxt="n")
eaxis(1); p.epsC()

re3 <- abs(relEAsy[s > 4000, "k=3"])
table(re3 * 2^52)
##    0  0.5    1
## 6953  124  225

plot(abs(relEAsy[,"k=4"]) ~ s, subset = s > 300 & s < 2000, log="xy", type="l", axes=FALSE)
eaxis(1);eaxis(2); p.epsC()

re4 <- abs(relEAsy[s > 1500, "k=4"])
table(re4 * 2^52)
##    0  0.5    1
## 7154  145  185

plot(abs(relEAsy[,"k=5"]) ~ s, subset = s > 200 & s < 1000, log="xy", type="l", axes=FALSE)
eaxis(1);eaxis(2); p.epsC()
r. <- pretty(sqrt( par("xaxp")[1:2] ), 10)
(rlab <- as.expression(lapply(r., function(rr) substitute(R^2, list(R = rr+0)))))
axis(3, at=r.^2, labels = rlab, col = 4, col.axis = 4, mgp = c(1.25,.5,0))

re5 <- abs(relEAsy[s > 700, "k=5"])
table(re5 * 2^52)
##    0  0.5    1
## 7285  141  199

stopifnot(exprs = {
    abs(relEAsy[s > 1e17,]) <= 2^-52 # even for k=0
    re1 <= 2^-52
    re2 <= 2^-52
    re3 <= 2^-52
    re4 <= 2^-52
    re5 <= 2^-52
})


### R code from vignette source '../vignettes/qnorm-asymp.Rnw' ====================
###
### code chunk number 2: qnormLog-compute
#########################################
qs <- 2^seq( 0, 29, by=1/256) # => s >= 1.84
lp <- pnorm(qs, lower.tail=FALSE, log.p=TRUE)
s <- -lp # = -pnorm(..) = -log(1 - Phi(qs)) > 0
##
qnrm    <- qnorm (-s, lower.tail=FALSE, log.p=TRUE)
qnrm405 <- qnormR(-s, lower.tail=FALSE, log.p=TRUE, version= "4.0.x") # R <= 4.0.5
qnrm410 <- qnormR(-s, lower.tail=FALSE, log.p=TRUE, version= "2020-10-17")
Rver <- shortRversion()
if(getRversion() <= "4.0.5") { # our qnormR(.., version="4.0.x")
    cat(sprintf("%s, \"4.0.5\",\n   all.equal(*, tol=0): %s;  identical(): %s\n", Rver,
                all.equal(qnrm, qnrm405, tolerance=0), identical(qnrm, qnrm405)))
    stopifnot(all.equal(qnrm, qnrm405, tolerance = 1e-12))
} else if(getRversion() < "4.3") { # our qnormR(*, version="2020-10-17") matches:
    cat(sprintf("%s, \"4.1.0\",\n   all.equal(*, tol=0): %s;  identical(): %s\n", Rver,
                all.equal(qnrm, qnrm410, tolerance=0), identical(qnrm, qnrm410)))
    ## see TRUE twice, for R 4.2.2, Linux{x86_64-pc-linux-gnu}  *and*  Windows{x86_64-w64-mingw32/x64}
                                        # M1mac(aarch64-apple-darwin20, R 4.2.1 ptchd): 2.675587e-16
    stopifnot(all.equal(qnrm, qnrm410, tolerance = 1e-12))
} else { # R version >= 4.3.x
    qnrm43 <- qnormR(-s, lower.tail=FALSE, log.p=TRUE, version = "2022")
    cat(sprintf("%s, >= 4.3.x,\n   all.equal(*, tol=0): %s;  identical(): %s\n", Rver,
                all.equal(qnrm, qnrm43, tolerance=0), identical(qnrm, qnrm43)))
    rE6 <- qnorm(-1e6, log.p=TRUE)/-1414.2077829910174  - 1
    cat(sprintf("  rE(-1e6) = %g\n", rE6))
    if(abs(rE6) < 7e-16) # have R-devel with new 2022 code:
        stopifnot(all.equal(qnrm, qnrm43, tolerance = 1e-14))
}

source(system.file("extraR", "qnorm-asymp-utils.R", package="DPQ"))
if(!exists("r0", mode="numeric"))
    q("no")
if(!dev.interactive(orNone=TRUE)) {
    dev.off()
    pdf("pqnorm-qnormAsy2.pdf", height=8*sqrt(2), width=8) # ~ 'A4'
}

###################################################
### code chunk number 27: plot-qnormAsy2
###################################################
mult.fig(5, main = "qnormAsymp(*, k) approximations in the 5 cutpoint regions")
r0 <- c(27, 55, 109, 840, 36000, 6.4e8) # <-- cutoffs  <--> in ../R/norm_f.R
# use k =  5   4    3    2      1       0    e.g.  k = 0  good for r >= 6.4e8
for(ir in 2:length(r0))
    p.qnormAsy2(r0[ir], k = 5 +2-ir, cex.main = .90)

