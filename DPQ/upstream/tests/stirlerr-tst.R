#### Testing  stirlerr()
#### =================== {previous 2nd part of this, now -->>> ./bd0-tst.R <<<---
require(DPQ)
for(pkg in c("Rmpfr", "DPQmpfr"))
    if(!requireNamespace(pkg)) {
        cat("no CRAN package", sQuote(pkg), " ---> no tests here.\n")
        q("no")
    }
require("Rmpfr")

source(system.file(package="DPQ", "test-tools.R", mustWork=TRUE))
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS(), matplotB()
##_ options(conflicts.policy = list(depends.ok=TRUE, error=FALSE, warn=FALSE))
require(sfsmisc) # masking  'list_' *and* gmp's factorize(), is.whole()
##_ options(conflicts.policy = NULL)

## plot1cuts() , etc:  ---> ../inst/extraR/relErr-plots.R <<<<<<<
##			    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(system.file(package="DPQ", "extraR", "relErr-plots.R", mustWork=TRUE))

do.pdf <- TRUE # (manually)
do.pdf <- !dev.interactive(orNone = TRUE)
do.pdf
if(do.pdf) {
    pdf.options(width = 9, height = 6.5) # for all pdf plots {9/6.5 = 1.38 ; 4/3 < 1.38 < sqrt(2) [A4]
    pdf("stirlerr-tst.pdf")
}

showProc.time()
(doExtras <- DPQ:::doExtras())
(noLdbl <- (.Machine$sizeof.longdouble <= 8)) ## TRUE when --disable-long-double or 'M1mac' ..
(M.mac <- grepl("aarch64-apple", R.version$platform)) # Mac with M1, M2, ... proc

abs19 <- function(r) pmax(abs(r), 1e-19) # cut  |err| to positive {for log-plots}

options(width = 100, nwarnings = 1e5, warnPartialMatchArgs = FALSE)



##=== Really, dpois_raw()  and  dbinom_raw()  *both* use stirlerr(x)  for "all"  'x > 0'
##            ~~~~~~~~~~~       ~~~~~~~~~~~~             ===========              ===== !

## below, 6 "it's okay, but *far* from perfect:" ===>  need more terms in stirlerr() [
## April 20: MM added more terms up to S10; 2024-01: up to S12 ..helps a little only
x <- lseq(1/16, 6, length=2048)
system.time(stM <- DPQmpfr::stirlerrM(Rmpfr::mpfr(x,2048))) # 1.7 sec elapsed
plot(x, stirlerr(x, use.halves=FALSE) - stM,      type="l", log="x", main="absolute Error")
plot(x, stirlerr(x, use.halves=FALSE) / stM - 1,  type="l", log="x", main="relative Error")
plot(x, abs(stirlerr(x, use.halves=FALSE) / stM - 1), type="l", log="xy",main="|relative Error|")
drawEps.h(-(52:50))
## lgammacor() does *NOT* help, as it is  *designed* for x >= 10 ... (but is interesting there)
##
## ==> Need another chebyshev() or rational-approx. for x in [.1, 7] or so !!

##=============>  see also ../Misc/stirlerr-trms.R  <===============
##                         ~~~~~~~~~~~~~~~~~~~~~~~

cutoffs <- c(15,35,80,500) # cut points, n=*, in the stirlerr() "algorithm"
##
n <- c(seq(1,15, by=1/4),seq(16, 25, by=1/2), 26:30, seq(32,50, by=2), seq(55,1000, by=5),
       20*c(51:99), 50*(40:80), 150*(27:48), 500*(15:20))
st.n <- stirlerr(n, "R3")# rather use.halves=TRUE; but here , use.halves=FALSE
plot(st.n ~ n, log="xy", type="b") ## looks good now (straight line descending {in log-log !}
nM <- mpfr(n, 2048)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose
all.equal(asNumeric(st.nM), st.n)# TRUE
all.equal(st.nM, as(st.n,"mpfr"))# .. difference: 3.381400e-14 was 1.05884.........e-15
all.equal(roundMpfr(st.nM, 64), as(st.n,"mpfr"), tolerance=1e-16)# (ditto)


## --- Look at the direct formula -- why is it not good for n ~= 5 ?
##
## Preliminary Conclusions :
## 1. there is *some* cancellation even for small n (how much?)
## 2. lgamma1p(n) does really not help much compared to lgamma(n+1) --- but a tiny bit in some cases

###  1. Investigating  lgamma1p(n)  vs  lgamma(n+1)  for n < 1 =============================================

##' @title Relative Error of lgamma(n+1) vs lgamma1p() vs MM's stirlerrD2():
##' @param n numeric, typically n << 1
##' @param precBits
##' @return relative error WRT  mpfr(n, precBits)
##' @author Martin Maechler
relE.lgam1 <-  function(n, precBits = if(doExtras) 1024 else 320) {
    M_LN2PI <- 1.837877066409345483 # ~ log(2*Const("pi",60)); very slightly more accurate than log(2*pi)
    st  <- lgamma(n +1)  - (n +0.5)*log(n)  + n  - M_LN2PI/2
    st. <- lgamma1p(n)   - (n +0.5)*log(n)  + n  - M_LN2PI/2      # "lgamma1p"
    st2 <- lgamma(n) + n*(1-(l.n <- log(n))) + (l.n - M_LN2PI)/2  # "MM2"
    st0 <- -(l.n + M_LN2PI)/2                                     # "n0"
    nM <- mpfr(n, precBits)
    stM <- lgamma(nM+1)  - (nM+0.5)*log(nM) + nM - log(2*Const("pi", precBits))/2
    ## stM <- roundMpfr(stM, 128)
    cbind("R3"      = asNumeric(relErrV(stM, st))
        , "lgamma1p"= asNumeric(relErrV(stM, st.))
        , "MM2"	    = asNumeric(relErrV(stM, st2))
        , "n0"      = asNumeric(relErrV(stM, st0))
          )
}

n <- 2^-seq.int(1022, 1, by = -(if(doExtras) 1/4 else 9/4))
relEx <- relE.lgam1(n)
showProc.time()

## Is *equivalent* to 'new' stirlerr_simpl(n version = *)   [not for <mpfr> though, see 'relEmat']:
(simpVer <- eval(formals(stirlerr_simpl)$version))
if(is.null(simpVer)) { warning("got wrong old version of package 'DPQ':")
    print(packageDescription("DPQ"))
    stop("invalid outdated version package 'DPQ'")
}
stir.allS <- function(n) sapply(simpVer, function(v) stirlerr_simpl(n, version=v))
stirS <- stir.allS(n) # matrix
nM <- mpfr(n, 256)  # "high" precision = 256 should suffice!
stirM  <- stirlerr(nM)
releS <- asNumeric(relErrV(stirM, stirS))
          all.equal(relEx, releS, tolerance = 0) # see TRUE on Linux
stopifnot(all.equal(relEx, releS, tolerance = if(noLdbl) 2e-15 else 1e-15))
simpVer3 <- simpVer[simpVer != "lgamma1p"] # have no mpfr-ified lgamma1p()!
## stirlerr_simpl(<mpfr>, *) :
stirM2 <- sapplyMpfr(simpVer3, function(v) stirlerr_simpl(nM, version=v))

## TODO ?:
## apply(stirM2, 2, function(v) all.equal(v, stirS, check.class=FALSE))
## releS2 <- asNumeric(relErrV(stirM2, stirS))
##           all.equal(relEx, releS, tolerance = 0) # see TRUE on Linux
## stopifnot(all.equal(relEx, releS, tolerance = 1e-15))
relEmat <- matrix(NA, ncol(stirM2), ncol(stirS),
                 dimnames = list(simpVer3, colnames(stirS)))
for(j in seq_len(ncol(stirM2)))
    for(k in seq_len(ncol(stirS)))
        relEmat[j,k] <- asNumeric(relErr(stirM2[,j], stirS[,k]))
relEmat
round(-log10(relEmat), 2) # well .. {why? / expected ?}

cols <- c("gray30", adjustcolor(c(2,3,4), 1/2));  lwd <- c(1, 3,3,3)
stopifnot((k <- length(cols)) == ncol(relEx), k == length(lwd))
matplot(n, relEx, type = "l", log="x", col=cols, lwd=lwd, ylim = c(-1,1)*4.5e-16,
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
mtext("really small errors are dominated by small (< 2^-53) errors of log(n)")
## very interesting: there are different intervals  <---> log(n) Qpattern !!
## -- but very small difference, only for n >~= 1/1000  but not before
drawEps.h(negative=TRUE) # abline(h= c(-4,-2:2, 4)*2^-53, lty=c(2,2,2, 1, 2,2,2), col="gray")
legend("topleft", legend = colnames(relEx), col=cols, lwd=3)

## zoomed in a bit:
n. <- 2^-seq.int(400,2, by = -(if(doExtras) 1/4 else 9/4))
relEx. <- relE.lgam1(n.)
matplot(n., relEx., type = "l", log="x", col=cols, lwd=lwd, ylim = c(-1,1)*4.5e-16,
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
drawEps.h(negative=TRUE)
legend("topleft", legend = colnames(relEx.), col=cols, lwd=3)

##====> Absolute errors (and look at  "n0") --------------------------------------
matplot(n., abs19(relEx.), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(4e-17, 5e-16),
        main = quote(abs(relErr(stirlerr_simpl(n, '*')))))
drawEps.h(); legend("top", legend = colnames(relEx.), col=cols, lwd=3)
lines(n., abs19(relEx.[,"n0"]), type = "o", cex=1/4, col=cols[4], lwd=2)

## more zooom-in
n.2 <- 2^-seq.int(85, 50, by= -(if(doExtras) 1/100 else 1/4))
stirS.2 <- sapply(c("R3", "lgamma1p", "n0"), function(v) stirlerr_simpl(n.2, version=v))
releS.2 <- asNumeric(relErrV(stirlerr(mpfr(n.2, 320)), stirS.2))

matplot(n.2, abs19(releS.2), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(4e-17, 5e-16),
        main = quote(abs(relErr(stirlerr_simpl(n, '*')))))
drawEps.h(); legend("top", legend = colnames(releS.2), col=cols, lwd=3)
abline(v = 5e-17, col=(cb <- adjustcolor("skyblue4", 1/2)), lwd=2, lty=3)
axis(1, at=5e-17, col.axis=cb, line=-1/4, cex = 3/4)

matplot(n.2, abs19(releS.2), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(4e-17, 5e-16),
        xaxt="n", xlim = c(8e-18, 1e-15), ## <<<<<<<<<<<<<<<<<<< Zoom-in
        xlab = quote(n), main = quote(abs(relErr(stirlerr_simpl(n, '*')))))
eaxis(1); drawEps.h(); legend("top", legend = colnames(releS.2), col=cols, lwd=3)
abline(v = 5e-17, col=(cb <- adjustcolor("skyblue4", 1/2)), lwd=2, lty=3)
mtext('stirlerr_simpl(*, "n0")  is as good as others for n <= 5e-17', col=adjustcolor(cols[3], 2))
##     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ===> "all *but* "n0" approximations for "larger" small n:
n2 <- 2^-seq.int(20,0, length.out = if(doExtras) 1000 else 200)
relEx2 <- relE.lgam1(n2)[,c("R3", "lgamma1p", "MM2")] # "n0" is "bad": |relE| >= 2.2e-6 !
cols <- c("gray30", adjustcolor(c(2,4), 1/2));  lwd <- c(1, 3,3)
stopifnot((k <- length(cols)) == ncol(relEx2), k == length(lwd))
matplot(n2, relEx2, type = "l", log="x", col=cols, lwd=lwd, ylim = c(-3,3)*1e-15, xaxt="n",
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
eaxis(1, sub10=c(-3,0)); drawEps.h(negative=TRUE)
legend("topleft", legend = colnames(relEx2), col=cols, lwd=3)
##==> "MM" is *NOT* good for n < 1  *but*

## "the same" -- even larger small n:
n3 <- seq(.01, 5, length = if(doExtras) 1000 else 200)
relEx3 <- relE.lgam1(n3)[,c("R3", "lgamma1p", "MM2")] # "no" is "bad" ..
stopifnot((k <- length(cols)) == ncol(relEx3), k == length(lwd))

matplot(n3, relEx3, type = "l", col=cols, lwd=lwd,
        main = "relative errors of direct (approx.) formula for stirlerr(n), small n")
legend("topleft", legend = colnames(relEx3), col=cols, lwd=3)
drawEps.h(negative=TRUE)

matplot(n3, abs19(relEx3), type = "l", col=cols, lwd=lwd,
        log="y", ylim = 2^-c(54, 44), yaxt = "n", ylab = quote(abs(relE)), xlab=quote(n),
        main = "|relative errors| of direct (approx.) formula for stirlerr(n), small n")
eaxis(2, cex.axis=0.9); legend("topleft", legend = colnames(relEx3), col=cols, lwd=3)
drawEps.h()
lines(n3, smooth.spline(abs(relEx3)[,1], df=12)$y, lwd=3, col=cols[1])
lines(n3, smooth.spline(abs(relEx3)[,2], df=12)$y, lwd=3, col=adjustcolor(cols[2], 1/2))
lines(n3, smooth.spline(abs(relEx3)[,3], df=12)$y, lwd=4, col=adjustcolor(cols[3], offset = rep(.2,4)))
## ===>  from  n >~= 1,  "MM2" is definitely better up to n = 5 !!

## Check  log() only :
plot(n, asNumeric(relErrV(log(mpfr(n, 256)), log(n))), ylim = c(-1,1)*2^-53,
     log="x", type="l", xaxt="n") ## ===> indeed --- log(n) approximation pattern !!
eaxis(1) ; drawEps.h(negative=TRUE)
showProc.time()

## =========== "R3"  vs  "lgamma1p" -------------------------- which is better?

## really for the very small n, all is dominated by -(n+0.5)*log(n);  and lgamma1p() is unnecessary!
i <- 1:20; ni <- n[i]
lgamma1p(ni)
- (ni +0.5)*log(ni)  + ni

## much less extreme:
n2 <- lseq(2^-12, 1/2, length = if(doExtras) 1000 else 200)
relE2 <- relE.lgam1(n2)[,-4]

cols <- c("gray30", adjustcolor(2:3, 1/2)); lwd <- c(1,3,3)
matplot(n2, relE2, type = "l", log="x", col=cols, lwd=lwd)
legend("topleft", legend=colnames(relE2), col=cols, lwd=2, lty=1:3)
drawEps.h(negative=TRUE)

matplot(n2, abs19(relE2), type = "l", log="xy", col=cols, lwd=lwd, ylim = c(6e-17, 1e-15),
        xaxt = "n"); eaxis(1, sub10=c(-2,0))
legend("topleft", legend=colnames(relE2), col=cols, lwd=2, lty=1:3)
drawEps.h()
## "MM2" is *worse* here, n < 1/2
for(j in 1:3) lines(n2, smooth.spline(abs(relE2[,j]), df=10)$y, lwd=3,
                    col=adjustcolor(cols[j], 1.5, offset = rep(-1/4, 4)))
## "lgammap  very slightly better in [0.002, 0.05] ...
## "TODO": draw 0.90-quantile curves {--> cobs::cobs() ?} instead of mean-curves?

## which is better?  ... "random difference"
d.absrelE <- abs(relE2[,"R3"]) - abs(relE2[,"lgamma1p"])
plot (n2, d.absrelE, type = "l", log="x", # no clear picture ...
      main = "|relE_R3|  -  |relE_lgamma1p|", axes=FALSE, frame.plot=TRUE)
eaxis(1, sub10=c(-2,1)); eaxis(2); axis(3, at=max(n2)); abline(v = max(n2), lty=3, col="gray")
## 'lgamma1p' very slightly better:
lines(n2, smooth.spline(d.absrelE, df=12)$y, lwd=3, col=2)

## not really small n at all == here see, how "bad" the direct formula gets for 1 < n < 10 or so
n3 <- lseq(2^-14, 2^2, length = if(doExtras) 800 else 150)
relE3 <- relE.lgam1(n3)[, -4]

matplot(n3, relE3, type = "l", log="x", col=cols, lty=1, lwd = c(1,3),
        main = quote(rel.lgam1(n)), xlab=quote(n))

matplot(n3, abs19(relE3),  type = "l", log="xy", col=cols, lwd = c(1,3), xaxt="n",
        main = quote(abs(rel.lgam1(n))), xlab=quote(n), ylim = c(2e-17, 4e-14))
drawEps.h(); eaxis(1, sub10=c(-2,3))
legend("topleft", legend=colnames(relE3), col=cols, lwd=2)
## very small difference --- draw the 3 smoothers :
for(j in 1:3) {
    ll <- lowess(log(n3), abs19(relE3[, j]), f= 1/12)
    with(ll, lines(exp(x), y, col=adjustcolor(cols[j], 1.5), lwd=3))
}
## ==> lgamma1p(.) very slightly in n ~ 10^-4 -- 10^-2 --- but not where it matters: n ~ 0.1 -- 1 !!
##     "MM2"  gets best from  n >~ 1 !
abline(v=1, lty=3, col = adjustcolor(1, 3/4))
showProc.time()


###  2. relErr( stirlerr(.) ) ============================================================

##' Very revealing plot showing the *relative* approximation error of stirlerr(<dblprec>)
##'
p.stirlerrDev <- function(n, precBits = if(doExtras) 2048L else 512L,
                          stnM = stirlerr(mpfr(n, precBits), use.halves=use.halves, verbose=verbose),
                          abs = FALSE,
                          ## cut points, n=*, in the stirlerr() algorithm; "FIXME": sync with ../R/dgamma.R <<<<
                          scheme = c("R3", "R4.4_0"),
                          cutoffs = switch(match.arg(scheme)
                                         , R3   = c(15, 35, 80, 500)
                                         , R4.4_0 = c(4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.7,
                                                      6.1, 6.5, 7,  7.9, 8.75, 10.5, 13,
                                                      20, 26, 60, 200, 3300, 17.4e6)
                                           ## {FIXME: need to sync} <==> ../man/stirlerr.Rd <==> ../R/dgamma.R
                                           ),
                          use.halves = missing(cutoffs),
                          direct.ver = c("R3", "lgamma1p", "MM2", "n0"),
                          verbose = getOption("verbose"),
                          type = "b", cex = 1,
                          col = adjustcolor(1, 3/4), colnB = adjustcolor("orange4", 1/3),
                          log = if(abs) "xy" else "x",
                          xlim=NULL, ylim = if(abs) c(8e-18, max(abs(N(relE)))))
{
    op <- par(las = 1, mgp=c(2, 0.6, 0))
    on.exit(par(op))
    require("Rmpfr"); require("sfsmisc")
    st <- stirlerr(n, scheme=scheme, cutoffs=cutoffs, use.halves=use.halves, direct.ver=direct.ver,
                   verbose=verbose)
    relE <- relErrV(stnM, st) # eps0 = .Machine$double.xmin
    N <- asNumeric
    form <- if(abs) abs(N(relE)) ~ n else N(relE) ~ n
    plot(form, log=log, type=type, cex=cex, col=col, xlim=xlim, ylim=ylim,
         ylab = quote(relErrV(stM, st)), axes=FALSE, frame.plot=TRUE,
         main = sprintf("stirlerr(n, cutoffs) rel.error [wrt stirlerr(Rmpfr::mpfr(n, %d))]",
                        precBits))
    eaxis(1, sub10=3)
    eaxis(2)
    mtext(paste("cutoffs =", deparse1(cutoffs)))
    ylog <- par("ylog")
    ## FIXME:  improve this ---> drawEps.h() above
    if(ylog) {
        epsC <- c(1,2,4,8)*2^-52
        epsCxp <- expression(epsilon[C],2*epsilon[C], 4*epsilon[C], 8*epsilon[C])
    } else {
        epsC <- (-2:2)*2^-52
        epsCxp <- expression(-2*epsilon[C],-epsilon[C], 0, +epsilon[C], +2*epsilon[C])
    }
    dy <- diff(par("usr")[3:4])
    if(diff(range(if(ylog) log10(epsC) else epsC)) > dy/50) {
        lw <- rep(1/2, 5); lw[if(ylog) 1 else 3] <- 2
        abline(  h=epsC, lty=3, lwd=lw)
        axis(4, at=epsC, epsCxp, las=2, cex.axis = 3/4, mgp=c(3/4, 1/4, 0), tck=0)
    } else ## only x-axis
        abline(h=if(ylog) epsC else 0, lty=3, lwd=2)
    abline(v = cutoffs, col=colnB)
    axis(3, at=cutoffs, col=colnB, col.axis=colnB,
         labels = formatC(cutoffs, digits=3, width=1))
    invisible(relE)
} ## p.stirlerrDev()

showProc.time()

n <- lseq(2^-10, 1e10, length = if(doExtras) 4096 else 512)
n <- lseq(2^-10, 5000, length = if(doExtras) 4096 else 512)
## store "expensive" stirlerr() result, and re-use many times below:
nM <- mpfr(n, if(doExtras) 2048 else 512)
st.nM <- stirlerr(nM, use.halves=FALSE) ## << on purpose

p.stirlerrDev(n=n, stnM=st.nM, use.halves = FALSE) # default cutoffs= c(15, 40, 85, 600)
p.stirlerrDev(n=n, stnM=st.nM, use.halves = FALSE, ylim = c(-1,1)*1e-12) # default cutoffs= c(15, 40, 85, 600)

## show the zoom-in region in next plot
yl2 <- 3e-14*c(-1,1)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

if(do.pdf) { dev.off() ; pdf("stirlerr-relErr_1.pdf") }

## drop n < 7:
p.stirlerrDev(n=n, stnM=st.nM, xlim = c(7, max(n)), use.halves=FALSE) # default cutoffs= c(15, 40, 85, 600)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

## The first plot clearly shows we should do better:
## Current code is switching to less terms too early, loosing up to 2 decimals precision
if(FALSE) # no visible difference {use.halves = T / F }:
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2, use.halves = FALSE)
p.stirlerrDev(n=n, stnM=st.nM, ylim = yl2, use.halves = TRUE)# exact at n/2 (n <= ..)
abline(h = yl2, col=adjustcolor("tomato", 1/4), lwd=3, lty=2)

showProc.time()


if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-1.pdf") }

### ~19.April 2021: "This is close to *the* solution" (but see 'cuts' below)
cuts <- c(7, 12, 20, 26, 60, 200, 3300)
##        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st.   <- stirlerr(n=n , cutoffs = cuts, verbose=TRUE)
st.nM <- stirlerr(n=nM, cutoffs = cuts, use.halves=FALSE) ## << on purpose
relE <- asNumeric(relErrV(st.nM, st.))
head(cbind(n, relE), 20)
## nice printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(relE, 4))
      , quote=FALSE)

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts)
## and zoom in:
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim = yl2)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, ylim = yl2/20)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-2.pdf") }

## zoom in ==> {good for n >= 10}
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", ylim = 2e-15*c(-1,1),
              cutoffs = cuts)## old default cutoffs = c(15,35, 80, 500)

if(do.pdf) { dev.off(); pdf("stirlerr-relErr_6-fin-3.pdf") }
showProc.time()


##-- April 20: have more terms up to S10 in stirlerr() --> can use more cutoffs
n <- n5m <- lseq(1/64, 5000, length = if(doExtras) 4096 else 512)
nM <- mpfr(n, if(doExtras) 2048L # a *lot* accuracy for stirlerr(nM,*)
	      else 512L)
ct10.1 <- c(            5.4, 7.5, 8.5, 10.625, 12.125, 20, 26, 60, 200, 3300)# till 2024-01-19
ct10.2 <- c(            5.4, 7.9, 8.75,10.5  , 13,     20, 26, 60, 200, 3300)
cuts <-
ct12.1 <- c(5.22, 6.5,  7.0, 7.9, 8.75,10.5  , 13,     20, 26, 60, 200, 3300)
##        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5.25 is "too small" but the direct formula is already really bad there, ...
st.nM <- roundMpfr(stirlerr(nM, use.halves=FALSE, ## << on purpose;
                            verbose=TRUE), precBits = 128)
## NB: for x=xM <mpfr>; `cutoffs` are *not* used.
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0")
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0", ylim = c(-1,1)*3e-15)
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0", ylim = c(-1,1)*1e-15)

p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", scheme = "R4.4_0", abs=TRUE)
axis(1,at= 2:6, col=NA, col.axis=(cola <- "lightblue"), line=-3/4)
abline(v = 2:6, lty=3, col=cola)
if(FALSE)## using exact values sferr_halves[] *instead* of MPFR ones: ==> confirmation they lay on top
lines((0:30)/2, abs(stirlerr((0:30)/2, cutoffs=cuts, verbose=TRUE)/DPQ:::sferr_halves - 1), type="o", col=2,lwd=2)

if(FALSE) ## nice (but unneeded) printout :
print(cbind(n       = format(n, drop0trailing = TRUE),
            stirlerr= format(st.,scientific=FALSE, digits=4),
            relErr  = signif(relE, 4))
      , quote=FALSE)

showProc.time()


## ========== where should the cutoffs be ? ===================================================

.stirl.cutoffs <- function(scheme)
    eval(do.call(substitute, list(formals(stirlerr)$cutoffs, list(scheme = scheme))))
drawCuts <- function(scheme, axis=NA, lty = 3, col = "skyblue", ...) {
    abline(v = (ct <- .stirl.cutoffs(scheme)), lty=lty, col=col, ...)
    if(is.finite(axis)) axisCuts(side = axis, at = ct, col=col, ...)
}
axisCuts <- function(scheme, side = 3, at = .stirl.cutoffs(scheme), col = "skyblue", line = -3/4, ...)
    axis(side, at=at, labels=formatC(at), col.axis = col, col=NA, col.ticks=NA, line=line, ...)
mtextCuts <- function(cutoffs, scheme, ...) {
    if(!missing(scheme)) cutoffs <- .stirl.cutoffs(scheme)
    mtext(paste("cutoffs =", deparse1(cutoffs)), ...)
}


if(do.pdf) { dev.off(); pdf("stirlerr-tst_order_k.pdf") }

mK <- 20L # := max(k)
## order = k = 1:mK  terms in series approx:
k <- 1:mK
n <- 2^seq(1, 28, by=1/16)
nM <- mpfr(n, 1024)
stnM <- stirlerr(nM) # the "true" values

stirlOrd <- sapply(k, function(k.) stirlerr(n, order = k.))
stirlO_lgcor <- cbind(stirlOrd, sapply(5:6, function(nal) lgammacor(n, nalgm = nal)))
relE <- asNumeric(stirlO_lgcor/stnM -1) # "true" relative error

## use a "smooth" but well visible polette :
palROBG <- colorRampPalette(c("red", "darkorange2", "blue", "seagreen"), space = "Lab")
palette(adjustcolor(palROBG(mK+2), 3/4))
##                            -- 2 lgamcor()'s

(tit.k  <- substitute(list(    stirlerr(n, order=k) ~~"error",  k == 1:mK),  list(mK = mK)))
(tit.kA <- substitute(list(abs(stirlerr(n, order=k) ~~"error"), k == 1:mK),  list(mK = mK)))
lgammacorTit <- function(...) mtext("+ lgammacor(x, 5) [bad] + lgammacor(x, 6) [good]", col=1, ...)

matplotB(n, relE, cex=2/3, ylim = c(-1,1)*1e-13, col=k,
         log = "x", xaxt="n", main = tit.k)
lgammacorTit()
eaxis(1, nintLog = 20)
drawCuts("R4.4_0")

## zoom in (ylim)
matplotB(n, relE, cex=2/3, ylim = c(-1,1)*5e-15, col=k,
        log = "x", xaxt="n", main = tit.k)
lgammacorTit()
eaxis(1, nintLog = 20); abline(h = (-2:2)*2^-53, lty=3, lwd=1/2)
drawCuts("R4.4_0", axis = 3)

## log-log  |rel.Err|  -- "linear"
matplotB(n, abs19(relE), cex=2/3, col=k, ylim = c(8e-17, 1e-7), log = "xy", main=tit.kA)
mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
lgammacorTit(line=-1)
drawCuts("R4.4_0", axis = 3)

## zoom in -- still "large" {no longer carry the  lgammacor() .. }:
n2c <- 2^seq(2, 8, by = 1/(if(doExtras) 256 else 32))
nMc <- mpfr(n2c, 1024)
stnMc <- stirlerr(nMc) # the "true" values
stirlOrc <- sapply(k, function(k.) stirlerr(n2c, order = k.))
relEc <- asNumeric(stirlOrc/stnMc -1) # "true" relative error

matplotB(n2c, relEc, cex=2/3, ylim = c(-1,1)*1e-13, col=k,
         log = "x", xaxt="n", main = tit.k)
eaxis(1, sub10 = 2)
drawCuts("R4.4_0", axis=3)

## log-log  |rel.Err|  -- "linear"
matplotB(n2c, abs19(relEc), cex=2/3, col=k, ylim = c(8e-17, 1e-3), log = "xy", main=tit.kA)
mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
drawCuts("R4.4_0", axis = 3)


## zoom into the critical n region
nc <- seq(3.5, 11, by = 1/(if(doExtras) 128 else 16))
ncM <- mpfr(nc, 256)
stncM <- stirlerr(ncM) # the "true" values
stirlO.c <- sapply(k, function(k) stirlerr(nc, order = k))
relEc <- asNumeric(stirlO.c/stncM -1) # "true" relative error


## log-log  |rel.Err|  -- "linear"
matplotB(nc, abs19(relEc), cex=2/3, col=k, ylim = c(2e-17, 1e-8),
        log = "xy", xlab = quote(n), main = quote(abs(relErr(stirlerr(n, order==k)))))
mtext(paste("k =", deparse(k))) ; drawEps.h(lwd = 1/2)
lines(nc, abs19(asNumeric(stirlerr_simpl(nc, "R3" )/stncM - 1)), lwd=1.5, col=adjustcolor("thistle", .6))
lines(nc, abs19(asNumeric(stirlerr_simpl(nc, "MM2")/stncM - 1)), lwd=4,   col=adjustcolor(20, .4))
## lines(nc, abs19(asNumeric(stirlerr_simpl(nc,"MM2")/stncM - 1)), lwd=3, col=adjustcolor("purple", 2/3))
legend(10^par("usr")[1], 1e-9, legend=paste0("k=", k), bty="n", lwd=2,
       col=k, lty=1:5, pch= c(1L:9L, 0L, letters)[seq_along(k)])
drawCuts("R4.4_0", axis=3)

## Zoom-in [only]
matplotB(nc, abs19(relEc), cex=2/3, col=k, ylim = c(4e-17, 1e-11), xlim = c(4.8, 6.5),
        log = "xy", xlab = quote(n), main = quote(abs(relErr(stirlerr(n, order==k)))))
mtext(paste("k =", deparse(k))) ; drawEps.h(lwd = 1/2)
lines(nc, abs19(asNumeric(stirlerr_simpl(nc, "R3" )/stncM - 1)), lwd=1.5, col=adjustcolor("thistle", .6))
lines(nc, abs19(asNumeric(stirlerr_simpl(nc, "MM2")/stncM - 1)), lwd=4,   col=adjustcolor(20, .4))

k. <- k[-(1:6)]
legend("bottomleft", legend=paste0("k=", k.), bty="n", lwd=2,
       col=k., lty=1:5, pch= c(1L:9L, 0L, letters)[k.])
drawCuts("R4.4_0", axis=3)

showProc.time()

##--- Accuracy of "R4.4_0" -------------------------------------------------------

for(nc in list(seq(4.75, 28, by=1/(if(doExtras) 512 else 64)), # for a bigger pix
               seq(4.75,  9, by=1/(if(doExtras) 1024 else 128))))
{
  ncM <- mpfr(nc, 1024)
  stncM <- stirlerr(ncM) # the "true" values
  stirl.440 <- stirlerr(nc, scheme = "R4.4_0")
  stirl.3   <- stirlerr(nc, scheme = "R3")
  relE440 <- asNumeric(relErrV(stncM, stirl.440))
  relE3   <- asNumeric(relErrV(stncM, stirl.3  ))
  ##
  plot(nc, abs19(relE440), xlab=quote(n), main = quote(abs(relErr(stirlerr(n, '"R4.4_0"')))),
       type = "l", log = "xy", ylim = c(4e-17, 1e-13))
  mtextCuts(scheme="R4.4_0", cex=4/5)
  drawCuts("R4.4_0", lty=2, lwd=2, axis=4)
  drawEps.h()
  if(max(nc) <= 10) abline(v = 5+(0:20)/10, lty=3, col=adjustcolor(4, 1/2))
  if(TRUE) { # but just so ...
      c3 <- adjustcolor("royalblue", 1/2)
      lines(nc, pmax(abs(relE3), 1e-18), col=c3)
      title(quote(abs(relErr(stirlerr(n, '"R3"')))), adj=1, col.main = c3)
      drawCuts("R3", lty=4, col=c3); mtextCuts(scheme="R3", adj=1, col=c3)
  }
  addOrd <- TRUE
  addOrd <- dev.interactive(orNone=TRUE)
  if(addOrd) {
      if(max(nc) >= 100) {
          i <- (15 <= nc & nc <= 85) # (no-op in [4.75, 7] !)
          ni <- nc[i]
      } else { i <- TRUE; ni <- nc }
      for(k in 8:17) lines(ni, abs19(asNumeric(relErrV(stncM[i], stirlerr(ni, order=k)))), col=adjustcolor(k, 1/3))
      title(sub = "stirlerr(*, order k = 8:17)")
 }
} ## for(nc ..)

if(FALSE)
lines(nc, abs19(relE440))# *re* draw!

showProc.time()


palette("Tableau")

## Focus more:
stirlerrPlot <- function(nc, k, res=NULL, legend.xy = "left", full=TRUE, precB = 1024,
                         ylim = c(3e-17, 2e-13), cex = 5/4) {
    stopifnot(require("Rmpfr"), require("graphics"))
    if(is.list(res) && all(c("nc", "k", "relEc","splarelE") %in% names(res))) { ## do *not* recompute
        list2env(res, envir = environment())
    } else { ## compute
        stopifnot(is.finite(nc), nc > 0, length(nc) >= 100, k == as.integer(k), 0 <= k, k <= 20)
        ncM <- mpfr(nc, precB)
        stncM <- stirlerr(ncM) # the "true" values
        stirlO.c <- sapply(k, function(k) stirlerr(nc, order = k))
        relEc <- asNumeric(stirlO.c/stncM -1) # "true" relative error
        ## log  |rel.Err|  -- "linear"
        ## smooth() on log-scale {and transform back}:
        splarelE <- apply(log(abs19(relEc)), 2, function(y) exp(smooth.spline(y, df=4)$y))
        ## the direct formulas (default "R3", "MM2"):
        arelEs0 <- abs(asNumeric(stirlerr_simpl(nc       )/stncM - 1))
        arelEs2 <- abs(asNumeric(stirlerr_simpl(nc, "MM2")/stncM - 1))
    }
    pch. <- c(1L:9L, 0L, letters)[k]
    if(full)
        matplotB(nc, abs19(relEc), col=k, pch = pch., cex=cex, ylim=ylim, log = "y",
                 xlab = quote(n), main = quote(abs(relErr(stirlerr(n, order==k)))))
    else ## smooth only
        matplotB(nc, splarelE, col=adjustcolor(k,2/3), pch=pch., lwd=2, cex=cex, ylim=ylim, log = "y",
                 xlab = quote(n), main = quote(abs(relErr(stirlerr(n, order==k)))))
    mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
    legend(legend.xy, legend=paste0("k=", k), bty="n", lwd=2, col=k, lty=1:5, pch = pch.)
    abline(v = 5+(0:20)/10, lty=3, col=adjustcolor(10, 1/2))
    drawCuts("R4.4_0", axis=3)
    if(full) {
        matlines(nc, splarelE, col=adjustcolor(k,2/3), lwd=4)
        lines(nc, pmax(arelEs0, 1e-19), lwd=1.5, col=adjustcolor( 2, 0.2))
        lines(nc, pmax(arelEs2, 1e-19), lwd=1,   col=adjustcolor(10, 0.2))
    }
    lines(nc, smooth.spline(arelEs0, df=12)$y, lwd=3, col= adjustcolor( 2, 1/2))
    lines(nc, smooth.spline(arelEs2, df=12)$y, lwd=3, col= adjustcolor(10, 1/2))
    invisible(list(nc=nc, k=k, relEc = relEc, splarelE = splarelE, arelEs0=arelEs0, arelEs2=arelEs2))
}

rr1 <- stirlerrPlot(nc = seq(4.75, 9.0, by=1/(if(doExtras) 1024 else 256)),
                    k = 7:20)
stirlerrPlot(res = rr1, full=FALSE, ylim = c(8e-17, 1e-13))
if(interactive())
stirlerrPlot(res = rr1)

rr <- stirlerrPlot(nc = seq(5, 6.25, by=1/(if(doExtras) 2048 else 256)), k = 9:18)
stirlerrPlot(res = rr, full=FALSE, ylim = c(8e-17, 1e-13))

showProc.time()


palette("default")

if(do.pdf) { dev.off(); pdf("stirlerr-tst_order_k-vs-k1.pdf") }

##' Find 'cuts', i.e.,  a region  c(k) +/- s(k) i.e. intervals [c(k) - s(k), c(k) + s(k)]
##' where c(k) is such that relE(n=c(k), k) ~= eps)

##' 1. Find the c1(k)  such that |relE(n, k)| ~= c1(k) * n^{-2k}
findC1 <- function(n, ks, e1 = 1e-15, e2 = 1e-5, res=NULL,
                   precBits = if(doExtras) 1024 else 256, do.plot = TRUE, ...)
{
    if(is.list(res) && all(c("n", "ks", "arelE") %in% names(res))) { ## do *not* recompute, take from 'res':
        list2env(res, envir = environment())
    } else { ## compute
        stopifnot(require("Rmpfr"),
                  is.numeric(ks), ks == (k. <- as.integer(ks)), length(ks <- k.) >= 1,
                  length(e1) == 1L, length(e2) == 1L, is.finite(c(e1,e2)), e1 >= 0, e2 >= e1,
                  0 <= ks, ks <= 20, is.numeric(n), n > 0, is.finite(n), length(n) >= 100)
        nM <- mpfr(n, precBits)
        stirM <- stirlerr(nM) # the "true" values
        stirlOrd <- sapply(ks, function(k) stirlerr(n, order = k))
        arelE <- abs(asNumeric(stirlOrd/stirM -1)) # "true" relative error
    }
    arelE19 <- pmax(arelE, 1e-19)
    ## log  |rel.Err|  -- "linear"
    ## on log-scale {and transform back; for linear fit, only use values inside [e1, e2]
    if(do.plot) # experi
        matplotB(n, arelE19, log="xy", ...)
    ##  matplot(n, arelE19, type="l", log="xy", xlim = c(min(n), 20))

    ## re-compute these, as they *also* depend  on (e1, e2)
    c1 <- vapply(seq_along(ks), function(i) {
        k <- ks[i]
        y <- arelE19[,i]
        iUse <- e1 <= y & y <= e2
        if(sum(iUse) < 10) stop("only", sum(iUse), "values in [e1,e2]")
            ## .lm.fit(cbind(1, log(n[iUse])), log(y[iUse]))$coefficients
        ## rather, we *know* the error is  c* n^{-2k} , i.e.,
        ##      log |relE| = log(c) - 2k * log(n)
        ## <==> c = exp( log|relE| + 2k * log(n))
        exp(mean(log(y[iUse]) + 2*k * log(n[iUse])))
    }, numeric(1))
    if(do.plot) {
        drawEps.h()
        for(i in seq_along(ks))
            lines(n, c1[i] * n^(-2*ks[i]), col=adjustcolor(i, 1/3), lwd = 4, lty = 2)
    }
    invisible(list(n=n, ks=ks, arelE = arelE, c1 = c1))
} ## findC1()

c1.Res <- findC1(n = 2^seq(2, 26, by = 1/(if(doExtras) 128 else 32)), ks = 1:18)
(s.c1.fil <- paste0("stirlerr-c1Res-", myPlatform(), ".rds"))
saveRDS(c1.Res, file = s.c1.fil)


if(!exists("c1.Res")) {
    c1.Res <- readRDS(s.c1.fil)
    ## re-do the "default" plot of findC1():
    findC1(res = c1.Res, xaxt="n"); eaxis(1, sub10=2)
}

## the same, zoomed in:
findC1(res = c1.Res, xlim = c(4, 40), ylim = c(2e-17, 1e-12))
ks <- c1.Res$ks; pch. <- c(1L:9L, 0L, letters)[ks]
legend("left", legend=paste0("k=", ks), bty="n", lwd=2, col=ks, lty=1:5, pch = pch.)

## smaller set : larger e1 :
c1.r2 <- findC1(res = c1.Res, xlim = c(4, 30), ylim = c(4e-17, 1e-13), e1 = 4e-15)
legend("left", legend=paste0("k=", ks), bty="n", lwd=2, col=ks, lty=1:5, pch = pch.)

print(digits = 4,
      cbind(ks, c1. = c1.Res$c1, c1.2 = c1.r2$c1,
            relD = round(relErrV(c1.r2$c1, c1.Res$c1), 4)))

c1.. <- c1.r2[c("ks", "c1")] # just the smallest part is needed here:
(s.c1.fil <- paste0("stirlerr-c1r2-", myPlatform(), ".rds"))
saveRDS(c1.., file = s.c1.fil) # was "stirlerr-c1.rds"

## 2. Now, find the n(k) +/- se.n(k)  intervals
## Use these c1 from above

##' Given c1-results relErr(n,k);  |relE(n,k)  ~= c1 * n^{-2k} , find  n such that
##'  |relE(n,k)| line {in log-log scale} cuts  y = eps, i.e., n* such that |relE(n,k)| <= eps for all n >= n*
n.ep <- function(eps, c1Res, ks = c1Res$ks, c1 = c1Res$c1, ...) {
    stopifnot(is.finite(eps), length(eps) == 1L, eps > 0,
              length(ks) == length(c1), is.numeric(c1), is.integer(ks), ks >= 1)
    ## n: given k, the location where the  |relE(n,k)| line {in log-log} cuts  y = eps
    ##     |relE(n,k)  ~= c1 * n^{-2k}          <==>
    ##  log|relE(n,k)| ~= log(c1) - 2k* log(n)  <==>
    ## c := mean{  exp( log|relE(n,k)| + 2k* log(n) ) }  ------- see findC1()
    ## now,  solve for  n :
    ##     c1 * n^{-2k}      == eps
    ## log(c1) - 2k* log(n)  == log(eps)
    ##      log(n)  ==     (log(eps) - log(c1)) / (-2k)	<==>
    ##        n     ==  exp((log(c1) - log(eps)) / 2k)
    exp((log(c1) - log(eps))/(2*ks))
}

## get c1..
if(!exists("c1..")) c1.. <- readRDS("stirlerr-c1.rds")

ne2 <- n.ep(2^-51, c1Res = c1..) ## ok
ne1 <- n.ep(2^-52, c1Res = c1..)
ne. <- n.ep(2^-53, c1Res = c1..)

form <- function(n) format(signif(n, 3), scientific=FALSE)
data.frame(k = ks, ne2 = form(ne2), ne1 = form(ne1), ne. = form(ne.),
           cutoffs = form(rev(.stirl.cutoffs("R4.4_0")[-1])))

##    ------- Linux F 36/38 x86_64 (nb-mm5|v-lynne)
##  k        ne2         ne1         ne.     cutoffs
##  1 8660000.00 12200000.00 17300000.00 17400000.00
##  2    2150.00     2560.00     3040.00     3700.00
##  3     159.00      178.00      200.00      200.00
##  4      46.60       50.80       55.40       81.00
##  5      23.40       25.10       26.90       36.00
##  6      15.20       16.10       17.10       25.00
##  7      11.50       12.10       12.70       19.00
##  8       9.43        9.85       10.30       14.00
##  9       8.20        8.53        8.86       11.00
## 10       7.42        7.68        7.95        9.50
## 11       6.89        7.11        7.34        8.80
## 12       6.53        6.72        6.92        8.25
## 13       6.27        6.44        6.62        7.60
## 14       6.10        6.25        6.41        7.10
## 15       5.98        6.12        6.26        6.50 * (not used)
## 16       5.90        6.03        6.16        6.50 *   "   "
## 17       5.85        5.97        6.10        6.50 *   "   "
## 18       5.83        5.95        6.06        6.50 << used all the way down to 5.25

## ok --- correct order of magnitude ! --- good!


## 2b.  find *interval* around the  'n(eps)' values

## -- Try simply
d.k <- ne. - ne1
## interval
int.k <- cbind(ne1 - d.k,
               ne1 + d.k)
## look at e.g.
data.frame(k=ks, `n(k)` = form(ne1), int = form(int.k))
##   k        n.k.       int.1       int.2
##   1 12200000.00  7180000.00 17300000.00
##   2     2560.00     2070.00     3040.00
##   3      178.00      156.00      200.00
##   4       50.80       46.20       55.40
##   5       25.10       23.30       26.90
##   6       16.10       15.20       17.10
##   7       12.10       11.40       12.70
##   8        9.85        9.41       10.30
##   9        8.53        8.19        8.86
##  10        7.68        7.41        7.95
##  11        7.11        6.88        7.34
##  12        6.72        6.52        6.92
##  13        6.44        6.27        6.62
##  14        6.25        6.10        6.41
##  15        6.12        5.98        6.26
##  16        6.03        5.90        6.16
##  17        5.97        5.85        6.10
##  18        5.95        5.83        6.06

##' as function {well, *not* computing c1.k from scratch
nInt <- function(k, c1.k, ep12 = 2^-(52:53)) {
    if(length(k) == 1L) { # special convention to call for *one* k, with c1.k vector
        stopifnot(k == (k <- as.integer(k)), k >= 1, length(c1.k) >= k)
        c1.k <- c1.k[k]
    }

    ## see n.ep() above
    n_ <- function(eps, k, c1) {
        stopifnot(is.finite(eps), length(eps) == 1L, eps > 0,
                  length(k) == length(c1), is.numeric(c1), is.integer(k), k >= 1)
        exp((log(c1) - log(eps))/(2*k))
    }

    ne1 <- n_(ep12[1], k, c1.k)
    ne. <- n_(ep12[2], k, c1.k)
    d.k <- ne. - ne1
    stopifnot(d.k > 0)
    ## interval: {"fudge" 0.5 / 2.5} from results -- also 'noLdbl' gives *quite* different pic!:
    odd <- k %% 2 == 1
    cbind(ne1 - ifelse( odd & noLdbl, 1.5, 0.5) * d.k,
          ne1 + ifelse(!odd         , 8, 2.5) * d.k)
}

nInt(k= 1, c1..$c1)
nInt(k= 2, c1..$c1)
nInt(k=18, c1..$c1)

nints.k <- nInt(ks, c1..$c1)
## for printing
form(as.data.frame( nints.k ))

## (-.5 , +2.5)             ## originally ( -1, +1)
## 1  9710000.00 24900000.00    # 7180000.00 17300000.00
## 2     2320.00     3770.00    #    2070.00     3040.00
## 3      167.00      232.00    #     156.00      200.00
## 4       48.50       62.30    #      46.20       55.40
## 5       24.20       29.60    #      23.30       26.90
## 6       15.70       18.50    #      15.20       17.10
## 7       11.70       13.60    #      11.40       12.70
## 8        9.63       10.90    #       9.41       10.30
## 9        8.36        9.36    #       8.19        8.86
## 10       7.54        8.36    #       7.41        7.95
## 11       7.00        7.68    #       6.88        7.34
## 12       6.62        7.21    #       6.52        6.92
## 13       6.36        6.88    #       6.27        6.62
## 14       6.17        6.64    #       6.10        6.41
## 15       6.05        6.48    #       5.98        6.26
## 16       5.96        6.36    #       5.90        6.16
## 17       5.91        6.28    #       5.85        6.10
## 18       5.89        6.24    #       5.83        6.06

## 3. Then for each of the intervals, compare  order  k  vs k+1
## --  ==> optimal cutoff  { how much platform dependency ?? }

### Here, compute only
find1cuts <- function(k, c1,
                      n = nInt(k, c1), # the *set* of n's  or the 'range'
                      len.n = if(doExtras) 1000 else 250,
                      precBits =  if(doExtras) 1024 else 256, nM = mpfr(n, precBits),
                      stnM = stirlerr(nM),
                      stirlOrd = sapply(k+(0:1), function(.k.) stirlerr(n, order = .k.)),
                      relE = asNumeric(stirlOrd/stnM -1), # "true" relative error for the {k, k+1}
                      do.spl=TRUE, df.spline = 9, # df = 5 gives 2 cutpoints for k==1
                      do.low=TRUE, f.lowess = 0.2,
                      do.cobs = require("cobs", quietly=TRUE), tau = 0.90)
{
    ## check  relErrV( stirlerr(n, order=k )  vs
    ##                 stirlerr(n, order=k+1)
    if(length(n) == 2L)
        if(n[1] < n[2]) n <- seq(n[1], n[2], length.out = len.n)
        else stop("'n' must be *increasing")
    force(relE)
    y <- abs19(relE)
    ## NB: all smoothing --- as in stirlerrPlot() above -- should happen in log-space
    ## (log(abs19(relEc)), 2, function(y) exp(smooth.spline(y, df=4)$y))
    ly <- log(y) # == log(abs19(relE)) == log(max(|r|, 1e-19))
    if(do.spl) {##
        s1 <- exp(smooth.spline(ly[,1], df=df.spline)$y)
        s2 <- exp(smooth.spline(ly[,2], df=df.spline)$y)
    }
    if(do.low) { ## lowess
        s1l <- exp(lowess(ly[,1], f=f.lowess)$y)
        s2l <- exp(lowess(ly[,2], f=f.lowess)$y)
    }
    ## also use cobs() splines for the 90% quantile !!
    EE <- environment() # so can  set do.cobs to FALSE in case of error
    if(do.cobs) { ## <==> require("cobs") # yes, this is in tests/
        cobsF <- function(Y) cobs(n, Y, tau=tau, nknots = 6, lambda = -1,
                                  print.warn=FALSE, print.mesg=FALSE)
        ## sparseM::chol(<matrix.csr>) now gives error when it gave warning {about singularity}
        cobsF <- function(Y) {
            r <- tryCatch(cobs(n, Y, tau=tau, nknots = 6, lambda = -1,
                               print.warn=FALSE, print.mesg=FALSE),
                          error = identity)
            if(inherits(r, "error")) {
                assign("do.cobs", FALSE, envir = EE) # and return
                list(fitted = FALSE)
            }
            else
                r
        }
        cs1 <- exp(cobsF(ly[,1])$fitted)
        cs2 <- exp(cobsF(ly[,2])$fitted)
    }
    smooths <- list(spl = if(do.spl ) cbind(s1, s2 ),
                    low = if(do.low ) cbind(s1l,s2l),
                    cobs= if(do.cobs) cbind(cs1,cs2))
    ## diffL <- list(spl = if(do.spl ) s2 -s1,
    ##                   low = if(do.low ) s2l-s1l,
    ##                   cobs= if(do.cobs) cs2-cs1)
### FIXME: simplification does not always *work* -- (i, n.) are not always ok
### ------ notably within R-devel-no-ldouble  {hence probably macOS M1 ... ..}
    sapply(smooths, function(s12) { d <- s12[,2] - s12[,1]
        ## typically a (almost or completely) montone increasing function, crossing zero *once*
        ## compute cutpoint:
        ## i := the first n[i] with d(n[i]) >= 0
        i <- which(d >= 0)[1]
        if(length(i) == 1L && !is.na(i) && i > 1L) {
            i_ <- i - 1L  # ==> d(n[i_]) < 0
            ## cutpoint must be in [n[i_], n[i]] --- do linear interpolation
            n. <- n[i_] - (n[i] - n[i_])* d[i_] / (d[i] - d[i_])
        } else {
            if(length(i) != 1L) i <- -length(i)
            n. <- NA_integer_
        }
        c(i=i, n.=n.)
    }) -> n.L

    list(k=k, n=n, relE = unname(relE), smooths=smooths, i.n = n.L)
} ## find1cuts()

k. <- 1:15
system.time(
    ## Failed on lynne [2024-06-04, R 4.4.1 beta]  with
    ## Error in .local(x, ...) : insufficient space  ---> SparseM :: chol(<..>)
    ## ==> now we catch this inside find1cuts():
    resL  <- lapply(setNames(,k.), function(k) find1cuts(k=k, c1=c1..$c1))
) ## -- warnings, notably from cobs() not converging
## needs  12 sec (!!)  user  system elapsed =
(s.find15.fil <- paste0("stirlerr-find1_1-15_", myPlatform(), ".rds"))
## now we catch cobs() errors {from SparseM::chol},
## okCuts <- !inherits(resL, "error")
## if(okCuts) {
saveRDS(resL, file = s.find15.fil) # was  "stirlerr-find1_1-15.rds"
## } else  traceback()
## 11: stop(mess)
## 10: .local(x, ...)
## 9: chol(e, tmpmax = tmpmax, nsubmax = nsubmax, nnzlmax = nnzlmax)
## 8: chol(e, tmpmax = tmpmax, nsubmax = nsubmax, nnzlmax = nnzlmax)
## 7: rq.fit.sfnc(Xeq, Yeq, Xieq, Yieq, tau = tau, rhs = rhs, control = rqCtrl)
## 6: drqssbc2(x, y, w, pw = pw, knots = knots, degree = degree, Tlambda = if (select.lambda) lambdaSet else lambda,
##        constraint = constraint, ptConstr = ptConstr, maxiter = maxiter,
##        trace = trace - 1, nrq, nl1, neqc, niqc, nvar, tau = tau,
##        select.lambda = select.lambda, give.pseudo.x = keep.x.ps,
##        rq.tol = rq.tol, tol.0res = tol.0res, print.warn = print.warn)
## 5: cobs(n, Y, tau = tau, nknots = 6, lambda = -1, print.warn = FALSE,
##        print.mesg = FALSE) at stirlerr-tst.R!udBzpT#32
## 4: cobsF(ly[, 1]) at stirlerr-tst.R!udBzpT#34
## 3: find1cuts(k = k, c1 = c1..$c1) at #1
## 2: FUN(X[[i]], ...)
## 1: lapply(setNames(, k.), function(k) find1cuts(k = k, c1 = c1..$c1))

ok1cutsLst <- function(res) {
    stopifnot(is.list(res), sapply(res, is.list)) # must be list of lists
    cobsL <- lapply(lapply(res, `[[`, "smooths"), `[[`, "cobs")
    vapply(cobsL, is.array, NA)
}
(resLok <- ok1cutsLst(resL))
##     1     2     3     4     5     6  ......    15
## FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  ......  TRUE


if(FALSE) {
## e.g.  in R-devel-no-ldouble:
(r1 <- find1cuts(k=1, c1=c1..$c1))$i.n # list --- no longer ok {SparseM::chol -> insufficient space}
(r2 <- find1cuts(k=2, c1=c1..$c1))$i.n # "good"
(r3 <- find1cuts(k=3, c1=c1..$c1))$i.n # 3-vector: only 3x  'i' == 1
}

## if(okCuts) {
mult.fig(15, main = "stirlerr(n, order=k) vs order = k+1")$old.par -> opar
invisible(lapply(resL[resLok], plot1cuts))
## plus the "last" ones {also showing that k=15 is worse here anyway than k=17}
str(r17 <- find1cuts(k=17, n = seq(5.1, 6.5, length.out = if(doExtras) 1500 else 400), c1=c1..$c1))
plot1cuts(r17) # no-ldouble is *very* different than normal:  k *much better* than k+1

(s.find17.fil <- paste0("stirlerr-find1_17_", myPlatform(), ".rds"))
saveRDS(r17, file = s.find17.fil) # was  "stirlerr-find1_17.rds"
## }

## ditto  for tau = 0.8  {quantile for cobs()}:
system.time(
    resL.8 <- lapply(setNames(,k.), function(k) find1cuts(k=k, c1=c1..$c1, tau = 0.80))
) # warnings from cobs  {again 12.1 sec}
(resL8ok <- ok1cutsLst(resL.8))

mult.fig(15, main = "stirlerr(n, order=k) vs order = k+1 -- tau = 0.8")
invisible(lapply(resL.8[resL8ok], plot1cuts))
## plus "last" one
str(r17.8 <- find1cuts(k=17, n = seq(5.1, 6.5, length.out = if(doExtras) 1500 else 400), c1=c1..$c1, tau = 0.80))
plot1cuts(r17.8)
par(opar)

n.mn  <- sapply(resL,   `[[`, "i.n", simplify = "array")
n.mn8 <- sapply(resL.8, `[[`, "i.n", simplify = "array")
## of course, only the cobs part differs:
##   ... when I later see errors here ... what's going on??
str(n.mn) # all NULL  for "no-double" !!
dim(n.mn["n.",,])
str(n.mn8["n.", "cobs",]) # (all NA if {cobs} n.a.)

sessionInfo()

form(data.frame(k = k., cutoffs = rev(.stirl.cutoffs("R4.4_0")[-1])[k.],
                t(n.mn ["n.",,]), cobs.80 = n.mn8["n.", "cobs",]))
##     k     cutoffs         spl         low        cobs     cobs.80
## 1   1 17400000.00 15600000.00 15800000.00 15800000.00 15700000.00
## 2   2     3700.00          NA          NA          NA          NA == from tau=0.90 I'd use ~ 5000
## 3   3      200.00      200.00      201.00      206.00      205.00  205
## 4   4       81.00          NA          NA       86.60       86.20   86
## 5   5       36.00       27.20       27.10       27.00       27.10   27
## 6   6       25.00       23.50          NA          NA          NA   23.5
## 7   7       19.00       12.80       12.80       12.80       12.80   12.8
## 8   8       14.00          NA          NA       12.30       12.80   12.3
## 9   9       11.00        8.91        8.95        8.89        8.86    8.9
## 10 10        9.50        9.69          NA        9.72          NA  --skip--
## 11 11        8.80        7.26        7.31        7.32        7.29    7.3
## 12 12        8.25        8.16          NA        8.10        7.76  --skip--
## 13 13        7.60        6.51        6.53        6.51        6.52    6.52
## 14 14        7.10        7.43          NA        7.47        7.43  --skip--
## 15 15        6.50        6.10        6.10        6.08        6.09    6.10
## 16                                                                 --skip--
## 17                                                                   5.25 or something all the way dow

##---- In the  no-ldouble case --- this is very different:
##  k     cutoffs        spl        low       cobs    cobs.80
##  1 17400000.00 8580000.00 8580000.00 8370000.00 8370000.00
##  2     3700.00         NA         NA    5890.00    6180.00
##  3      200.00     159.00     159.00     160.00     159.00
##  4       81.00      87.10      87.10      84.80      86.00
##  5       36.00      23.50      23.50      23.50      23.50
##  6       25.00      23.70      23.70         NA         NA
##  7       19.00      11.50      11.50      11.50      11.50
##  8       14.00         NA         NA      12.90      13.00
##  9       11.00       8.20       8.20       8.21       8.21
## 10        9.50         NA         NA         NA       9.69
## 11        8.80       6.84       6.84       6.85       6.83
## 12        8.25         NA         NA       8.27       7.95
## 13        7.60         NA         NA         NA         NA
## 14        7.10         NA         NA       7.43       7.40
## 15        6.50         NA         NA         NA         NA


## M1 macOS is similar to  no-ldouble : "=="  if equaln
## M1 mac | aarch64-apple-darwin20 | macOS Ventura 13.3.1 | R-devel (2024-01-29 r85841) ; LAPACK version 3.12.0
##     k     cutoffs        spl        low       cobs    cobs.80
## 1   1 17400000.00 8580000.00 8580000.00 8370000.00 8370000.00  ==
## 2   2     3700.00         NA         NA    5900.00         NA  ~=
## 3   3      200.00     159.00     159.00     160.00     159.00  ==
## 4   4       81.00      87.10      87.10      84.80      86.00  ==
## 5   5       36.00      23.50      23.50      23.50      23.50
## 6   6       25.00      23.70      23.70         NA         NA  ==
## 7   7       19.00      11.50      11.50      11.50      11.50
## 8   8       14.00         NA         NA      12.90      13.00  ==
## 9   9       11.00       8.20       8.20       8.21       8.21
## 10 10        9.50         NA         NA         NA       9.69  ==
## 11 11        8.80       6.84       6.84       6.85       6.83
## 12 12        8.25         NA         NA       8.27       7.95  ==
## 13 13        7.60         NA         NA         NA         NA  ==
## 14 14        7.10         NA         NA       7.43       7.40  ==
## 15 15        6.50         NA         NA         NA         NA  ==




## ========== Try a slightly better direct formula ======================================================

## after some trial error:  gamma(n+1) = n*gamma(n) ==> lgamma(n+1) = lgamma(n) + log(n)
## stirlerrD2 <- function(n) lgamma(n) + n*(1-(l.n <- log(n))) + (l.n - log(2*pi))/2


palette("default")

if(do.pdf) { dev.off(); pdf("stirlerr-tst_others.pdf") }

n <- n5m # (1/64 ... 5000) -- goes with st.nM
p.stirlerrDev(n=n, stnM=st.nM, cex=1/4, type="o", cutoffs = cuts, abs=TRUE)
axis(1,at= 2:6, col=NA, col.axis=(cola <- "lightblue"), line=-3/4)
abline(v = 2:6, lty=3, col=cola)
i.n <- 1 <= n & n <= 15 ; cr2 <- adjustcolor(2, 0.6)
lines(n[i.n], abs(relErrV(st.nM[i.n], stirlerr_simpl(n[i.n], "MM2" ))), col=cr2, lwd=2)
legend(20, 1e-13, legend = quote(  relErr(stirlerr_simpl(n, '"MM2"'))), col=cr2, lwd=3, bty="n")


n <- seq(1,6, by= 1/(if(doExtras) 200 else 64))
stM <- stirlerr(mpfr(n, if(doExtras) 512 else 256), use.halves = FALSE)
relE <- asNumeric(relErrV(stM, cbind(stD = stirlerr_simpl(n), stD2 = stirlerr_simpl(n, "MM2"))))
signif(apply(abs(relE), 2, summary), 4)
##               stD      stD2
## Min.    3.093e-17 7.081e-18
## 1st Qu. 2.777e-15 1.936e-15
## Median  8.735e-15 5.238e-15
## Mean    1.670e-14 1.017e-14  .. well "67% better"
## 3rd Qu. 2.380e-14 1.408e-14
## Max.    1.284e-13 9.382e-14

c2 <- adjustcolor(1:2, 0.6)
matplotB(n, abs19(relE), type="o", cex=3/4, log="y", ylim = c(8e-17, 1.3e-13), yaxt="n", col=c2)
eaxis(2)
abline(h = 2^-53, lty=1, col="gray")
smrelE <- apply(abs(relE), 2, \(y) lowess(n, y, f = 0.1)$y)
matlines(n, smrelE, lwd=3, lty=1)
legend("topleft", legend = expression(stirlerr_simpl(n), stirlerr_simple(n, "MM2")),
       bty='n', col=c2, lwd=3, lty=1)
drawEps.h(-(53:48))


## should we e.g., use interpolation spline through sfserr_halves[] for n <= 7.5
## -- doing the interpolation on the  log(1 - 12*x*stirlerr(x)) vs  log2(x)  scale -- maybe ?
stirL <- curve(1-12*x*stirlerr(x, cutoffs = cuts, verbose=TRUE), 1/64, 8, log="xy", n=if(doExtras) 2048 else 512)
## just need "true" values for x = 2^-(6,5,4,3,2) in addition to those we already have at x = 1/2, 1.5, 2, 2.5, ..., 7.5, 8
xM <- mpfr(stirL$x, if(doExtras) 2048 else 256)
stilEM <- roundMpfr(1 - 12*xM*stirlerr(xM, verbose=TRUE), 128)
relEsml <- relErrV(stilEM, stirL$y)
plot(stirL$x, relEsml) # again: stirlerr() is "limited.." for ~ [2, 6]

## The function to interpolate:
plot(asNumeric(    (stilEM)) ~ x, data = stirL, type = "l", log="x")
plot(asNumeric(sqrt(stilEM)) ~ x, data = stirL, type = "l", log="x")
plot(asNumeric( log(stilEM)) ~ x, data = stirL, type = "l", log="x")

y <- asNumeric( log(stilEM))
x <- stirL$x
spl <- splinefun(log(x), y)

plot(asNumeric(log(stilEM)) ~ x, data = stirL, type = "l", log="x")
lines(x, spl(log(x)), col=2)
summary(rE <- relErrV(target = y, current = spl(log(x)))) # all 0 {of course: interpolation at *these* x}

ssmpl <- smooth.spline(log(x), y)
str(ssmpl$fit)# only 174 knots, 170 coefficients
methods(class="smooth.spline.fit") #-> ..
str(yp <- predict(ssmpl$fit, x=log(x)))
lines(x, yp$y, col=adjustcolor(3, 1/3), lwd=3) # looks fine
## but
summary(reSpl <- relErrV(target = y, current = yp$y))
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## -6.35e-10 -7.55e-11  1.10e-12  0.00e+00  7.32e-11  5.53e-10
## which is *NOT* good enough, of course ....



showProc.time()
