####--------------- bd0()  & ebd0() ------------------------------------------------
##
## was 2nd part of ./stirlerr-tst.R

require(DPQ)
for(pkg in c("Rmpfr", "DPQmpfr"))
    if(!requireNamespace(pkg)) {
        cat("no CRAN package", sQuote(pkg), " ---> no tests here.\n")
        q("no")
    }

n0 <- numeric()
stopifnot(identical(n0, dpois_raw(1, n0))) # gave an error before 2025-05-14

require("Rmpfr")
options(warnPartialMatchArgs = FALSE)

source(system.file(package="DPQ", "test-tools.R", mustWork=TRUE))
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()
##_ options(conflicts.policy = list(depends.ok=TRUE, error=FALSE, warn=FALSE))
require(sfsmisc) # masking  'list_' *and* gmp's factorize(), is.whole()
##_ options(conflicts.policy = NULL)

pks <- c("sfsmisc", "DPQ", "Rmpfr", "DPQmpfr")
sapply(lapply(setNames(,pks), packageVersion), format)

showProc.time()
(doExtras <- DPQ:::doExtras() && !grepl("valgrind", R.home()))

## ebd0 constants: the column sums of "bd0_scale":  log(n / 1024) for all these n
## ---- according to the *comments* in the C code -- so here I test that at least the *sums* are correct
bd0.n <- c(2048,2032,2016,2001,1986,1971,1956,1942,1928,1913,1900,1886,1872,1859,
           1846,1833,1820,1808,1796,1783,1771,1759,1748,1736,1725,1713,1702,1691,
           1680,1670,1659,1649,1638,1628,1618,1608,1598,1589,1579,1570,1560,1551,
           1542,1533,1524,1515,1507,1498,1489,1481,1473,1464,1456,1448,1440,1432,
           1425,1417,1409,1402,1394,1387,1380,1372,1365,1358,1351,1344,1337,1331,
           1324,1317,1311,1304,1298,1291,1285,1279,1273,1266,1260,1254,1248,1242,
           1237,1231,1225,1219,1214,1208,1202,1197,1192,1186,1181,1176,1170,1165,
           1160,1155,1150,1145,1140,1135,1130,1125,1120,1116,1111,1106,1101,1097,
           1092,1088,1083,1079,1074,1070,1066,1061,1057,1053,1049,1044,1040,1036,
           1032,1028,1024)

stopifnot(
    all.equal(bd0.n,
              1024 * exp(colSums(DPQ:::logf_mat)))
) ## on lynne (64-bit, Fedora 32, 2021) they are even *identical*
identical(bd0.n, 1024 * exp(colSums(DPQ:::logf_mat))) # amazingly to me

## Also, the numbers themselves decrease monotonely,
## their differences are close to, but *not* monotone:
diff(bd0.n) # -16 -16 -15 -15 -15 -15 -14 -14 -15 -13 -14 ...
            #                             ^^^^^^^^^^^^^^  (etc)

do.pdf <- TRUE
do.pdf <- !dev.interactive(orNone = TRUE)
do.pdf
if(do.pdf) { pdf.options(width = 9, height = 6.5)# for all {9/6.5 = 1.38 ; 4/3 < 1.38 < sqrt(2) [A4]
    pdf("diff-bd0_tab.pdf")
}

plot(diff(bd0.n), type="b")
c2 <- adjustcolor(2, 1/2)
par(new=TRUE)
plot(diff(bd0.n, differences = 2), type="b", col=c2, axes=FALSE, ann=FALSE)
axis(4, at=-1:2, col=c2, col.axis=c2)

showProc.time()


## use functionality originally in ~/R/MM/NUMERICS/dpq-functions/15628-dpois_raw_accuracy.R
## now -- require(Rmpfr)


## transition till DPQmpfr exports this *and* that version is on CRAN, to ease maintainer("DPQ"):
##                      vvvvvvvvvvvvvvvvvvvvvvvvvvv DPQmpfr 0.3-3 has, but does *not* export dpoisEr()
if(file.exists(ff <- "~/R/Pkgs/DPQmpfr/R/dpoisEr.R")) withAutoprint({ #-------------

source(ff)
str(dpoisEr) # ==>  prBits = 1536 __large__ default: found that  256 was too small for lambda = 1e100  (??)
## but should  prBits  not depend on  lambda

##-- ----- *or* move to vignette ../vignettes/log1pmx-etc.Rnw  <<<<<<<<<<<<<<<

##-------- small lambda --- is  dpois_simpl0() good ?

range(dpE40 <- dpoisEr(40.25, x=0:200)) # integer only: dpois(x, ..) is 0 for non-int !!!
## -2.442882e-16  3.645529e-16  was -4.401959e-16  3.645529e-16
str(attributes(dpE40))
p.dpoisEr(dpE40) # showing errors to be small, as w/ range(..) above

## dpois_simpl0() uses "old" direct formula on original scale: factorial(x)
stopifnot(factorial(170) < Inf,
          factorial(171) == Inf)
xS <- 0:170 # the full range of "sensible" x values for dpois_simpl0
range(dpE40simpl <- dpoisEr(40.25, x=xS, dpoisFUN = dpois_simpl0))
str(attributes(dpE40simpl))

## -1.299950e-13  1.118291e-13
p.dpoisEr(dpE40simpl)
## --> suprising: errors are *very* small  up to  x <= 49, then in in the order of 1e-13
## zoom in [y- range only]
p.dpoisEr(dpE40simpl, ylim = c(-1,1)*4e-16)

## zoom into small x --- integer x only:
range(dpE40simpl2 <- dpoisEr(40.25, x=0:49, dpoisFUN = dpois_simpl0))
## -2.889661e-16  2.076597e-16
p.dpoisEr(dpE40simpl2) #   --- almost all in [-eps, +eps]

## zoom into small x  and use non-integer x:
range(dpE40simpl2d <- dpoisEr(40.25, x=seq(0, 49, by=1/8), dpoisFUN = dpois_simpl0))
## [1] -3.861877e-14  3.080750e-14  == Oops !  blown up to
p.dpoisEr(dpE40simpl2d)

})

### MM: moved much of this to   Rmpfr  vignette:
###
" ~/R/D/R-forge/Rmpfr/pkg/vignettes/gamma-inaccuracy.Rnw "
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## and the R code actually to the R script, also pkg {Rmpfr}
" ~/R/Pkgs/Rmpfr/vignettes/gamma-inaccuracy_R/plot-factErr-def.R "
##=======>  gamma(x) itself suffers from the fact that  exp(y) has a *large* relative error,
##          -------- when  |y| ~ 100 or so, more specifically, the
## relative error of   exp(y) =  |y| * {rel.err(y)} , since
##   exp(((1+ eps)*y) = exp(y) * exp(eps*y) >=  exp(y) (1 + eps*y)  and indeed,
## the inaccuracy of y (i.e. eps)  is blown up by a factor |y|  which is not small here!


## close to over-/underflow -------

### Large lambda == np == M -------

if(do.pdf) { dev.off(); pdf("bd0-ebd0.pdf") }

##-- TODO ----- *or* move to vignette ---> ../vi???????????

LL <- 1e20
dput(x1 <- 1e20 - 2e11) # 9.99999998e+19

(P1 <-         dpois     (x1,       LL)) # was 3.989455e-11; now 5.520993e-98
(P1m <- Rmpfr::dpois(mpfr(x1, 128), LL)) # 5.52099285934214335003128935..e-98
## However -- the ebd0() version
(P1e <- dpois_raw(x1, LL, version="ebd0_v1"))## was 3.989455e-11, but now good!
asNumeric(relErr(P1m, P1))  # 3.218894e-14
asNumeric(relErr(P1m, P1e)) # 3.218894e-14 -- the same, as R's dpois() now *does* ebd0
stopifnot(exprs = {
    all.equal(P1 , 5.520992859342e-98, tol=1e-12)
    all.equal(P1e, P1, tol=1e-12)
    all.equal(P1m, P1, tol=1e-12)
})

options(digits = 9)

## indeed:  regular  bd0()  works "ok" ... and ebd0() now does, too
(bd.1 <- bd0(x1, LL, verbose=2))
## bd0(1e+20, 1e+20): T.series w/ 2 terms -> bd0=200
## [1] 200
(bd.1M <- bd0(x1, mpfr(LL, 128), verbose=2)) # we checked, 128 is sufficient
## bd0(1e+20, 1e+20): T.series w/ 3 terms -> bd0=200
## ---> 199.9999919413334091607468236761591740489
    asNumeric(bd.1 / bd.1M - 1)# -1.82e-17 -- suggests bd0() is really accurate here
stopifnot(abs(bd.1 / bd.1M - 1) < 3e-16,
          all.equal(199.999991941333, bd.1, tolerance=1e-14))

(ebd1 <- sum((ebd.1 <- ebd0(x1, LL, verbose=TRUE))))# fixed since  June 6, 2021
asNumeric(relErr(bd.1M, ebd1)) # 1.603e-16  ... ebd0() slightly *less* accurate than bd0() !!


showProc.time()

### Large x  -- small  np == M ------------------------------------

mpfrPrec <- 1024
mpfrPrec <- 256

yy <-   bd0 (1e307, 10^(-5:1), verbose=TRUE)
yhl <- ebd0 (1e307, 10^(-5:1), verbose=TRUE)
yhlC<- ebd0C(1e307, 10^(-5:1)) ## newly returns data.frame
stopifnot(yy == Inf, yhl$yh == Inf, yhlC == yhl)
yM <- bd0(mpfr(1e307, mpfrPrec), 10^(-5:1))
roundMpfr(range(yM), 12) ##  7.0363e+309 7.1739e+309 -- *are* larger than  DBL_MAX


### Now  *BOTH*  x and lambda are large : ---------------------------------------

##' lseq()--like sequence of length ~= n2ex  around  "round(x)" to use as `np` in  bd0(x, np)
np_bd0 <- function(x, rnd.N = 256, ex2M = 2, n2ex = 256) {
    stopifnot(is.numeric(x), x >= 0, length(x) == 1L)
    l2x <- round(rnd.N*log2(x))/rnd.N
    ## 2 ^ seq()  i.e.,  lseq()--like:
    2^(l2x + seq(-ex2M, ex2M, by=1/(n2ex/2/ex2M)))
}

## (FIXME?? Small loss for ebd0, see below) <<< ???
##  is bd0(<mpfr>, *) really accurate --- (currently, it seems even precBits = 128 is fine)
##  it uses it's own convergent series approximation for |x-np| < .. ????

##' Compute relative errors (wrt MPFR) for "many" bd0() versions
##' @title
##' @param x, np  main arguments for computing bd0*(x, np) =  D_0(x, np) == D_0(x, M)
##' @param mpfrPrec positive integer, typically >= 64, `precBits` bit precision for mpfr-number arithmetic
##' @param delta
##' @param tol_logcf
##' @param ...
##' @param chkVerb  logical, TRUE by default, indicating if some checking + verbosity should happen
##' @param keepMpfr logical, FALSE by default, indication if full precision mpfr result should be storted
##' @return a \code{list()}, ..
##' @author Martin Maechler
bd0ver <- function(x, np = np_bd0(x), mpfrPrec, delta = 0.1, tol_logcf = 1e-14, ...,
                   chkVerb=TRUE, keepMpfr=FALSE) {                              ### passed to log1pmx()
    stopifnot(length(mpfrPrec <- as.integer(mpfrPrec)) == 1,
              !is.na(mpfrPrec), mpfrPrec >= 64,
              x >= 0, np >= 0)
    yy   <-  cbind(bd0      = bd0       (x, np, delta = delta),
                   bd0.l    = bd0_l1pm  (x, np, tol_logcf=tol_logcf, ...),
                   bd0.p1l1 = bd0_p1l1  (x, np, tol_logcf=tol_logcf, ...),
                   bd0.p1.d = bd0_p1l1d (x, np, tol_logcf=tol_logcf, ...),
                   bd0.p1.d1= bd0_p1l1d1(x, np, tol_logcf=tol_logcf, ...))
    yhl  <- ebd0 (x, np)
    yhlC <- ebd0C(x, np)
    if(chkVerb) {
        yhl.  <- ebd0 (x, np, verbose=TRUE)
        yhlC. <- ebd0C(x, np, verbose=TRUE)
        stopifnot(identical(yhl., yhl),
                  identical(yhlC., yhlC))
    }
    epsC <- .Machine$double.eps
    aeq0 <- all.equal(yhl, yhlC, tol = 0)
    aeq4 <- all.equal(yhl, yhlC, tol = 4*epsC)
    if(!isTRUE(aeq4)) warning("the C and R versions of ebd0() differ:", aeq4)
    stopifnot(is.whole(yhl [["yh"]]),
              is.whole(yhlC[["yh"]]))
    yM  <-  bd0(mpfr(x, mpfrPrec),
                mpfr(np,mpfrPrec), verbose=chkVerb)# more accurate ! (?? always ??)
    relE <- relErrV(target = yM, # the mpfr one
                    cbind(ebd0 = yhl [["yh"]] + yhl [["yl"]],
                          ebd0C= yhlC[["yh"]] + yhlC[["yl"]],
                          yy))
    relE <- structure(asNumeric(relE), dim=dim(relE), dimnames=dimnames(relE))
    ## return:
    list(x=x, np=np, delta = delta, bd0=yy, ebd0=yhl, ebd0C=yhlC,
         bd0M=if(keepMpfr) yM, # <- expensive
         aeq0=aeq0, aeq4=aeq4, relE = relE)
}

x. <- 1e307
cbind(log10.M=300:308, ebd0(x., 10^(300:308))) # yl are all 0 ; and bd0(x.,x.) = 0

bd0v.7  <- bd0ver(x., mpfrPrec = 128, chkVerb = FALSE)

bd0v.10 <- bd0ver(x., mpfrPrec = 1024)
stopifnot( all.equal(bd0v.7, bd0v.10, tol=0),
          bd0v.7$aeq0, # even tol=0 equality !
          bd0v.7$aeq4 )
## ==> 256 bit gives the *same* (asNumeric() - double-prec accuracy) as 1024 bits !
## so, at least here
rm(bd0v.10)
showProc.time()

##' Plot the result of bd0ver() -- i.e., relative errors for "many" bd0() versions
##' @title
##' @param bd0v result of  bd0ver()
##' @param dFac
##' @param log a string, "", "x", "y", or "xy"
##' @param type plot `type`, typically "b" or "l"
##' @param add logical
##' @param col.add
##' @param smooth
##' @param f.lowsm lowess() smoothing parameter when adding smooths'
##' @return
##' @author Martin Maechler
p.relE <- function(bd0v, dFac = if(max(np) >= 8e307) 1e10 else 1,
                   log = "x", type="b",
                   add = FALSE, col.add = adjustcolor(k+2, 2/3),
                   smooth = TRUE, f.lowsm = 1/16)
{
    stopifnot(length(x <- bd0v$x) == 1 # for now
            , is.numeric(x), is.numeric(np <- bd0v$np), length(np) > 1
            , is.numeric(dFac), dFac > 0, length(dFac) == 1
            , is.matrix(relE <- bd0v$relE)
            , (k <- ncol(relE)) >= 1
            , sum(iOk <- local({ y <- bd0v$bd0
                ay <- apply(y, 1L, \(ro) any(is.finite(ro) & ro != 0)); ay })) > 1
              )
    np. <- np[iOk]/dFac
  if(add) { ## add only lines for  relE[,"bd0"]
    relE <- relE[iOk, "bd0"] ##--- only (with varying delta, typically)
    abs <- par("ylog")
    if(abs) relE <- abs(relE)
    lines(np., relE, col = col.add, lwd=3)
  } else { ## full
    txtRE <- "relative Errors"
    if(abs <- grepl("y", log)) {
        txtRE <- paste0("|", txtRE, "|")
        relE <- abs(relE)
        yli <- pmax(2^-54, range(relE[iOk,], finite = TRUE))
    } else
        yli <- range(relE[iOk,], finite = TRUE)
    ## */dFac : otherwise triggering axis() error
    ## log - axis(), 'at' creation, _LARGE_ range: invalid {xy}axp or par; nint=5
    ## 	 axp[0:1]=(1e+299,1e+308), usr[0:1]=(7.28752e+298,inf); i=9, ni=1
    pc <- 1:k
    matplot(np., relE[iOk,], type=type, log=log, pch=pc, col=1+pc,
            main = paste(txtRE, "WRT  bd0(<mpfr-accurate>)"),
            xlim = range(np)/dFac, # show full range
            ylim = yli,
            xlab = paste0("np[iOk]", if(dFac != 1) sprintf("/ dFac,  dFac=%g",dFac)),
            ## could use  sfsmisc::pretty10exp(1e10, drop.1=TRUE)
            xaxt="n"); eaxis(1, sub10=3)
    mtext(sprintf("bd0(x, np),  x = %g", x))
    if(k >= 2) legend("topright", colnames(relE), pch=pc, lty=1:2, col=1+pc, bty="n")
    if(abs) {
        abline(h  = 2^(-54:-51), lty = 3, lwd = c(1,1,2,1), col="gray")
        axis(4, at= 2^(-54:-51), las = 1, col.axis="gray", tick = FALSE, cex = 3/4, hadj = +1/2,
             expression(epsilon[C]/4, epsilon[C]/2, epsilon[C], 2*epsilon[C]))
    }
  } # not add
    if(abs && smooth) { # add smoothed relE
        if(add) {
            smRE <- lowess(np., relE[iOk], f = f.lowsm)$y
            lines(np., smRE, col = adjustcolor(col.add, 2/3), lwd = 2.5, lty = "dashed")
            lines(np., smRE, col = adjustcolor("gray44",  1/2), lwd = 4)#lty = 1
        } else {
            smRE <- apply(relE[iOk,], 2L, function(y) lowess(np., y, f = f.lowsm)$y)
            matlines(np., smRE, col = adjustcolor(1+pc,   2/3), lwd = 2.5, lty = "dashed")
            matlines(np., smRE, col = adjustcolor("gray44", 1/2), lwd = 4, lty = 1)
        }
    }
    colD <- if(add) "skyblue" else 2
    jO <-   if(add) 5 else 1
    delta <- bd0v$delta
    inR <- abs(x - np) <= delta * (x + np)
    ##			  -----
    rngTayl <- range(np[inR])/dFac; f.rng <- format(rngTayl, digits = 4)
    message(sum(inR), " np[]/dFac values inside delta-range = [", f.rng[1],", ", f.rng[2],"]")
    abline(v = rngTayl, col = adjustcolor(colD, 2/3), lwd = 2, lty = 3)
    u4 <- par("usr")[3:4] %*% c(jO,64-jO)/64; yUp <- if(abs) 10^u4 else u4
    arrows(rngTayl[1], yUp, rngTayl[2], yUp,  code = 3,
           lwd = 2, col = adjustcolor(colD, 2/3))
    text(rngTayl[1], yUp, substitute(delta == D, list(D = delta)), col=colD, adj = c(5/4, 3/4))
  if(!add) {
    rug(np[!iOk]/dFac, col=2)
    axis(1, at=x/dFac, quote(x), col=2, col.axis=2, lwd=2, line=-1)
  }
}

if(do.pdf) { dev.off(); pdf("p.relE_bd0ver.pdf") }

system.time(bd0v.7.d.40 <- bd0ver(x., mpfrPrec = 128, delta = 0.40, chkVerb = FALSE)
            ) ## larger delta --> some "long taking" Taylor sums

p.relE(bd0v.7)
p.relE(bd0v.7, log = "xy")
p.relE(bd0v.7.d.40, add=TRUE, col.add = adjustcolor("steelblue", .8))



## ==> NOTE:  a whole small (extreme) range where  bd0() is *better* than ebd0() !!!
## correct #{digits}
with(bd0v.7, cbind(log2.lam = log2(np), np, round(-log10(abs(relE)), 1))) ## around 2^[1018, 1021]


with(bd0v.7, stopifnot(yhl[["yl"]] == 0)) # which is not really good and should maybe change !
## Fixed now : both have 4 x Inf and then are equal {but do Note relE difference above!}
stopifnot(all.equal(bd0v.7$ebd0[["yh"]],
                    bd0v.7$bd0[,"bd0"], tolerance = 4 * .Machine$double.eps))
showProc.time()


## bd0() and  ebd0()  are _Inf_  for first three lambda's .. but they *must* be as truly > DBL_MAX
## At the 4th value, Llam = 2^993, bd0() no longer overflows; _FIXME_ ebd0() shld *not* overflow: yM = 1.7598e+308
summary(Llam <- 2^c(990:1023, 1024 - 1e-12))
bd0M <- bd0ver(x., Llam, mpfrPrec = 256, keepMpfr=TRUE)
with(bd0M, data.frame(log2.L = log2(np), bd0 = bd0, ebd0, bd0M. = format(bd0M, digits=8)))

matplot(log2(Llam), with(bd0M, cbind(bd0 = bd0/x., yh=ebd0[["yh"]]/x., asNumeric(bd0M/x.))),
        type="o", ylab = "*bd0*(x., L) / x.", pch=1:3,
        main= paste("ebd0(x., Lam) and bd0(*) for x=",format(x.)," *and* larg Lam"))
abline(h=0, lty=3, col=adjustcolor("gray20", 1/2))
axis(1, at=log2(x.), labels="log2(x.)", line=-1, tck=-1/16, col=2, col.axis=2)
legend("topright", c("bd0()", "ebd0()", "MPFR bd0()"), bty="n", lty=1:3, pch=1:3, col=1:3)
dMax <- .Machine$double.xmax
abline(h = dMax / x., col=4, lty=3) ; ux <- par("usr")[1:2]
text(c(ux %*% c(3,1))/4, dMax/x., pos=3,
     sprintf("bd0(.) > DBL_MAX = %.5g", dMax), col=4)
showProc.time()

bd0.2 <- bd0ver(x., mpfrPrec = 512, delta = 0.25, keepMpfr=TRUE)
p.relE(bd0.2)#, dFac=1)             #===========
p.relE(bd0.2, log = "xy")

### --- 2025-05 (inspired by ~/R/MM/NUMERICS/dpq-functions/dbinom_Lrg-bug.R ):

## dbinom_raw(x=1.2e+308, n=1.72e+308, p=0.2, q=0.8, give_log=1): -->
str(bd0_fns <- sfsmisc::list_(bd0, ebd0, bd0_l1pm, bd0_p1l1d, bd0_p1l1d1))
stopifnot(sapply(bd0_fns, is.function)) # fails e.g. when have 'bd0' matrix
## x ~= M = np = 3.44e307
## tail(x34.7 <- 1e307*seq(1, 18, by=1/32)) # the last is  Inf ==> have fixed all bd0() functions
## for speed, interactive
tail(x34.7 <- 1e307*seq(1, 18, by=1/ 4))

lBM <- lapply(bd0_fns, \(BD0) BD0(x34.7, 3.44e+307)) ## gave if( <NaN> ) error !
## because of the following error [because it called log1pmx(NaN, ..) giving NaN]; now ok: Inf
bd0_l1pm(Inf, 3.44e+307)
stopifnot(identical(Inf, bd0_l1pm(Inf, 3.44e+307)),
          identical(-Inf, log1pmx(Inf))) # now log1pmx(Inf) |--> -Inf correctly

mBM <- cbind(do.call(cbind, lBM[-2]), ebd0 = with(lBM$ebd0, yl+yh))
cbind(x=x34.7, mBM) # after fixes: too early overflow to  Inf *only* happens for ebd0()
## and  *bd0*(Inf, <finite>)  \--> Inf now
stopifnot(tail(mBM,1) == Inf)

if(!doExtras) # gets too expensive
    quit("no") ## FIXME: do not quit: use less precision, etc

options(digits = 6, width = 130) # width: for tables


## zoom in more:
L.3 <- 2^c(seq(1019, 1021, by=1/128))
if(FALSE) ## too slow
 system.time(bd0.3 <- bd0ver(x., L.3, mpfrPrec = 1024, chkVerb=FALSE))# 10.4 sec (was only 7.3s !?)
system.time(bd0.3 <- bd0ver(x., L.3, mpfrPrec =  256, chkVerb=FALSE)) #  2.25 s
p.relE(bd0.3) # up to 1e-11  rel.error !!
p.relE(bd0.3, log="xy")
system.time(bd0.3.d25 <- bd0ver(x., L.3, delta = .25, mpfrPrec =  256, chkVerb=FALSE)) #  2.14 s
p.relE(bd0.3.d25, add = TRUE)

## different x :
system.time(bd0.2.2e307 <- bd0ver(2e307,  mpfrPrec = 256, delta = 0.15, chkVerb=FALSE)) # 1.44 s
p.relE(bd0.2.2e307, log="xy")
system.time(bd0.2.2e307.d.50 <- bd0ver(2e307,  mpfrPrec = 256, delta = 0.50, chkVerb=FALSE)) # 1.44 s
p.relE(bd0.2.2e307.d.50, add=TRUE)

## less large x .. still same problem:  ebd0() is worse than bd0()
system.time(bd0.2.2e305 <- bd0ver(2e305, np=np_bd0(2e305, ex2M = 7),# <- enlarge x-range
                                  mpfrPrec = 256, delta = 0.40, chkVerb=FALSE)) # 1.9 s
p.relE(bd0.2.2e305, log="xy") # interesting: quite different behavior on the left & right !

## less large x .. still same problem:  ebd0() is worse than bd0()
x <- 1e250
bd0.2.1e250 <- bd0ver(x, mpfrPrec = 256, delta = 0.3, chkVerb=FALSE)
p.relE(bd0.2.1e250, log="xy")
## the pd0_l*() and pd0_p*()  log1pmx() and logcf() using versions are quite *stable*
## *and* differ on the left and right  of  np = x

x <- 1e120
bd0.2.1e120 <- bd0ver(x, mpfrPrec = 256, delta = 0.25, chkVerb=FALSE)
## p.relE(bd0.2.1e120) # still rel.E -8e-13
p.relE(bd0.2.1e120, log = "xy", f.lowsm = 1/8)

x <- 1e20
bd0.2.1e20 <- bd0ver(x, mpfrPrec = 256, delta = 0.15, chkVerb=FALSE)
## p.relE(bd0.2.1e20) #  rel.E -8e-12
p.relE(bd0.2.1e20, log = "xy") #  rel.E -8e-12
bd0.2.1e20.d40 <- bd0ver(x, mpfrPrec = 256, delta = .40, chkVerb=FALSE)
p.relE(bd0.2.1e20.d40, add=TRUE)
## again: quite different behavior left and right of np = x

apply(bd0.2.1e20$relE, 2, quantile)
##              ebd0        ebd0C          bd0
## 0%   -1.92048e-12 -1.92048e-12 -8.57215e-16
## 25%  -4.95045e-16 -4.95045e-16 -1.25294e-16
## 50%  -9.57122e-17 -9.57122e-17 -1.52898e-17
## 75%   2.86728e-16  2.86728e-16  1.12103e-16
## 100%  8.79272e-13  8.79272e-13  7.24720e-16

x <- 1e14
bd0.2.1e14 <- bd0ver(x, mpfrPrec = 256, delta = 0.25, chkVerb=FALSE)
## p.relE(bd0.2.1e14) #  rel.E -2.7e-12
p.relE(bd0.2.1e14, log = "xy")

apply(bd0.2.1e14$relE, 2, quantile)
##              ebd0        ebd0C          bd0
## 0%   -2.68404e-13 -2.68404e-13 -7.64861e-16
## 25%  -1.13517e-16 -1.13517e-16 -1.11198e-16
## 50%   3.98896e-17  3.98896e-17  1.82322e-17
## 75%   2.28878e-16  2.28878e-16  1.43407e-16
## 100%  1.05107e-12  1.05107e-12  9.37795e-16


#      = not soo large ===> see that ebd0() is *better* than bd0() -- outside Taylor range!
x <- 1e9
bd0.2.1e9     <- bd0ver(x, mpfrPrec = 256, chkVerb=FALSE)
bd0.2.1e9.d25 <- bd0ver(x, mpfrPrec = 256, delta = 0.25, chkVerb=FALSE)
## p.relE(bd0.2.1e9,     log = "xy") #  rel.E -2.7e-12
p.relE(bd0.2.1e9.d25, log = "xy") #  bd0() still better _inside_ delta-range; but *worse* outside
bd0.2.1e9.d50 <- bd0ver(x, mpfrPrec = 256, delta = 0.50, chkVerb=FALSE)

bd0.1e9.RE <- cbind(bd0.2.1e9$relE,
                    bd0.d.25 = bd0.2.1e9.d25$relE[,"bd0"],
                    bd0.d.50 = bd0.2.1e9.d50$relE[,"bd0"])
apply(bd0.1e9.RE, 2, quantile) |> print(digits = 3)
##           ebd0     ebd0C       bd0     bd0.l   bd0.p1l bd0.p1.d1  bd0.d.25  bd0.d.50
## 0%   -5.65e-13 -5.65e-13 -5.61e-15 -8.33e-16 -5.37e-16 -7.59e-16 -7.14e-16 -6.34e-16
## 25%  -1.54e-16 -1.54e-16 -1.61e-16 -1.12e-16 -7.88e-17 -1.26e-16 -1.41e-16 -1.38e-16
## 50%  -1.63e-17 -1.63e-17 -4.60e-18 -3.02e-17  3.53e-17  2.80e-18 -2.03e-17 -2.38e-17
## 75%   1.22e-16  1.22e-16  1.74e-16  8.25e-17  1.63e-16  1.63e-16  1.21e-16  1.09e-16
## 100%  2.33e-13  2.33e-13  5.04e-15  5.98e-16  1.13e-15  1.13e-15  8.91e-16  6.08e-16

apply(abs(bd0.1e9.RE), 2, quantile) |> print(digits = 3)
##          ebd0    ebd0C      bd0    bd0.l  bd0.p1l bd0.p1.d1 bd0.d.25 bd0.d.50
## 0%   6.14e-19 6.14e-19 6.14e-19 2.94e-19 2.94e-19  1.24e-18 6.14e-19 6.14e-19
## 25%  5.68e-17 5.68e-17 7.08e-17 5.20e-17 5.65e-17  7.52e-17 5.90e-17 5.73e-17
## 50%  1.34e-16 1.34e-16 1.70e-16 1.04e-16 1.33e-16  1.49e-16 1.30e-16 1.30e-16
## 75%  4.14e-16 4.14e-16 3.63e-16 1.75e-16 2.39e-16  2.61e-16 2.41e-16 2.23e-16
## 100% 5.65e-13 5.65e-13 5.61e-15 8.33e-16 1.13e-15  1.13e-15 8.91e-16 6.34e-16

## ===> bd0(  delta = 0.50)  is "best" here

##     =
x <- 1e6
bd0.2.1e6    <- bd0ver(x, mpfrPrec = 256, chkVerb=FALSE)
bd0.2.1e6.30 <- bd0ver(x, mpfrPrec = 256, delta = .30, chkVerb=FALSE)
## p.relE(bd0.2.1e6) #  rel.E +- 4e-12
p.relE(bd0.2.1e6,    log="xy", f.lowsm = 1/8)
p.relE(bd0.2.1e6.30, log="xy", f.lowsm = 1/8)

apply(abs(bd0.2.1e6.30$relE), 2, quantile)
## is not ideal summary, as "left" and "right" differ so much


x <- 1e3
bd0.2.1e3 <- bd0ver(x, mpfrPrec = 256, chkVerb=FALSE)
## p.relE(bd0.2.1e3) #  rel.E +- 4e-12
p.relE(bd0.2.1e3, log="xy") #  rel.E +- 4e-12
bd0.2.1e3.d30 <- bd0ver(x, mpfrPrec = 256, delta = 0.30, chkVerb=FALSE)
p.relE(bd0.2.1e3.d30, add = TRUE, f.lowsm = 1/8)
## --> ebd0() is superior in the very "outskirts"
## --> delta = 0.30 is clearly still too small

## - {number of correct digits}:
## matplot(L.x, log10(abs(bd0.2.1e3$relE[,-1])), type="l")
## abline(v = range(np[abs(x - np) <= 0.1 * (x + np)]), col=adjustcolor(2, 2/3), lwd=2, lty=3)

##  now extend to the full "lambda" / np range [slowish !]
str(np.x <- seq(.5, 4*x, length.out=1001))
bd0..1e3 <- bd0ver(x, np.x, mpfrPrec = 256, chkVerb=FALSE)
## this is slow:
## p.relE(bd0..1e3, log="") #  rel.E +- 4e-12
bd0..1e3.d.40 <- bd0ver(x, np.x, delta = 0.40, mpfrPrec = 256, chkVerb=FALSE)


if(do.pdf) { dev.off(); pdf("p.relE_bd0ver_x=1000.pdf") }

p.relE(bd0..1e3, log="y") #  rel.E +- 4e-12
p.relE(bd0..1e3.d.40, add = TRUE, f.lowsm = 1/8) #  rel.E +- 4e-12
mtext(shortRversion(), adj=1, cex=2/3); mtext('p.relE(bd0..1e3, log="y")', adj=0, cex=2/3)

apply(bd0..1e3$relE, 2, quantile)
## even here, bd0(x, np)  is more accurate around  np ~= x
## but in the flanks,  ebd0 is better :
matplot(np.x, log10(abs(bd0..1e3$relE[,-1])), type="l")

if(do.pdf) dev.off()

### all the above till    `` quit("no") ''  is *only* run  if(doExtras)

