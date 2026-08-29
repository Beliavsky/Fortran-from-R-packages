library(DPQ)

###---> Automatically find places where qgamma() is not so precise (PR#2214) :
###     For PR#2214, had  '1e-8' below and found quite a bit
##  see /u/maechler/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/qgamma-ex.R ..

## FIXME: Timing ! --- partly these matplot() partly get quite slow ~?
source(system.file(package="DPQ", "test-tools.R", mustWork=TRUE))
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()
showProc.time()

(doExtras <- DPQ:::doExtras())
(sdir <- system.file("safe", package="DPQ")) ## save directory (to read from)

### Nowadays finds cases in a special region for really small p and cutoff 1e-11 :
set.seed(47)
n <- if(doExtras) 100 else 32
res <- cbind(p=1,df=1,rE=1)[-1,]
for(M in 1:(if(doExtras) 20 else 10))
for(p in runif(n)) for(df in rlnorm(n)) {
    r <- 1- pchisq(qchisq(p, df),df)/p
    if(abs(r) > 1e-11) res <- rbind(res, c(p,df,r))
}

### use df in U[0,1]: finds two cases with bound 1e-11
for(p in runif(n)/2) for(df in runif(n)) {
    qq <- qchisq(p, df)
    if(qq > 0 && p > 0) {
        r <- 1- pchisq(qq, df) / p
        if(abs(r) > 1e-11) res <- rbind(res, c(p,df,r))
    }
}

### now df very close to 0 : ==> finds more cases
for(p in sort(c(runif(64)/2, exp(-(1+rlnorm(256))))))
    for(df in 2^-rlnorm(256, mean=2, sdlog=1.5)) {
    qq <- qchisq(p, df)
    if(qq > 0 && p > 0) {
        r <- 1- pchisq(qq, df) / p
        if(abs(r) > 1e-11) res <- rbind(res, c(p,df,r))
    }
}
showProc.time()

require(graphics)
if(!dev.interactive(orNone=TRUE)) pdf("qgamma-appr.pdf")
eaxis <- sfsmisc::eaxis

showProc.time()
## if(nrow(res) > 0) {
cat("Found inaccurate examples where  pchisq(qchisq(p, df),df) != p\n")
## sort in p, then df:
res <- res[order(res[,"p"], res[,"df"]), ]
rE <- res[,"rE"]
if(nrow(res) > 20) { hist(rE, breaks = 30); rug(rE) }
plot(res[,1:2])##--> quite interesting : all along one curve
## p <= 1/2  and df <= 1 (about) !!
res <- cbind(res, nDig = round(-log10(abs(rE)), 1))
print(res, digits=12)

if(requireNamespace("scatterplot3d")) {
    scatterplot3d::scatterplot3d(res[,1:3], type ='h') ## quite interesting:
    ## the inaccurate (p,df) points are on nice monotone curve !!!
    ## this is *less* revealing
    scatterplot3d::scatterplot3d(res[,c("p","df","nDig")], type ='h')
}
rL <- res[abs(res[,'rE']) > 1e-9,]
rL <- rL[order(rL[,1],rL[,2]),]
rL
plot(rL[,1:2], type = "b", main = "inaccurate pchisq/qchisq pairs")

plot(rL[,1:2], type = "b", log = "x", ylim = range(0, rL[,"df"]),
     xaxt = "n",
     main = "inaccurate pchisq/qchisq pairs"); abline(h = 0, lty=2)
## aha -- a perfect line !!
lines(res[,1:2], col = adjustcolor(1, 0.5))
eaxis(1); axis(1, at = 1/2)

d <- as.data.frame(res)
plot  (df ~ log(p), data = d, type = "b", cex=1/4, col="gray")
points(df ~ log(p), data = as.data.frame(rL), col=2, cex = 1/2)

summary(fm <- lm (df ~ log(p), data = d, weights = -log(abs(rE))))
## R^2 = 0.9998

p0 <- 2^seq(-50,-1, by=1/8)
dN <- data.frame(p  = p0,
                 df = predict(fm, newdata = data.frame(p = p0)))
rE <- with(dN, 1- pchisq(qchisq(p, df),df)/p)
dN <- cbind(dN, rE = rE, nDig = round(-log10(abs(rE)), 1))
print(dN, digits=10)

## } ## only when we find inaccurate regions
showProc.time()


## Oops: another  qgamma() / qchisq() problem:  mostly NaN's == all solved now
curve(qgamma(x, 20), 1e-16,  1e-10, log='x')
curve(qgamma(x, 20), 1e-300, .99 , log='xy') # and add the critical region from above:
           abline(v=c(1e-16, 1e-10), col="light blue")
curve(qgamma(x, 20), 1e-26,  1e-07, log='x')
##-> now using  log=TRUE in same region:
curve(qgamma(x,      20, log=TRUE), -38, -16)## no problem!!
curve(qgamma(exp(x), 20), add=TRUE, col="green3", n=2001)
## had problem here, but no longer !

##--> Further fix for qgamma: when 'x' is very small: use "log=TRUE of log(x)"!

## had bug (gave NaN), but no longer:
(q_12 <- qgamma(1e-12, 20))
all.equal(1e-12, pgamma(q_12, 20), tol=0)# show rel.err (Lnx 64-bit: 4.04e-16)
stopifnot(
    all.equal(1e-12, pgamma(q_12, 20), tolerance = 1e-14)
)


## --- Nice graphic : --- but amazingly *S..L..O..W*

p.qgammaSml <- function(from= 1e-110, to = 1e-5, ylim = c(0.4, 1000),
                        n = 201, k.lab = 3,
                        a1 = c(10, seq(10.1,20, by=.2), 21:105),
                        a2 = seq(110,330, by=10),
                        a3 = seq(350,1600, by=50))
{
    ## Purpose: nice qgamma() lines  ``for small x'' aka p
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 22 Mar 2004, 14:23
    x <- exp(seq(log(from), log(to), length = n))

    op <- par(las=1, lab = c(10,10, 7), xaxs = "i", mex = 0.8)
    on.exit(par(op))
    plot(x, qgamma(x, a1[1]), log="xy", ylim=ylim, type='l', xaxt = "n",
         main = paste("qgamma(x, a) for very small x, a in [",
         formatC(a1[1]),", ",formatC(max(a1,a2,a3)),"] - log-log", sep=''),
         sub = R.version.string)
    lab.x <- pretty(log10(c(from,to)), 20)
    axis(1, at=10^lab.x, lab = paste("10^",formatC(lab.x),sep=''))
    if(is.nan(qgamma(1e-12, 20)))
        text(1e-60, 20, "all  NaN", cex = 2)
    if(!is.finite(qgamma(1e-140, 155)))
        text(1e-240, 5, "all +Inf", cex = 2)

    lines.txt <- function(a.s, col = par("col")) {
        col <- rep(col, length=length(a.s))
        for(i in seq(along=a.s)) {
            qx <- qgamma(x, (a <- a.s[i]))
            if(i %% k.lab == 0 &&
               any(ifi <- is.finite(qx) & qx >= ylim[1])) {
                ik <- (i%%(2*k.lab))/k.lab # = 0 or 1
                j <- quantile(which(ifi), c(.02,(1:3)/4+ ik/10, .98))
                ## "segments" around the labels :
                i0 <- 1
                for(jj in j) {
                    ii <- i0:(jj-1)
                    i2 <- jj + -1:1
                    lines(x[ii], qx[ii], col=col[i])
                    lines(x[i2], qx[i2], col=col[i], type = 'c')
                    i0 <- jj+1
                }
                text(x[j], qx[j], formatC(a), col= "gray40", cex = 0.8)
            }
            else
                lines(x, qx, col=col[i])

        }
    }
    oo <- options(warn = -1)
    lines.txt(a1[-1])
    lines.txt(a2, col= 2)
    lines.txt(a3, col= rainbow(length(a3), .8, .8,
                               start = (max(a3)-min(a3))/(1+max(a3))))
    invisible(options(oo))
}

showProc.time()

p.qgammaSml()
p.qgammaSml(1e-300)
p.qgammaSml(1e-300,1e-50, a2= seq(100,360, by=4), a3=seq(350,1500, by=10))

showProc.time()

## The "upper" problematic corner:
p.qgammaSml(1e-19, 1e-3, a2=NULL,a3=NULL, ylim=c(.1,20))
p.qgammaSml(1e-19, 1e-3, a2=seq(1,12, by=.04), ylim=c(.1,20),a3=NULL,k.lab=10)
## now shows the problem (quite well):
## could it be in  pgamma()'s inaccuracy, leading to qgamma() bias ?
aa <- c(seq(9, 22, by=0.25),seq(22.3,40,by=0.4))
caa <- formatC(range(aa))
sfsmisc::mult.fig(2)
curve(pgamma(x, sh=aa[1]), 0.5, 20, log = 'xy', ylim = c(1e-60, .2),
      main = sprintf("pgamma(x, a) for a in [%s,%s]", caa[1],caa[2]))
for(sh in aa) curve(pgamma(x, sh), add = TRUE, col=2)
abline(h=c(1e-15), col="light blue", lty=2)

curve(pgamma(x, sh=aa[1]), 0.5, 20, log = 'xy', ylim = c(1e-15, .8),
      main = sprintf("pgamma(x, a) for a in [%s,%s]", caa[1],caa[2]))
for(sh in aa) curve(pgamma(x, sh), add = TRUE, col=2)
## the "border curve" between "Pearson" and "Continued fraction (upper tail)"
## in  pgamma.c :
curve(pgamma(max(1,x), x), add = TRUE, col=4)
## ==> pgamma() is perfect here {series expansion up to eps_C accuracy}!

aa <- c(seq(9, 22, by=0.25),seq(22.3,40.4,by=0.4))
p.qgammaSml(1e-24, 1e-5, a1=aa, a2=NULL,a3=NULL, ylim=c(.8,8))
## -------- save the above?
aa1 <- c(aa,seq(40.5,90, by=0.5))
p.qgammaSml(1e-60, 1e-5, a1=aa1, a2=NULL,a3=NULL, ylim=c(.9, 16))
aa2 <- c(aa1, seq(91,150, by= 1))
p.qgammaSml(1e-90, 1e-5, a1=aa2, a2=NULL,a3=NULL, ylim=c(.9, 35))
aa3 <- c(aa2, seq(150,250, by= 2), seq(253, 400, by=5))
p.qgammaSml(1e-200, 1e-5, a1=aa3, a2=NULL,a3=NULL, ylim=c(.9, 100))
p.qgammaSml(1e-200, 1e-5, a1=aa3, a2=NULL,a3=NULL, ylim=c(.9, 200),k.lab=9e9)
p.qgammaSml(1e-60,  1e-5, a1=aa3, a2=NULL,a3=NULL, ylim=c(.9, 200),k.lab=9e9)

showProc.time()

## lower a \> 10

curve(qgamma(x, 19),    1e-14,    1e-9, log='x')
curve(qgamma(x, 18),    1e-14,    1e-9, log='x')
curve(qgamma(x, 15),    1e-11,    5e-9, log='x')
curve(qgamma(x, 13),    5e-10,    1e-8, log='x')
curve(qgamma(x, 11),     1e-8,    5e-8, log='x')
curve(qgamma(x, 10.5), 4.2e-8,    6e-8, log='x')
curve(qgamma(x, 10.3),   6e-8,    7e-8, log='x')
curve(qgamma(x, 10.2), 7.1e-8,  7.6e-8, log='x')
curve(qgamma(x, 10.15),7.7e-8,  7.9e-8, log='x')
curve(qgamma(x, 10.14),7.88e-8,7.92e-8, log='x',n=10001)

## no more problems for smaller a!!  here:
curve(qgamma(x, 10.13), 1e-10, 5e-4, log='x',n=20001)
curve(qgamma(x, 10.12), 1e-10, 5e-4, log='x',n=20001)
curve(qgamma(x, 10.1), 1e-10,  5e-4, log='x',n=20001)

showProc.time()

##--- the "+Inf" / premature "0" case:
curve(qgamma(x, 155, log=TRUE), -1500, 0, log='y', n=2001,col=2)
curve(qgamma(x, 1e3, log=TRUE), -1500, 0, log='y', n=2001,col=2)
## now works, but slowly and with kink
curve(qgamma    (x, 1e5, log=TRUE), -3e5,  0, log='y', n=2001,col=2,lwd=3)
curve(qgammaAppr(x, 1e5, log=TRUE), add = TRUE, n=2001, col="blue",lwd=.4)
## --- curves are almost "identical"
## ===> the kink *does* come from the initial approx... hmm

## still "identical"
curve(qgamma    (x, 1e4, log=TRUE), -3e4,  0, log='y', n=2001,col=2)
curve(qgammaAppr(x, 1e4, log=TRUE), add = TRUE, n=2001, col="tomato3")

## now see some difference (approx. has kink at ~ -165)
curve(qgamma    (x, 100, log=TRUE), -200,  0, log='y', n=2001,col=2)
curve(qgammaAppr(x, 100, log=TRUE), add = TRUE, n=2001, col="tomato3")
##
(kk <- 100 * 2/1.24)# 161.29
curve(qgamma    (x, 100, log=TRUE), -1.1*kk,  -.95*kk, log='y', n=2001,col=2)
curve(qgammaAppr(x, 100, log=TRUE), add = TRUE, n=2001, col="tomato3")
abline(v = -kk, col='blue', lty=2)# exactly: kink is at  a * 2 / 1.24 = a / .62
curve(qgammaAppr(x - 100/.62, 100,log=TRUE), -1e-3, +1e-3)

showProc.time()

p.qgammaLog <- function(alpha, xl.f = 1.5, xr.f = 0.4, n = 2001)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 30 Mar 2004, 18:44
    kk <- -alpha / .62 # = (alpha * 2) / (-1.24)
    curve(qgamma(x, alpha, log=TRUE), xl.f*kk,  xr.f*kk, log='y',
          n=n, col=2, lwd=3.6, lty = 4,
          main= paste("qgamma(x, alpha=",formatC(alpha,digits=10),", log = TRUE)"))
    lines(kk, qgamma(kk, alpha, log=TRUE), type = 'h', lty = 3)
    curve(qgamma    (exp(x), alpha),      add = TRUE, col="orange", n=n, lwd= 2)
    curve(qgammaAppr(x, alpha, log=TRUE), add = TRUE, col=3, n=n,lwd = .4)
}

showProc.time()

p.qgammaLog(25)
p.qgammaLog(16)# ~ [-25, -20]
p.qgammaLog(12,   1.2, 0.8)# small problem remaining
p.qgammaLog(11,   1.2, 0.8)# even smaller
p.qgammaLog(10.5, 1.1, 0.9)# even smaller
p.qgammaLog(10.25, 1.1, 0.9)# even smaller
## 2019-08: __nothing__ visible from here on:
p.qgammaLog(10.18, 1.02, 0.98)# even smaller
p.qgammaLog(10.15, 1.02, 0.98)# even smaller
p.qgammaLog(10.14, 1.001, 0.999)# even smaller
p.qgammaLog(10.139, 1.0002, 0.9998)#
p.qgammaLog(10.138, 1.0002, 0.9998)#
p.qgammaLog(10.137, 1.00001, 0.99999)#
p.qgammaLog(10.13699, 1.0000001, 0.9999999)#
p.qgammaLog(10.1369899, 1.0000001, 0.9999999)#
p.qgammaLog(10.1369894, 1.0000001, 0.9999999)#
p.qgammaLog(10.1369893, 1.0000001, 0.9999999)# even smaller at -16.34998

showProc.time()

##-- here is the boundary --- for 64-bit AMD Opteron ---
##                        and for 32-bit AMD Athlon

p.qgammaLog(10.1369892, 1.0000001, 0.9999999)# no more
p.qgammaLog(10.136989, 1.0000001, 0.9999999)#
p.qgammaLog(10.136988, 1.0000001, 0.9999999)#
p.qgammaLog(10.136985, 1.0000001, 0.9999999)#
p.qgammaLog(10.13698, 1.0000001, 0.9999999)#
p.qgammaLog(10.13697, 1.0000001, 0.9999999)#
p.qgammaLog(10.13695, 1.0000001, 0.9999999)#
p.qgammaLog(10.1368, 1.000001, 0.999999)#
p.qgammaLog(10.1365, 1.000001, 0.999999)#
p.qgammaLog(10.136, 1.000001, 0.999999)#
p.qgammaLog(10.125, 1.1, 0.9)# --- see it now
p.qgammaLog(10, 1.2, 0.8)
p.qgammaLog(9)

showProc.time()

## For large alpha: show difference to see problem better
## ---> for alpha  >= 10,  the x problem starts *roughly* at x = -0.8*alpha
##

sfsmisc::mult.fig(2)
curve(qgammaAppr(x, 5, log=TRUE), - 8.1, -8, n=2001)
curve(qgammaAppr(x- 5/.62, 5, log=TRUE), -1e-15, 0)

## is the kink from pgamma() ? : no: this looks fine,
curve(pgamma(x, 1e5, log=TRUE), 1, 2e5, log='x', n=2001,col=2)
## and this does too:
curve( dgamma(x, 1e5),           .5e5, 2e5); par(new=TRUE)
curve( dgamma(x, 1e5, log=TRUE), .5e5, 2e5, col=2, yaxt="n")
axis(4,col.axis=2); par(new=TRUE)
curve( pgamma(x, 1e5), .5e5, 2e5, n=2001, col=3); par(new=TRUE)
curve( pgamma(x, 1e5, log=TRUE), .5e5, 2e5, n=2001, col=4); par(new=TRUE)
curve(-pgamma(x, 1e5, log=TRUE,lower=FALSE), .5e5, 2e5, n=2001, col=4)
## all looking nice


x <- 10^seq(2,6, length=4001)
qx <- qgamma(pgamma(x, 1e5, log=TRUE), 1e5, log=TRUE)
plot(x, qx, type ='l', col=2, asp = 1); abline(0,1, lty=3)

showProc.time()

###-------------  Approximations of  qgamma() ------
##

## source("/u/maechler/R/MM/NUMERICS/dpq-functions/qchisqAppr.R")
##--> qchisqAppr()
##--> qchisqWH [ = Wilson Hilferty ]
##--> qchisqKG [ = Kennedy & Gentle's improvements "a la AS 91" ]
## dyn.load("/u/maechler/R/MM/NUMERICS/dpq-functions/qchisq_appr.so")

## Consider the two different implementations of
##  lgamma1p(a) := lgamma(1+a) == log(gamma(1+a) == log(a*gamma(a))  "stable":

lseq <- if(requireNamespace("sfsmisc"))  sfsmisc::lseq  else
	function(from, to, length) exp(seq(log(from), log(to), length.out = length))

if(require("Rmpfr")) { ##---------------- MPFR numbers -------------------------

.mpfr.all.eq <- Rmpfr::all.equal
AllEq <- function(target, current, ...)
    .mpfr.all.eq(target, current, ...,
                 formatFUN = function(x, ...) Rmpfr::format(x, digits = 9))

print(gammaE <- Const("gamma",200)); pi. <- Const("pi",200)
print(a0 <- (gammaE^2 + pi.^2/6)/2)
print(psi2.1 <- -2*zeta(mpfr(3,200)))# == psigamma(1,2) =~ -2.4041138
print(a1 <- (psi2.1 - gammaE*(pi.^2/2 + gammaE^2))/6)

x <- lseq(1e-30, 0.8, length = if(doExtras) 1000 else 125)
x. <- mpfr(x, 200)
xct. <- log(x. * gamma(x.)) ## using  MPFR  arithmetic .. no overflow ...
xc2. <- log(x.) + lgamma(x.)##  (ditto)
print(AllEq(xct., xc2., tol = 0)) # 3.15779......e-57
xct <- as.numeric(xct.)
stopifnot(exprs = {
    AllEq(xct., xc2., tol = 1e-45)
    AllEq(xct , xc2., tol = 1e-15)
    all.equal(lgamma1p(x), lgamma1p(x, tol= 1e-16), tol=0)
    ## -> no difference; i.e., default tol = 1e-14 seems fine enough!
})
showProc.time()

m.appr <- cbind(log(x*gamma(x)), lgamma(1+x), log(x) + lgamma(x),
                lgamma1p.(x, k=1, cut=3e-6),
                lgamma1p.(x, k=2, cut=1e-4),
                lgamma1p.(x, k=3, cut=8e-4),
                lgamma1p(x))#, tol= 1e-14), # = default

eMat <- m.appr - xct # absolute error
## Relative errors:
str(reMat. <- m.appr/xct. - 1)
str(reMat <- as(reMat., "array")) # as(., "matrix") fails in older versions

matplot(x,      eMat , log="x",  type="l", lty=1) #-> problematic  log(x) + lgamma(x) for "large"
matplot(x, abs( eMat), log="xy", type="l", lty=1) #-> but good for small; lgamma1p is much better
matplot(x, abs(reMat), log="xy", type="l", lty=1)
abline(v= 3.47548562941137e-08, col = "gray80", lwd=3)#<- the cutoff value of  lgamma1p()
##---> should use earlier cutoff!
## zoom in:

matplot(x, abs(reMat), log="xy", type="l", col=1:7, lty=1,
        lwd=2, xlim=c(8e-9, 1e-3), ylim = c(1e-18, 1e-7), axes=FALSE, frame=TRUE,
        main = expression(lgamma1p(x) == log(Gamma(x+1)) ~~~ "approximations"
                          ~~~ abs(rel.Err(.))))
eaxis(1); eaxis(2)
abline(v= 3.47548562941137e-08, col = "gray80", lwd=3)#<- the cutoff value of  lgamma1p()
abline(h= c(1,2,4)*.Machine$double.eps, lty=3, col="skyblue")
legend("topright", col=1:7, lty=1,lwd=2,
       c("log(x*gamma(x))", "lgamma(1+x)", "log(x) + lgamma(x)",
         "lgamma1p.(x, k=1, c=3e-6)",
         "lgamma1p.(x, k=2, c=1e-4)",
         "lgamma1p.(x, k=3, c=8e-4)",
         "lgamma1p(x)"), bty="n", ncol=2)
abline(v = c(3e-6, 1e-4, 8e-4), col=4:6, lty=2, lwd=1/2)

## FIXME: do the same for the  lgaamma1p_series()

    ## rm(x., xct., xc2., reMat., eMat, AllEq)
    detach("package:Rmpfr")
showProc.time()

} ## end if( MPFR ) ----------------------------------------------------------------


## ../R/qchisqAppr.R -- talks about the "small shape" qgamma() approximation
## -----------------  --> .qgammaApprBnd() :
curve(.qgammaApprBnd, 1e-18, 1e-15, col=2)
abline(h=0, col="gray70", lty=2)
eps.c <- .Machine$double.eps
axis(3,at=(1:3)* eps.c,
     label=expression(epsilon[c], 2*epsilon[c], 3*epsilon[c]))
(rt.b <- uniroot(.qgammaApprBnd, c(1,3)*eps.c, tol=1e-12))
rt.b$root ## 3.954775e-16
rt.b$root / eps.c ## 1.781072
##==> for a < 1.781*eps, bnd > 0 ==> we have  log(p) < bnd  for all p
## otherwise, we should effectively 'test'
curve(.qgammaApprBnd, 1e-16, 1e-10, log="x", col=2)
showProc.time()


## source  ("/u/maechler/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/qgamma-fn.R")
## ##--> qchisqAppr.R() -- which has 'kind = ' argument!
## ##--> qgamma.R()

p.qchi.appr <-
    function(x, qm= { m <- cbind(qchisq(x, df, log=TRUE),
                                 sapply(knds, function(kind)
                                        qchisqAppr.R(x,df,log=TRUE,kind=kind)))
                      colnames(m) <- c("True", "default", knds[-1])
                      m },
             df,
             knds = list(NULL,"chi.small", "WH", "p1WH", "df.small"),
             call = match.call(), main = deparse(call), log = "y", ...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 25 Mar 2004, 22:08

    col <- c(2,1,3:6)
    lty <- c(1,3,1,1,1,1)
    lwd <- c(2,2,1,1,1,1)
    matplot(x, qm, col=col, lty=lty, lwd=lwd, log = log,
            type = 'l', main = main, ...)
    y0 <- c( .02, .98) %*% par('usr')[3:4]
    if(par("ylog")) y0 <- 10^y0
    legend(min(x), y0, colnames(qm), col=col, lty=lty, lwd=lwd)
    invisible(list(x=x, qm=qm, call=match.call()))
}


pD.chi.appr <- function(pqr, err.kind=c("relative", "absolute"),
                        type = "l", log = "y",
                        lwds = c(2, rep(1, k-1)),
                        cols = seq(along=lwds),
                        ltys = rep(1,k),
                        ...)
{
    ## Purpose: Plot Difference from "True" qchisq()
    ## ----------------------------------------------------------------------
    ## Arguments: pqr: a list as resulting from p.chi.appr()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 31 Mar 2004, 09:38
    err.kind <- match.arg(err.kind)
    if(!is.list(pqr) || !is.numeric(k <- ncol(pqr$qm)-1) || k <= 0)
        stop("invalid first argument 'pqr'")
    with(pqr, {
        err <- abs(if(err.kind == "relative")
                   (1- qm[,-1] / qm[,1]) else (qm[,-1] - qm[,1]))
        matplot(x, err, ylab = "",
                main = paste(err.kind,"Error from", deparse(call)),
                type= type, log= log, lty=ltys, col=cols, lwd=lwds, ...)
        legend(par("xaxp")[1], par("yaxp")[2], colnames(qm)[-1],
               lty=ltys, col=cols, lwd=lwds)
    })
}

## if(FALSE)		# you can manually set
##     do.pqchi <- TRUE    # before source()ing this file
## if(!exists("do.pqchi") || !is.logical(do.pqchi))
##     do.pqchi <-  interactive()

## if(do.pqchi) { #------------- FIXME look at speed up or cache indeed !?? <<<<<

## pqFile <- "/u/maechler/R/MM/NUMERICS/dpq-functions/pqchi.rda"
## ## ls -l  ................... 1325446 Nov  2  2009 pqchi.rda
## if(file.exists(pqFile)) {
##     attach(pqFile) ## it loads more than we create here __FIXME__
##     print(ls(2, all.names=TRUE))
## ##  [1] "pq.1"     "pq.25"    "pq.25."   "pq.31"    "pq.33"    "pq.33."   "pq.33.2"
## ##  [8] "pq.33.3"  "pq.33.4"  "pq.5"     "pq.5."    "pq1"      "pq1."     "pq1.95"
## ## [15] "pq1.95."  "pq1.95.2" "pq10"     "pq10."    "pq10.2"   "pq100"    "pq2"
## ## [22] "pq2."     "pq2.05"   "pq2.05."  "pq2.05.2" "pq2.5"    "pq2.5."   "pq2.5.2"
## ## [29] "pq200"    "pq2L"     "pq4"      "pq4."     "pq4.2"    "pqFile"
## }

showProc.time()
x <- seq(-300, -.01, length=501)

all.equal(qchisqAppr.R(x, 200, log=TRUE),
          qchisqAppr  (x, 200, log=TRUE),tol=0)
## 4.48 e-16 / TRUE (Opteron)

all.equal(qchisqAppr.R(x, 2, log=TRUE),
          qchisqAppr  (x, 2, log=TRUE),tol=0)
## 3.90 e-16 / TRUE (Opteron)

all.equal(qchisqAppr.R(x, 0.1, log=TRUE),
          qchisqAppr  (x, 0.1, log=TRUE),tol=0)
## 7.15 e-15 / 1.179e-8 !!!!! (Opteron)

pq200 <- p.qchi.appr(x = seq(-300, -.01, length=501), df = 200)
pq100 <- p.qchi.appr(x = seq(-160, -.01, length=501), df = 100)
## after (slow) computing, quickly repeat:
with(pq200, p.qchi.appr(x=x, qm=qm, call=call))
with(pq100, p.qchi.appr(x=x, qm=qm, call=call))

## this "hangs forever" -- before I introduced 'maxit' (for 'nu.small'):
pq10  <- p.qchi.appr(x = seq(-12, -.005, length=501), df = 10)
## want to see the jump:
pq10. <- p.qchi.appr(x = seq(-10, -4, length=501), df = 10)
pq10.2<- p.qchi.appr(x = seq(-8.5,-7.5, length=501), df = 10)
with(pq10.2, p.qchi.appr(x=x, qm=qm, call=call))


pq2.5  <- p.qchi.appr(x = seq(-3.4, -.01, length=501), df = 2.5)
pq2.5. <- p.qchi.appr(x = seq(-2.1, -1.8, length=901), df = 2.5)#the jump
## what about p1WH (which is fantastic for df=2)?
pq2.5.2<- p.qchi.appr(x = seq(-0.5, -1e-3, length=901), df = 2.5)
with(pq2.5, p.qchi.appr(x=x, qm=qm, call=call))
with(pq2.5.2, p.qchi.appr(x=x, qm=qm, call=call))
pD.chi.appr(pq2.5.2)# nothing special

pq2.05  <- p.qchi.appr(x = seq(-3.4, -.01, length=501), df = 2.05)
pq2.05. <- p.qchi.appr(x = seq(-2.5, -1.5, length=901), df = 2.05)#the jump
## ^^ the jump from  chi.small to WH is much too late here
## what about p1WH (which is fantastic for df=2)?
pq2.05.2<- p.qchi.appr(x = seq(-0.4, -1e-5, length=901), df = 2.05)
pD.chi.appr(pq2.05.2)                   # p1WH is starting to become better
                                        # and the jump (WH -> p1WH) is too late

with(pq2.05, p.qchi.appr(x=x, qm=qm, call=call))
with(pq2.05.2, p.qchi.appr(x=x, qm=qm, call=call))


## Here, all are 'ok' (but "nu.small"):
pq2L <- p.qchi.appr(seq(-300, -.01, length=201), df = 2)

pq2  <- p.qchi.appr(x = seq(-5, -.001, length=501), df = 2)
pq2. <- p.qchi.appr(x = seq(-2.5, -1, length=901), df = 2)
with(pq2., p.qchi.appr(x=x, qm=qm, call=call))

pq4   <- p.qchi.appr(x = seq(-8,    -0.01, length = 501), df = 4)
summary(warnings())
pq4.2 <- p.qchi.appr(x = seq(-0.1, -1e-04, length = 901), df = 4)

pq1.95  <- p.qchi.appr(x = seq(-3., -.01, length=501), df = 1.95)
pq1.95. <- p.qchi.appr(x = seq(-2.2, -1.5, length=901), df = 1.95)#the jump -1.57
## ^^ the jump from  chi.small to WH is *much* too late here
## what about p1WH (which is fantastic for df=2)?
pq1.95.2<- p.qchi.appr(x = seq(-0.4, -1e-7, length=901), df = 1.95)
pD.chi.appr(pq1.95.2)                   # p1WH is starting to become better
                                        # and the jump (WH -> p1WH) is too late

with(pq1.95, p.qchi.appr(x=x, qm=qm, call=call))
with(pq1.95.2, p.qchi.appr(x=x, qm=qm, call=call))


pq1  <- p.qchi.appr(x = seq(-4, -.001, length=501), df = 1)
pq1. <- p.qchi.appr(x = seq(-1.8, -.5, length=901), df = 1)
with(pq1., p.qchi.appr(x=x, qm=qm, call=call))

pq.5  <- p.qchi.appr(x = seq(-1.5, -.001, length=501), df = 0.5)
pq.5. <- p.qchi.appr(x = seq(-0.8, -0.2, length=901), df = 0.5, ylim=c(.04,.6))

pq.33  <- p.qchi.appr(x= seq(-0.9, -.001,length=501), df= 0.33)
pq.33. <- p.qchi.appr(x= seq(-0.4, -0.02,length=901), df= 0.33)
pq.33.2<- p.qchi.appr(x= seq(-0.4, -0.2, length=901), df= 0.33, ylim=c(.15,.60))
with(pq.33.2, p.qchi.appr(x=x, qm=qm, call=call,ylim=c(.15,.45)))

pq.33.3<- p.qchi.appr(x= seq(-0.4, -0.005, length=901), df= 0.33, ylim=c(.15, 4.00))
with(pq.33.3, p.qchi.appr(x=x, qm=qm, call=call))#,ylim=c(.15, 8)))

pq.33.4<- p.qchi.appr(x= seq(-0.003, -1e-6, length=901), df= 0.33,ylim=c(5,25))
with(pq.33.4, p.qchi.appr(x=x, qm=qm, call=call,ylim=c(5,25)))

## nu <= 0.32 is the "magic" border of  Best & Roberts

pq.31  <- p.qchi.appr(x = seq(-0.45,-.010, length=501), df = 0.31)
with(pq.31, p.qchi.appr(x=x, qm=qm, call=call))

pq.25 <- p.qchi.appr(x = seq(-0.3, -0.02, length=901), df = 0.25)

pq.1 <- p.qchi.appr(x = seq(-0.16, -0.01, length=901), df = 0.1)
with(pq.1, p.qchi.appr(x=x, qm=qm, call=call))
showProc.time()

## if(!file.exists(pqFile)) # don't overwrite for now (as it contains pq2L ,
## save(list=ls(pat="^pq"), file = pqFile)

##}## end if(do.pqchi){ only if interactive } ======================================

pD.chi.appr(pq2L, "abs")
pD.chi.appr(pq2L, "rel")
##  --> want only much smaller x-range:
pD.chi.appr(pq2,"abs")#--> fantastic  p1WH
pD.chi.appr(pq2)      # (ditto)

pD.chi.appr(pq4.2)# p1WH: only at very right
pD.chi.appr(pq4.2, xlim=c(-.016,0))# p1WH: only at very right


## no Newton step here, eg:
(qgA100 <- qgammaAppr(1e-100, 100))
(qg.100 <- qgamma    (1e-100, 100))
all.equal(qgA100, qg.100)
## too much different
dgamma(1e-100, 100, log = TRUE)# -23154.7  i.e.,  "non-log" is 0

qgamma.R(1e-100, 100, verbose = TRUE)#-> final Newton fails!

## but here, the final Newton iterations do work :
x <- qgamma.R(log(1e-100), 100, log = TRUE, verbose = TRUE)
pgamma(x, 100) # = 1e-100 ! perfect
showProc.time()

###--> Use this to devise an improved final Newton algorithm !!!


## From: Prof Brian Ripley <ripley@stats.ox.ac.uk>
## To: skylab.gupta@gmail.com
## Cc: R-bugs@biostat.ku.dk, r-devel@stat.math.ethz.ch
## Subject: Re: [Rd] qgamma inaccuracy (PR#12324)
## Date: Tue, 12 Aug 2008 20:50:50 +0100 (BST)

## This is a really extreme usage.  AFAICS the code works well enough down to
## shape=1e-10 or so, e.g.

qgamma(1e-10, 5e-11, lower.tail=FALSE)
## [1] 0.08237203
## in R 2.9.. .. 2.10.0 -- with an accuracy warning {which is *wrong*!}

## I would be interested to know what substantive problem you were trying to
## solve here that required such values.

## I am pretty sure that a completely different algorithm will be required.

## MM: It looks like this is basis for a new algo:
a <- 1e-14
gammaE <- 0.57721566490153286060651209008240243079 # Euler's gamma
curve(pgamma(x, a, lower.tail=FALSE)/a + log(x) + gammaE, 1e-300, 1e-1, log="",n=1000)
curve(pgamma(x, a, lower.tail=FALSE)/a + log(x) + gammaE - x, 1e-300, 1e-1, log="",n=1000)
## ==> Q = 1 - P = pgamma(x,a, lower=FALSE) ~= a*(-log(x) - gammaE + x - x^2/4)
## i.e.,  Q  ~= -a(log(x) + gammaE { -x + x^2/4 }
##      -Q/a - gammaE ~= log(x)  { -x + x^2/4 }
## ==> x ~=  exp(- (Q/a + gammaE))
## e.g., example below:
Q <- 1e-100; a <- 5e-101
## MM: Find inverse :
str(r.a <- uniroot(function(x) pgamma(x,a,lower.tail=FALSE) - Q,
                   int = c(0.01, 0.1), tol=1e-20))
dput(x0 <- r.a$root) ## 0.0823720296207203
(x1 <- exp(- (Q/a + gammaE)))## 0.07598528 .. not so good
qgammaApprSmallP(Q, a, lower.tail=FALSE)## ~= 0.07598528 -- the same!!

pgamma(x0, a, lower.tail=FALSE) ## 1.00000e-100.
pgamma(x1, a, lower.tail=FALSE) ## 1.03728e-100  ... hmm "close"
##

## MM: -- now look at the bigger picture
p.qg.2a <- function(l2x.min= -15, l2x.max = -100, n = 501,
                    do.offset = FALSE,
                    type = "o", log = "x", cex = 0.6, ...) {
    x.log <- any("x" == strsplit(log,"")[[1]])
    x <- if(x.log) 2^seq(l2x.min, l2x.max, length=n)
             else seq(2^l2x.min,2^l2x.max, length=n)
    if(do.offset)
        plot(x, qgamma(2*x, x, lower.tail=FALSE) - 0.0823720296206873,
             type=type, cex=cex, log=log, ...)
    else plot(x, qgamma(2*x, x, lower.tail=FALSE),
              type=type, cex=cex, log=log, ...)
}
p.qg.2a() # was "very bad" in R <= 2.10.0, now --> 0.082372... "perfect" smooth
## still a little ---but acceptable--- remaining inaccuracy ...zooming in:
p.qg.2a(-43,-55, do.offset=TRUE)
p.qg.2a(,-1024)
p.qg.2a(,-1024, log="", pch=".")## linear in x !!
## zoom in at the limit
p.qg.2a(-30,-1024, do.offset=TRUE, ylim = 1e-11*c(-1,1))
p.qg.2a(-33,-1024, do.offset=TRUE, ylim = 1e-12*c(-1,1))
p.qg.2a(-33,-1024, do.offset=TRUE, ylim = 1e-13*c(-1,1))

a <- 2^-(7:900)
qg <- qgamma(2*a, a, lower.tail=FALSE)
re <- 1-pgamma(qg, a, lower.tail=FALSE)/(2*a)
plot(a, re, log="x", type="b", col=2)
stopifnot(abs(re) < 2e-12) # but really, *should be a bit better

showProc.time()
## For completeness we may write that in due course, but for now (R 2.7.2) I
## suggest just issuing a warning for miniscule 'shape'.

## On Thu, 7 Aug 2008, skylab.gupta@gmail.com wrote:

## > Full_Name:
## > Version: 2.7.1 (2008-06-23)
## > OS: windows vista
## > Submission from: (NULL) (216.82.144.137)
## >
## >
## > Hello,
## >
## > I have been working with various probability distributions in R, and it seems
## > the gamma distribution is inaccurate for some inputs.

## > For example, qgamma(1e-100, 5e-101, lower.tail=FALSE) gives: 1.0. However, it
## > seems this is incorrect; I think the correct answer should be
## > 0.082372029620717283. When I check these numbers using pgamma, I get:

(qg <- qgamma(1e-100, 5e-101, lower.tail=FALSE))# 1 (wrong, originally)
## 0.08237203  now (2009-11-04), i.e. ok
pgamma(qg, 5e-101, lower.tail=FALSE)# now -> 1e-100 : ok

pgamma(0.082372029620717283,5e-101, lower.tail=FALSE)
## 1.0000000000000166e-100.

RE.pqgamma <- function(p, shape, lower.tail = TRUE, log.p = FALSE) {
    ## Relative Error of  pgamma(qgamma(*), ..):
    1 - pgamma(qgamma(p, shape, lower.tail=lower.tail, log.p=log.p),
               shape=shape, lower.tail=lower.tail, log.p=log.p) / p
}
RE.qpgamma <- function(q, shape, lower.tail = TRUE, log.p = FALSE) {
    ## Relative Error of  qgamma(pgamma(*), ..):
    1 - qgamma(pgamma(q, shape, lower.tail=lower.tail, log.p=log.p),
               shape=shape, lower.tail=lower.tail, log.p=log.p) / q
}

## Ok, how extreme can we get -- let a := alpha := shape  --> 0 :
x <- 1e-100
a <- 2^-(7:300)# is still "ok":
plot(a, (re <- RE.pqgamma(x, a, lower.tail=FALSE)), type="b", col=2, log="x")# oops!
a <- 2^-(7:400)
plot(a, (re <- RE.pqgamma(x, a, lower.tail=FALSE)), type="b", col=2, log="x")# oops!
## Oops!
## but, it is clear
qgamma(x, 2^-400, lower.tail=FALSE)## is exactly 0

## -> it goes to 0 quickly . .. zooming in:
curve(qgamma(1e-100, x, lower.tail=FALSE), 1e-110, 1e-70, log="xy", col=2, axes=FALSE)
eaxis(1);eaxis(2)

## from when on is it exactly 0:
uniroot(function(u) qgamma(1e-100, 2^u,lower.tail=FALSE)-1e-315, c(-400, -300))$root
## -341.6941
## => use
a <- 2^-(7:341)
plot(a, (re <- RE.pqgamma(x, a, lower.tail=FALSE)),
     type="l", col=2, log="x")# small glips
## zoom in:
curve(abs(RE.pqgamma(1e-100, x, lower.tail=FALSE)), 2^-341, 1e-90, log="xy", n=2000)
curve(RE.pqgamma(1e-100, x, lower.tail=FALSE), 1e-100, 10e-100, n=2000)
curve(RE.pqgamma(1e-100, x, lower.tail=FALSE), 4e-100, 6e-100, n=2000)
## Ok: at least here is a problem
RE.pqgamma(1e-100, 5e-100, lower.tail=FALSE)# -0.1538

## more general
curve(RE.pqgamma(x, 5*x, lower.tail=FALSE), 1e-100, 1e-20, n=10000, log="x")
## problem *everywhere* , starting quite early: (lesser problem at ~ 1e-16 !)
curve(RE.pqgamma(x, 5*x, lower.tail=FALSE), 1e-25, 1e-15, n=1000, log="x")
curve(RE.pqgamma(x, 5*x, lower.tail=FALSE), 1e-21, 10e-21,n=1000)
curve(RE.pqgamma(x, 5*x, lower.tail=FALSE), 2e-21, 6e-21, n=1000)
curve(RE.pqgamma(x, 5*x, lower.tail=FALSE), 4e-21, 4.5e-21, n=1000)
## and indeed, it's qgamma() that jumps here
curve(qgamma(x, 5*x, lower.tail=FALSE), 4e-21, 4.5e-21, ylim=c(.97, 1.1))

## well, looking at  pgamma(), finally reveals the buglet is there first:
## There's a jump at x = 1 !!
curve(pgamma(x, 1e-30, lower=FALSE), .9999, 1.0001)
curve(pgamma(x, 1e-20, lower=FALSE), .9999, 1.0001)
curve(pgamma(x, 1e-17, lower=FALSE), .9999, 1.0001)
curve(pgamma(x, 1e-15, lower=FALSE), .9999, 1.0001)
curve(pgamma(x, 1e-13, lower=FALSE), .9999, 1.0001)
curve(pgamma(x, 1e-12, lower=FALSE), .9999, 1.0001)# still clearly visible
curve(pgamma(x, 1e-11, lower=FALSE), .9999, 1.0001)# barely visible
## for larger alpha == shape, must zoom in more and more:
curve(pgamma(x, 1e-11, lower=FALSE), .99999, 1.00001)#
curve(pgamma(x, 1e-10, lower=FALSE), .999999, 1.000001)
curve(pgamma(x, 1e-9,  lower=FALSE), .9999999, 1.0000001)
curve(pgamma(x, 1e-8,  lower=FALSE), .99999999, 1.00000001)
curve(pgamma(x, 1e-7,  lower=FALSE), .999999999, 1.000000001)
curve(pgamma(x, 1e-6,  lower=FALSE), .9999999999, 1.0000000001)
curve(pgamma(x, 1e-5,  lower=FALSE), .99999999999, 1.00000000001)
curve(pgamma(x, 1e-4,  lower=FALSE), .999999999999, 1.000000000001)
## now we get close to noise level:
curve(pgamma(x, 1e-3,  lower=FALSE), .9999999999999, 1.0000000000001)
curve(pgamma(x, 1e-2,  lower=FALSE), .99999999999999, 1.00000000000001)

showProc.time()

del.pgamma <- function(a, eps = 1e-13)
{
  ## Purpose: *relative* jump size at x = 1 of pgamma(x, a, lower=FALSE)
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  5 Nov 2009, 16:08
    stopifnot(a > 0, length(a) == 1, eps > 0, length(eps) == 1)
    pp <- pgamma(1+c(-eps,eps), a, lower.tail = FALSE)
    (pp[2] - pp[1])*2/(pp[2] + pp[1])
}

a <- lseq(1e-300, 1e-3, length=400)
dpa <- sapply(a, del.pgamma)
plot(a, -dpa, log="xy", type="l", col=2, yaxt="n");eaxis(2)
## ok, it remains constant all the way to 1e-300
## --> focus
a <- lseq(1e-40, 1e-5, length=400)
dpa <- sapply(a, del.pgamma)
plot(a, -dpa, log="xy", type="l", col=2, axes=FALSE)
eaxis(1, at = 10^-seq(5,40, by=5))
eaxis(2)


xm <- .Machine$double.xmin
pgamma(xm, shape= 1e-20)# is "practically 1" --> *most* qgamma() will be exactly 0
## how close to 1 ?  ---> use upper_tail [possibly on log scale:]
pgamma(xm, shape= 1e-20, lower.tail=FALSE)           # 7.078192e-18
pgamma(xm, shape= 1e-20, lower.tail=FALSE, log=TRUE) # -39.48951

## Where is the 'boundary' (from which on  qgamma() must return 0, since it can't give
##  xm = 2.2e-308 ) ?
a <- 2^-(7:1030)
plot(a, (p <- pgamma(xm, a, lower.tail=FALSE, log=TRUE)),
     cex=.5,type ="l", col=2, log="x")
summary(lm. <- lm(p ~ log(a) + a + I(a^2))) ## coeff. of log(a) *is* 1
summary(lm. <- lm(p ~ offset(log(a)) + a + I(a^2)))

p.l <- pgamma(xm, a, lower.tail=FALSE, log=TRUE) - log(a)
dput(mean(tail(p.l))) ## 6.56218869790132
##=> pgamma(xm, a) ~= log(a) + 6.5621886979 +
## ok, to get this better, now need different a:
al <- seq(1e-300, 1e-3, length=200)
plot(al, (pl <- pgamma(xm, al, lower.tail=FALSE, log=TRUE)) - (log(al)+6.56218869790132),
     cex=.5,type ="l", col=2)
summary(lm.. <- lm(pl ~ offset(log(al) + 6.56218869790132) + 0 + al + I(al^2)))
confint(lm..)
coef(lm..)
##        al    I(al^2)
## -353.8587 20750.8205
##=> pgamma(xm, a) ~= log(a) + 6.5621886979 - 353.858745 * a

## > Similarly, for example: -- now (2009-11-04)  ok
qgamma(1e-100,0.005,lower.tail=FALSE)  # = 109.36757177 instead of 219.5937..
pgamma(109.36757177007101, 0.005, lower.tail=FALSE)# = 1.4787306506694758e-52.

## > This looks completely wrong. The correct value, I think, should be
## > 219.59373661415756. In fact,
pgamma(219.59373661415756, 0.005, lower.tail=FALSE)# = 9.9999999999999558e-101.
## >
## > In fact, when I do the following in R, the results are completely wrong,
## >
a <- 5*10^-(1:40)
z1 <- qgamma(1e-100,a,lower.tail=FALSE)
(y <- pgamma(z1,a,lower.tail=FALSE))
## The value of y that I get should be close to 1e-100, but they are not:
## [1] 1.000000e-100  1.871683e-51  1.478731e-52  1.444034e-53  1.440606e-54
## [6]  1.440264e-55  1.440230e-56  1.440226e-57  1.440226e-58  1.440226e-59
summary(abs(1 - y/1e-100))
plot(a, abs(1 - y/1e-100), log="xy", type="b")
stopifnot(abs(1 - y/1e-100) < 2e-13)# max (32b Linux, P.M) = 4.186e-14

## > The correct values of z1 should be:
z1true <- c(226.97154111939946, 222.15218724493326, 219.59373661415756,
            217.27485383840451, 214.98015408183574, 212.68797118872064,
            210.39614286838227, 208.10445550564617, 205.81289009100664,
            203.52144711679352)
all.equal(z1[1:10], z1true, tol=0)# 1.307e-16 on 32-bit (Pentium M); 1.686e-16 (64b Lnx, 2019)
stopifnot(all.equal(z1[1:10], z1true, tol = 1e-15))
showProc.time()
##>
##> With these values of z1true, we get the expected values:

(y <- pgamma(z1true, a, lower.tail=FALSE))
## [1] 1e-100 1e-100 1e-100 1e-100 1e-100 1e-100 1e-100 1e-100 1e-100 1e-100

## > I am using the precompiled binary version of R, under Windows Vista.
## > -----------
## >> version
## >    _
## > platform       i386-pc-mingw32
## > arch           i386
## > os             mingw32
## > system         i386, mingw32
## > status
## > major          2
## > minor          7.1
## > year           2008
## > month          06
## > day            23
## > svn rev        45970
## > language       R
## > version.string R version 2.7.1 (2008-06-23)
## > ------------
## >
## > So, it seems qgamma is inaccurate for small probability values in the upper
## > tail, when the shape parameter is also small.


###_-- MM:  Still wrong:

(xm <- 2^-1074.9999) # is less than .Machine $ double.xmin == the really x > 0

pgamma(xm, .00001)# 0.992589
qgamma(.99, .00001)##--> NaN  -- should give 0 or "xmin" or so
## FIXME -- ok, now

## but
curve(qgamma(x, .001, lower=FALSE), .4,     .8,   n=1001, log="y")
curve(qgamma(x, 1e-5, lower=FALSE), .002,   .2,   n=1001, log="xy")
curve(qgamma(x, 1e-7, lower=FALSE), 1e-5,   .04,  n=1001, log="xy")
curve(qgamma(x, 1e-12, lower=FALSE), 1e-12, 1e-2, n=1001, log="xy")

## or
curve(qgamma(x, 1e-121, lower=FALSE), 7e-119, 8e-119,
      n=2001, log="y", yaxt="n")
try( # reveals eaxis() bug ? -- for the *subnormal* numbers
 eaxis(2, at = 10^-seq(304,324, by=2))
)

curve(qgamma(x, .001, lower=FALSE), .4, .6, n=1001, log="y")
curve(qgamma(x, .001, lower=FALSE), .5, .55, n=1001, log="y")
## gives a *warning* from axis()  because of subnormal y-range  {was error, fixed in R-devel (2018-08)}
curve(qgamma(x, .001, lower=FALSE), .52, .53, n=1001, log="y")

curve(qgamma(x, .001, lower=FALSE)*1e100, .522, .526, n=1001, log="y")

showProc.time()
