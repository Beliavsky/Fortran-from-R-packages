
#### Experiments about  accurate computation of  log4p1p(p) := log(4*p*(1-p))
#### ------------------------------------------------------------------------
#### ../R/norm_f.R  uses this in some qnorm() approximations (around p~=0)

require(DPQ) # log1mexp
if(!dev.interactive(orNone=TRUE)) pdf("log4p1p-exp.pdf")

### log.p = TRUE: --------------------------------------------------
lp <- seq(-800, -1/4, by=1/4)
if(FALSE) ## no difference (between first two):
    lp <- seq(-1e5, -1, by=100)
c2 <- adjustcolor(2, 1/2)
c3 <- adjustcolor(3, 1/2)
plot(lp, log(4) + lp + log1mexp(-lp), type="l", lty=5)
## should be equivalent {as it always takes first branch in log1mexp()} -- but  faster
lines(lp, log(-4*expm1(lp)) + lp, type="l", col=4, lty=3, lwd=4)
## both above are better for lp < ~-700 than these two (which seem equivalent):
lines(lp, log(4*exp(lp)*1-exp(lp)),  col=c2, lwd=3)
lines(lp, log(4*exp(lp)*-expm1(lp)), col=c3, lwd=5)

### log.p = FALSE -----------------------------------------------------------------------------
p <- seq(0,1, by=2^-10)
matplot(p, cbind(log(4*p*(1-p)), log(4)+log(p)+log1p(-p), log1p(-4*(p-1/2)^2)), type="l")
abline(h=0,col="gray")
## relerr: problem only at p ~ 1/2
matplot(p, cbind(log(4*p*(1-p)), log(4)+log(p)+log1p(-p))/log1p(-4*(p-1/2)^2) -1, type="l")

if(!require("Rmpfr")) quit("no")
#---===============---====------------ Rmpfr needed from here -----------------------------

require("sfsmisc")# (*is* in strong dependencies of DPQ), for relErrV()

### log.p = TRUE: --------------------------------------------------
lpM <- mpfr(lp, 128) # (lp are exact)
reM <- relErrV(target = log(-4*expm1(lpM)) + lpM,
               current= cbind(log(4) + lp + log1mexp(-lp),
                              log(-4*expm1(lp)) + lp,
                              log(4*exp(lp)*-expm1(lp)),
                              log(4*exp(lp)*1-exp(lp))))

matplot(lp, asNumeric(reM), type="l") # last one is catastrophe
apply(abs(asNumeric(reM)), 2, summary)# the 3d one also under/overflows
matplot(lp, asNumeric(reM[,1:3]), type="l") # 3rd one (green) gets bad before overflow
matplot(lp, asNumeric(reM[,1:2]), type="l") #
## difference close to lp~0 i.e. p~1:
lp <- seq(-4, -1/128, by=1/128)# <- at first ==> see it's around ln(0.5) = -0.6931472

lp <- -(44500:46500)*2^-16
lpM <- mpfr(lp, 128) # (lp are exact)
reM <- asNumeric(relErrV(target = log(-4*expm1(lpM)) + lpM,
                         current= cbind(log(4) + lp + log1mexp(-lp)
                                      , log(-4*expm1(lp)) + lp
                                      # , log(4*exp(lp)*-expm1(lp))
                                        )))
matplot(lp, reM, type="l") # a "spike" around  log(p) = log(1/2) = -0.693147...

apply(abs(reM), 2, summary)
matplot(lp, abs(reM), type="l", log="y") #--> 2nd is slightly better!

## smooth on log scale, weight large errors  and re-transform:
psmspl <- function(x,y, ...) exp(predict(smooth.spline(x, y=log(y), w=sqrt(y), ...), x)$y)
smE <- apply(abs(reM), 2, psmspl, x=lp, df=28)
matlines(lp, smE, lwd = 4)
## -- the 2nd is slightly better in the worst region, i.e. close ln(2)



## log.p = FALSE -----------------------------------------------------------------------------
## =============
## see above
p <- seq(0,1, by=2^-10)
## matplot(p, cbind(log(4*p*(1-p)), ........)
pM <- mpfr(p, 128)
reM.<- relErrV(target = log(4*pM*(1-pM)),
               current= cbind(simp = log(4*p*(1-p)),     sumL = log(4)+log(p)+log1p(-p),
                              s2log= log(4*p)+log1p(-p), lg1p = log1p(-4*(p-1/2)^2)))
matplot(p, abs(asNumeric(reM.)), type="l")
## --> *spike* at  p=1/2, the 2nd is *clearly* the worst
apply(abs(asNumeric(reM.)), 2, summary)
##                 simp         sumL        s2log         lg1p
## Min.    0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
## 1st Qu. 1.908159e-17 9.297852e-17 5.160503e-17 1.908159e-17
## Median  3.872176e-17 2.318870e-16 1.269260e-16 3.908697e-17
## Mean    4.017137e-17 6.117178e-14 2.551975e-14 4.066526e-17
## 3rd Qu. 5.784552e-17 8.281096e-16 5.657260e-16 5.857990e-17
## Max.    1.093768e-16 2.425313e-11 4.851490e-12 1.117483e-16
##
## 'sumL' is the worst, 's2log' is close;
## drop "sumL" as worst -- zoom into the others
matplot(p, abs(asNumeric(reM.)[,-2]), type="l")# again 2nd, i.e. "s2log"

## only the 2 good ones,  simp, lg1p :
matplot(p, abs(asNumeric(reM.)[,c(1,4)]), type="l", log="y")
## if you look very carefully (and also by 'summary' table above):
## ==> the first -- simple direct formula -- is even very slightly better than the log1p

