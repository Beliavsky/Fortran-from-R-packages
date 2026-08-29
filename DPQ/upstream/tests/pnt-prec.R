### originally was  ~/R/MM/NUMERICS/dpq-functions/pnt-precision-problem.R
###                 - - - - - - - - - - - - - - - - - - - - - - - - - - -
##==> see also ./t-nonc-tst.R
##               ============
library(DPQ)

stopifnot(exprs = {
    require(graphics)
    require(sfsmisc) # lseq(), p.m(), mult.fig()
})

source(system.file(package="DPQ", "test-tools.R",
                   mustWork=TRUE))# ../inst/test-tools.R
## => showProc.time(), ...  list_() , loadList() ,  readRDS_() , save2RDS()

(doExtras <- DPQ:::doExtras())
## save directory (to read from):
(sdir <- system.file("safe", package="DPQ"))

if(!dev.interactive(orNone=TRUE)) pdf("pnt-precision-2.pdf")
op <- options(nwarnings = 1e5)

x <- 10^seq(2,12, by= .25)
px <- pt(x, df= 0.9, ncp = .01,
         lower.tail=FALSE, log=TRUE) ## 16 warnings
## full precision was not achieved in 'pnt'
plot(x,px, log="x", type="o") #- show catastrophic behavior
## R-devel 2008-11-04: slope is kept, but we have a jump still
##
## extend x  --> "the flattening" (at log(P[ ]) ~= -32) happens still
curve(pt(x, df= 0.9, ncp = .01, lower.tail=FALSE, log=TRUE),
      100, 1e20, log="x")

xs <- lseq(100, 1e6, len = 101)
pxs <- pt(xs, df= 0.9, ncp = .01, lower.tail=FALSE, log=TRUE)
coef(lm(pxs ~ log(xs)))
## (Intercept)      log(x)
##  -1.1484345  -0.8999986
## But (central) pt() has now such problem:
x <- 10^seq(2,17, by=1/4)
mp <- cbind(x,
            pnt = pt(x, df= 0.9, ncp = .01, lower.tail=FALSE, log=TRUE),
            pt  = pt(x, df= 0.9,            lower.tail=FALSE, log=TRUE))
mp
p.m(mp, type="l", log="x", xlab="x", main = "pt(x, df = 0.9,  lower = F, log=TRUE)")
legend("topright", c("ncp = 0.01", "central"), lty=1:2, col=1:2, bty="n")

## even closer: (ncp ~= 0) : --> even bigger problem : pnt(..) becomes -Inf,
## (even in R-devel 2008-11-04)
m2 <- cbind(x,
            pnt = pt(x, df=.9, ncp=1e-4,  lower.tail=FALSE, log=TRUE),
            pnt = pt(x, df=.9, ncp=1e-10, lower.tail=FALSE, log=TRUE),
            pt  = pt(x, df=.9,            lower.tail=FALSE, log=TRUE))
m2
p.m(m2, type="l", log="x", xlab="x", main = "pt(x, df = 0.9,  lower = F, log=TRUE)")
legend("topright", c("ncp = 1e-4", "ncp = 1e-10", "central"), lty=1:3, col=1:3, bty="n")

### Note:  The non-central F ( <==> non-central beta ) is sometimes *WORSE*:
p.t.vs.F <- function(df, ncp, x1 = 1, x2= 1e5, nout = 201, col = c("blue","red"))
{
    ## Purpose: Tail - Comparison of non-central t & non-central F
    ##    Abramowitz & Stegun 26.6.19, p.947

###>> This "non-sense"  non-central F <==> "two-tail" non-central t
###>>   i.e., there's NO direct representation *unless*  ncp = 0

    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Nov 2008, 17:12
    x <- lseq(x1, x2, length=nout)
    m <- cbind(pt=pt(x,           df=df, ncp= ncp  , log=TRUE,lower=FALSE),
               pf=pf(x^2, df1=1, df2=df, ncp= ncp^2, log=TRUE,lower=FALSE))
    matplot(x, m, type="l", log="x", ylab = "log(1 - F*(x))", lwd = 2, col=col,
            main = "upper tail prob.  p*(......, log=TRUE, lower.tail=FALSE)")
    legend("topright",
	   c(sprintf("pt(x,                 df = %g,  ncp=%g  )", df,ncp),
	     sprintf("pf(x^2, df1=1, df2= %g, ncp=%g ^2)",df,ncp)),
	   lty=1:2, col=col, lwd = 2, inset=.01)
    mtext(R.version.string, side=1, line=-1, adj=0.01, cex = 0.6)
    invisible(cbind(x=x, m))
}

p.t.vs.F(3.14, 2)
p.t.vs.F(  2, 20,, 1e9)
p.t.vs.F(0.2, 20,, 1e50)# pnt()_R-devel has small jump; 2.8.0 big jumps to horizontal

## but here (df << ncp),  pf() {ie pnbeta} is clearly better}:
p.t.vs.F(2, 200, 10,1e6)
p.t.vs.F(2, 200, 10,1e9)# breaks down eventually as well

## Small ncp, *however*  one of the two is WAY wrong ! -- "SYSTEMATIC ERROR" :
p.t.vs.F(3, 0.1) ## {and pf() breaks earlier for large x}
p.t.vs.F(3, 0.5) ## still {slightly}
p.t.vs.F(3, 1)   ## no visible problem
## but here, for x -> 0: pf() is systematically larger :
p.t.vs.F(3, 1, 1e-2,100)

showProc.time()

### ---------------------------------------------------------------------------------

## -> pntR() is R-level pnt() function
pt (1e10, df=.9, ncp=1e-10, lower.tail=FALSE, log=TRUE) # >> -Inf
pntR(1e10, df=.9, ncp=1e-10, lower.tail=FALSE, log=TRUE) #  "
## "finis: precision NOT reached"

## back to less extreme ncp:
x <- 2^(0:40)
px <- pt(x, df = 1., ncp = 1, lower.tail = FALSE, log = TRUE)## 14 warnings
plot(x,px, log="x", type="o") #- show "jump" {just small kink: R-devel 25.10.09}

## also the left tail (!) --- this is now (R-devel 25.Oct.09) worse
px <- pt(-x, df = 1.2, ncp = 1, log = TRUE)
plot(x, px, log="x", type="o") #- jump and wrong

## larger df: no jump, still wrong
px <- pt(x, df = 2, ncp = 1, lower.tail = FALSE, log = TRUE) ## 23 warnings
plot(x, px, log="x", type="o") #- show bad behavior
## does my one approximation help here ?
px. <- pntJW39(x, df = 2, ncp = 1, lower.tail = FALSE, log = TRUE)
lines(x,px., col = 2, lwd=2)
## no!  It's worse for large x!

px <- pt(x, df = 2, ncp = 2, lower.tail = FALSE, log = TRUE) ## 23 warnings
plot(x,px, log="x", type="o") #- show bad behavior

## The Mathematica example
## LogLogPlot[ 1 - CDF[NoncentralStudentTDistribution[3, 5], x], {x, 1, 10^6}]
## takes many minutes of computing !!!!
px <- pt(x, df = 3, ncp = 5, lower.tail = FALSE, log = TRUE)
plot(x, px, log="x", type="o") #- show bad behavior
## equivalent {different x-labels}:
plot(log(x), px, type="o")

## BTW:  Very slow convergence here {not in outer tail} -- can clearly be improved
pntR(20, df=0.9, ncp=30, lower.tail=FALSE, log=TRUE, errmax=1e-15)
pntR(20, df=  2, ncp=30, lower.tail=FALSE, log=TRUE, errmax=1e-15)

### "Solve this problem":  Find tail formula empirically
### ----------------------------------------
### log(1 - P(x)) ~= alpha - beta * log(x)
### ----------------------------------------

## Empirically:  alpha = g(ncp) ;  beta = h(df) << see below:  beta == df ( = nu) !!

summary(lm(px ~ log(x), subset=x > 7 & x < 35000))# seems clear, R^2 = 0.9999

showProc.time()

ptRTailAsymp <- function(df,ncp, f.x1 = 1e5, f.x0 = 1, nx = 1000,
                              do.plot=FALSE, verbose = do.plot,
                              ask = do.plot)
{
  ## Purpose: Right tail asymptotic of pt_noncentral()
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Oct 2008, 15:15
    stopifnot(length(df) == 1, is.numeric(df), df > 0)
    stopifnot(length(ncp) == 1, is.numeric(ncp), ncp >= 0)

    n.x0 <- round(.25 * nx) # at this index, we anchor the robust line
    ## Determine an x-range, where pt(*, lower.tail=FALSE, log=TRUE)
    ## is nearly linear
    alp <- .999
    x0. <- qt(alp, df=df, ncp=ncp)##<< it's absurd to use this here,
     ##     --  but unfortunately, qtAppr() is too unreliable here
    it <- 1
    while(!is.finite(x0.) && it < 20) {
        it <- it+1
        alp <- 1 - 2*(1-alp)
        x0. <- qt(alp, df=df, ncp=ncp)
    }
    stopifnot(is.finite(x0.))
    x0 <- x0.* f.x0
    x1 <- x0 * f.x1
    it <- 1; converged <- FALSE
    if(do.plot && ask) { op <- par(ask=TRUE); on.exit(par(op)) }
    repeat {
        xx <- lseq(x0, x1, len = nx)
        lx <- log(xx)
        px <- pt(xx, df=df, ncp=ncp, lower.tail=FALSE, log=TRUE)
        rob.slope <- median(diff(px), na.rm=TRUE) /mean(diff(lx))
        Out <- !is.finite(px)
        SSY <- sum((px[!Out] - mean(px[!Out]))^2)
        if(do.plot) {
            r <- lm.fit(cbind(1, lx[!Out]), px[!Out])
            if(it == 1) { ## first plot
                plot(px ~ lx, xlab = quote(log(x)),
                     xlim = extendrange(r=log(c(x0,x1)), f=0.2),
                     ylim = extendrange(px, f = 0.25),
                     main = sprintf("pt(exp(.), df=%g, ncp=%g, lower.tail=FALSE, log.p=TRUE)",
                                    df, ncp))
                mtext(sprintf("Search for [x0, x1] s.t. pt(*, log=TRUE) is *linear* in it; starting @ [%g, %g]",
                              log(x0), log(x1)))
            } else { # it >= 2 ... subsequent iterations
                lines(px ~ lx)
                lines(lx[!Out], r$fitted, col= "green2")
                abline(a = px[n.x0] - rob.slope*lx[n.x0],
                       b = rob.slope, col = "tomato", lwd=2)
            }
        }
        rob.fitted <- px[n.x0] + rob.slope*(lx - lx[n.x0])
        rob.resid <- px - rob.fitted

        R2.rob <- 1 - sum(rob.resid^2) / SSY
        cat(sprintf("robust R^2 = %12.8f\n", R2.rob))
        if(R2.rob < 0.99) {
            x1 <- x1/10
            if(chg0 <- x0 > x1/2) x0 <- x0/2
            message(c(sprintf("'robust' fit not ok\n ==> x1 := x1 / 10 = %g",
                              x1),
                      if(chg0)sprintf(" -> x0 changed too, to  x0= %g", x0)))
        }
        else { ## reasonable R2.rob
            if(all(sml.res <- abs(rob.resid) < 4*mad(rob.resid))) {
                message("converged in ", it, " iterations")
                converged <- TRUE
                break
            }
            ## else

            if (!sml.res[nx]) { ## largest x[] does NOT have small resid
                x1 <- exp(max(lx[sml.res]))
                if(verbose) cat(sprintf("new   x1: log(x1)= %12g\n", log(x1)))
            }
            else ## throw away small x *ONLY* after large x are fixed
                if (!sml.res[1]) { ## smallest x[] does NOT have small resid
                    x0 <- exp(min(lx[sml.res]))
                    if(verbose) cat(sprintf("new x0: log(x0)= %12g\n", log(x0)))
                }
            else ## no change at both ends.. hmm
                if(R2.rob > .99999 || it > 10) {
                    message("** R^2-pseudo-converged in ", it, " iterations")
                    converged <- NA
                    break
                }
        }
        if((it <- it+1) > 1000) {
            if(verbose) warning("too many iterations!")
            break
        }
    }## end{repeat}

    r <- lm.fit(cbind(1, lx), px)
    R2 <- 1 - sum(r$resid^2) / SSY
    cf <- r$coefficients
    if(do.plot) {
        legend("topright", c(sprintf("final int. [%g, %g]", log(x0), log(x1)),
                             sprintf("coef. = (%g, %g)", cf[1], cf[2])),
               pch=NA, bty="n")
        abline(v = log(c(x0,x1)), col=adjustcolor(1, 1/2), lty=2)
    }
##    browser()
    c(cf, df=df, ncp=ncp, converged=converged,
      x0=x0, x1=x1, alp=alp, x0.999=x0., R2 = R2, R2.rob=R2.rob)
}

p.tailAsymp <- function(r, add.central=FALSE, F.add = FALSE) {
    x01 <- r[c("x0","x1")]
    stopifnot(length(x01) == 2, is.numeric(x01))
    curve(pt(x, df=r["df"], ncp=r["ncp"], log=TRUE,lower.tail=FALSE),
          x01[1] / 100, x01[2]*100, log="x", ylab = "",
          main =
          sprintf(paste("pt(x, df=%g, ncp=%g, log",
                        if(prod(par("mfrow")> 2)) "..)"
                        else "=TRUE, lower.tail=FALSE)", sep=''),
                  r["df"], r["ncp"]))
    if(F.add) {
        ## This "non-sense"  non-central F <==> "two-tail" non-central t
        ## i.e., there's NO direct representation *unless*  ncp = 0
      colF <- "forestgreen"
      curve(pf(x^2, df1=1, df2=r["df"], ncp= r["ncp"]^2,
               log=TRUE,lower.tail=FALSE),
            add = TRUE, col = colF)
      x. <- 10^par("usr")[2]
      text(x., pf(x.^2,df1=1, df2=r["df"],ncp= r["ncp"]^2,
                  log=TRUE,lower.tail=FALSE),
           "pf(x^2, (1,df),  ncp^2, ...)", adj=c(1,-.2), col = colF)
    }
    if(add.central) {
      curve(pt(x, df=r["df"], log=TRUE,lower.tail=FALSE), add=TRUE,
            col = "midnightblue", lwd=2)
    }

    abline(v = x01, col ="gray")

    cL <- "tomato"
    mtext(sprintf("slope = %-10.6g", r[2]), line=-1  , adj=1, col=cL, cex=par("cex"))
    mtext(sprintf("intCpt= %-10.6g", r[1]), line=-2.2, adj=1, col=cL, cex=par("cex"))
    abline(a=r[1], b = r[2] * log(10), col = cL)
}

if(FALSE)
debug(ptRTailAsymp)
r35 <- ptRTailAsymp(df=3, ncp=5)
if(interactive()) x11(type = "Xlib")
r45 <- ptRTailAsymp(df=4, ncp=5, do.plot=TRUE)
r46 <- ptRTailAsymp(df=4, ncp=6,  do.plot=TRUE)
r410 <- ptRTailAsymp(df=4, ncp=10, do.plot=TRUE)
r420 <- ptRTailAsymp(df=4, ncp=20)
p.tailAsymp(r420)
p.tailAsymp(r410)
p.tailAsymp(r46)
p.tailAsymp(r45)

r440 <- ptRTailAsymp(df=4, ncp=40, do.plot=TRUE)## "wrong"
p.tailAsymp(r440)

r.810 <- ptRTailAsymp(df=.8, ncp=10, do.plot=TRUE)
p.tailAsymp(r.810)# tip top

r.518 <- ptRTailAsymp(df=.5, ncp=18, do.plot=TRUE)
p.tailAsymp(r.518)
## better (much faster convergence; better range
r.518 <- ptRTailAsymp(df=.5, ncp=18, do.plot=TRUE, f.x0=1/10,f.x1 = 100)
p.tailAsymp(r.518)

r.110 <- ptRTailAsymp(df=.1, ncp=10, do.plot=TRUE, f.x0=1e-6,f.x1 = 1e10)
p.tailAsymp(r.110) # easy simple

r71 <- ptRTailAsymp(df=7, ncp=1, do.plot=TRUE)
p.tailAsymp(r71)# looks good, even though slope = -6.956
## --converging to *central* t-distrib:
r7.01 <- ptRTailAsymp(df=7, ncp=.01, do.plot=TRUE)
p.tailAsymp(r7.01)# looks ok
r7.001 <- ptRTailAsymp(df=7, ncp=.001, do.plot=TRUE)
p.tailAsymp(r7.001, add.central=TRUE)# --- but the slope is not quite -7 ...
## hmm:
tt <- lseq(10,100, length=1000)
px <- pt(tt, df=7, log=TRUE,lower.tail=FALSE)
plot(px ~ log(tt), pch=".", main = "pt(tt, df=7, log=TRUE, lower.t=F)")
ll <- lm(px ~ log(tt))
summary(ll)
##              Estimate Std. Error t value Pr(>|t|)
## (Intercept)  4.595546   0.004271    1076   <2e-16 ***
## log(tt)     -6.930061   0.001214   -5707   <2e-16 ***
##             ^^^^^^^^^   .......
##--- significantly different from -7

tt <- lseq(100,1e7, length=1000) # much larger range:  t -> oo
px <- pt(tt, df=7, log=TRUE,lower.tail=FALSE)
plot(px ~ log(tt), pch=".",
         main = "pt(tt, df=7, log=TRUE, lower.t=F) -- larger range(tt)")
ll <- lm(px ~ log(tt))
summary(ll)##--> clearly: slope == -7  -- ok

## Left tail
tt <- rev(- lseq(100,1e7, length=1000))
px <- pt(tt, df=4.5, log=TRUE)
plot(px ~ log(-tt), pch=".",
         main = "pt(tt, df=7, log=TRUE, lower.tail=TRUE)")
ll <- lm(px ~ log(-tt))
summary(ll)##--> clearly: slope == -4.5  -- ok

### ====>  (empirical) Theorem:    slope == -df
###                                ============

### and same is true for left tail ---- see below
(rmat0 <- t(sapply(ls(patt="^r.*[0-9]$"), get)))
rmat <- rmat0[ rmat0[,"R2.rob"] > 0.9999 ,]
rmat[,2:3]

p.tailAsymp(r2.30 <- ptRTailAsymp(df= 2, ncp=30)) # very nice
## "Large" ncp  : (the search and resulting slope seem *wrong* ..)
p.tailAsymp(r2.50 <- ptRTailAsymp(df= 2, ncp=50, do.plot=TRUE))

showProc.time()

## Now do a larger set:

nd <- length(df. <- c(.1, .8, 1, 1.5, 2, 4, 10))
nc <- length(nc. <- c(.01, .1, 1, 2, 5))
## for indexing etc:
c.df <- formatC(df., width=1)
c.nc <- formatC(nc., width=1)
c.c <- outer(c.df, c.nc, paste, sep="_")
indR <- function(i.df, i.nc) paste(c.df[i.df], c.nc[i.nc], sep="_")

sfil1 <- file.path(sdir, "pnt-prec-sim1.rds")
if(!doExtras && file.exists(sfil1)) {

    Res <- readRDS_(sfil1)

} else { ## do run the simulation [always if(doExtras)] : ---------------

    r35 <- ptRTailAsymp(df=3, ncp=5)

    if(names(r35)[1] == "") names(r35)[1] <- "intercpt"
    Res <- matrix(NA_real_, nrow = length(r35), ## <-- length of output
                  ncol = nd*nc,
                  dimnames=list(names(r35), c.c))
    for(i.df in seq_along(df.)) {
        df <- df.[i.df]
        cat("\ndf = ", formatC(df)," : \n----------\n")
        for(i.nc in seq_along(nc.)) {
            cat(i.nc,"")
            ncp <- nc.[i.nc]
            r <- try(ptRTailAsymp(df=df, ncp=ncp))
            Res[, indR(i.df,i.nc) ] <- if(inherits(r, "try-error")) NA else r
        }
    }
    attr(Res, "version") <- list(R.version)
    save2RDS(Res, file=sfil1)
}##-- end{ run simulation } --------------------------------------------

dR <- as.data.frame(t(Res))

op <- mult.fig(mfrow=c(nd,nc),marP=-1)$old.par
op2 <- par(cex = par("cex") * 3/4) # or even more
for(i.df in seq_along(df.))
  for(i.nc in seq_along(nc.)) {
    r <- Res[, indR(i.df,i.nc) ]
    if(is.na(r["df"])) { frame() } else p.tailAsymp(r)
                                        # ------------
  }
par(op); par(op2)
showProc.time()


### ====>  Theorem:    slope == -df
###                    ============  [now have prove for ncp=0: central t]
### and same is true for left tail
##  but I see extra discontinuous behavior at x = 0 :

## df = 10  and various  ncp :
curve(pt(x, df=10, ncp=10, log=TRUE), -250,50, ylim=c(-50,0), n=10001)
title("pt(x, df=10, ncp= *, log=TRUE), -250,50, ..)  for various ncp")
##
curve(pt(x, df=10, ncp=0.01, log=TRUE), add=TRUE, n=10001, col=2)
curve(pt(x, df=10, ncp=0.0, log=TRUE), add=TRUE, n=10001, col=2)
##
curve(pt(x, df=10, ncp=0.1, log=TRUE), add=TRUE, n=10001, col="red2")
curve(pt(x, df=10, ncp=0.2, log=TRUE), add=TRUE, n=10001, col="red2")
curve(pt(x, df=10, ncp=0.5, log=TRUE), add=TRUE, n=10001, col="red2")
##
curve(pt(x, df=10, ncp=1, log=TRUE), add=TRUE, n=10001, col="blue2")
curve(pt(x, df=10, ncp=2, log=TRUE), add=TRUE, n=10001, col="blue2")
curve(pt(x, df=10, ncp=5, log=TRUE), add=TRUE, n=10001, col="blue2")


## same, df=5:  here we see increasing problem at x=0 for ncp -> large
curve(pt(x, df=5, ncp=100, log=TRUE), -5000,300, ylim=c(-50,0), n=10001)
curve(pt(x, df=5, ncp=0,   log=TRUE), add=TRUE, n=10001,
      col="blue",lty=3, lwd=3)
for(n. in c(1:8) / 10)
    curve(pt(x, df=5, ncp=n., log=TRUE), add=TRUE, n=2001, col="blue2")
for(n. in c(1:8))
    curve(pt(x, df=5, ncp=n., log=TRUE), add=TRUE, n=2001, col="red2")
for(n. in 10*c(1:8))
    curve(pt(x, df=5, ncp=n., log=TRUE), add=TRUE, n=2001, col="green2")

showProc.time()

1
## From: Michael E Meredith <meredith@easynet.co.uk>
## To: Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: Re: [Rd] Noncentral dt() with tiny 'x' values (PR#8874)
## Date: Thu, 18 May 2006 19:54:09 +0800

## Hello Martin,

## Thanks for the response. I see what you mean about the plots!

## One thing I should also mention is that I'm getting huge numbers of warnings:

## " full precision was not achieved in 'pnt'  "

## I've been hunting for a pbm in my code, but have now found I get it
## with 'power.t.test(*, sd = NULL)
# eg.
(ptt <- power.t.test(20, 1, power=0.8, sd=NULL, tol = 1e-8))
stopifnot(inherits(ptt, "power.htest"), is.list(ptt),
          all.equal(1.0999525, ptt$sd, tol = 1e-7))
## MM: This has been fixed for  R 2.4.0,  with NEWS entry
##
## o	pt() with a very small (or zero) non-centrality parameter could
## 	give an unduly stringent warning about 'full precision was not
## 	achieved'.  (PR#9171)



## I presume this is another related old buglet which is now appearing
## due to the improved reporting of warnings in R 2.3.0.

## Regards,        Mike.


## MM: The note below is
##
##     Fixed for R 2.3.1 --- with NEWS entry
##     -----     --------
## o	dt(x, df, ncp= not.0) no longer gives erratic values for
## 	|x| < ~1e-12.  (PR#8874)

## At 17:57 18/05/2006, you wrote:
## > >>>>> "MikeMer" == meredith  <meredith@easynet.co.uk>
## > >>>>>     on Thu, 18 May 2006 03:52:51 +0200 (CEST) writes:
## >
## >     MikeMer> Full_Name: Mike Meredith
## >     MikeMer> Version: 2.3.0
## >     MikeMer> OS: WinXP SP2
## >     MikeMer> Submission from: (NULL) (210.195.228.29)
## >
## >
## >
## >     MikeMer> Using dt() with a non-centrality parameter and
## > near-zero values for 'x' results
## >     MikeMer> in erratic output. Try this:

tst <- c(1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 0)
dt(tst,16,1)

## >     MikeMer> I get:  0.2381019 0.2385462 0.2296557 0.1851817
## > 0.6288373 3.8163916 (!!)
## >     MikeMer> 0.2382217
## >
## >I get quite different values (several '0', BTW), which just
## >confirms the erratic nature.
## >
## >As often, plots give even a clearer picture:


x <- lseq(1e-3, 1e-33, length= 301)

plot(x, dt(x, df=16, ncp=1),   type = "o", cex=.5, log = "x")
plot(x, dt(x, df=16, ncp=0.1), type = "o", cex=.5, log = "x")
plot(x, dt(x, df= 3, ncp=0.1), type = "o", cex=.5, log = "x")

## >
## >
## >     MikeMer> The 0.238 values are okay, the others nonsense, and
## >     MikeMer> they cause confusing spikes on plots of dt() vs 'x'
## >     MikeMer> if 'x' happens to include tiny values. (Other
## >     MikeMer> values of df and ncp also malfunction, but not all
## >     MikeMer> give results out by an order of magnitude!)
## >
## >I think almost all do, once you start looking at plots like the above.
## >
## >     MikeMer> I'm using the work-around dt(round(x,10),...), but
## >     MikeMer> dt() should really take care of this itself.
## >
## >or actually rather do something more smart; the cutoff at 1e-10
## >is quite crude.
## >
## >Note that this is not a new bug at all; but rather as old as
## >we have dt(*, ncp= .) in R.
## >
## >Thanks for reporting it!
## >Martin Maechler, ETH Zurich

## ===============================================
## Michael E Meredith
## Wildlife Conservation Society (WCS) Malaysia Program
## 7 Jalan Ridgeway, 93250 Kuching, Malaysia
## Fax: +60-82-252 799 Mobile: +60-19 888 1533
## email: meredith@easynet.co.uk  http://www.mered.org.uk
## Program website: http://www.wcsmalaysia.org

