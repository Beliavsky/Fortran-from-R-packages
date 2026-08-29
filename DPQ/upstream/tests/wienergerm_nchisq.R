### Testing numerical stability of h() {incl. Taylor around 0}
### ---------------------------    ===

library(DPQ)
## --> pchisq.W(), h0(), h1(), h2() === h(),  etc

mult.fig <- sfsmisc::mult.fig
## and (used via :: below) rrange <- sfsmisc::rrange

### TODO: Use many  stopifnot()  below <<<<<<<<<<<<< TODO >>>

## h() only defined in [0, 1]

if(!dev.interactive(orNone=TRUE)) pdf("wienergerm-h-gnt.pdf")

x <- seq(0, 1, length=101)
all.equal(h0(x), h1(x), tol= 5e-14)## TRUE (386 Linux)

## not only in [0,1] :
plot(h, -4, 1,  n= 2001)
abline(0, 1/6,col=2,lty=2); abline(h=c(0,1/2),v=0:1,col='gray', lty=3)
plot(h, 0, 1.1, n= 2001)## only in [0,1]
abline(0, 1/6,col=2,lty=3); abline(h=c(0,1/2),v=0:1,col='gray',lty=3)
## quite linear already:

## FIXME use h0,h1,h2
plot(h, 0, 0.1,  n= 2001);abline(0, 1/6,col=2,lty=3)
## very linear:
plot(h, 0, 0.01, n= 2001);abline(0, 1/6,col=2,lty=3)
## already breakdown for "h0"; ok for "h1"
plot(h, 0, 1e-4, n= 2001);abline(0, 1/6,col=2,lty=3)
## noise even for "h1" (perfekt for "h2"):
plot(h, 0, 1e-7, n= 2001);abline(0, 1/6,col=2,lty=3)
## Difference: looks very quadratic
x <- seq(0, 1e-7, length=2001)
plot(x, h2(x) - x/6, type='l', col=2)# look at y-scale!


### Look close to 0 :
y <- 2^seq(-50,-11, length=2001)
y <- 2^seq(-24,-11, length=2001)
op <- mult.fig(2)$old.par
for(log in c('', "x")) {
    plot(y,hnt(y,0)/y,  type='l', ylim = 1/6+c(-1,1)*1e-3, log=log)
    lines(y,hnt(y, 2)/y, col=3)
    if(log=='')title("hnt[d](y)/y  &  hnt(y,2)/y {only terms k=0,1,2}")
    else  title("at   y <- 2^ seq(-80,-10, length=2001)")
    ## lines(y,hnt(y,10), col=2)
}; par(op)

y <- 2^seq(-19,-11, length=2001)
op <- mult.fig(2)$old.par
for(log in c('', "x")) {
    plot (y,hnt(y,0)/y - 1/6, type='l', ylim = sfsmisc::rrange(hnt(y,0)/y - 1/6),
          log=log, xlab=NA)
    lines(y,hnt(y,4)/y - 1/6, col=3)
    if(log=='')title("hnt[d](y)/y -1/6 & hnt(y,4)/y - 1/6", xlab = "y")
    else  title("at   y <- 2^ seq(-19,-11, length=2001)", xlab = "y [log scale]")
    ## lines(y,hnt(y,10), col=2)
}; par(op)


### Wiener germ approx of  pchisq(x, df, ncp) -- Testing the Fortran / C code
###                        ------------------                ------------
### ==> ./wienergerm-pchisq-tst.R
###       =======================

nu.lam.expr <- function(df, ncp)
    ## a 'title' string
    substitute(list(nu == df, lambda == ncp),
               list(df = df, ncp = ncp))

### ---> from now on checking  C version only ==> pchisqW() <<<<<
###                            --------------     ---------

p.W.pchisq <- function(df, ncp, n = 513, x = NULL, do.log = TRUE,
                       cex = 0.5, s.col = "purple")
{
    ## Purpose: Plot "Wiener germ" approx. pnchisq()  problem around s=1 <==> x = df+ncp
    ## ---------------------------------------------------------------------- ==========
    ## Author: Martin Maechler, Date: 26 Jan 2004, 17:32

    x0 <- df + ncp
    if(is.null(x))
        x <- seq(0.5*min(df,ncp), 1.5*x0, length = n)## << better ?

    if(do.log) {
        op <- mult.fig(2, main="wienergerm : problem at x ~= ncp+df (df 'small')",
                 marP = c(0,0,0,1.5), quiet = TRUE)$old.par
        on.exit(par(op))
    }
    clr <- c("red", adjustcolor(c("forestgreen","midnightblue"), 3/4, 1/6))
    for(log in c('x', if(do.log)'xy')) {
        plot (x, pchisq (x, df=df, ncp=ncp), col = clr[1], log = log, cex = cex)
        abline(h = 0.5, col = "gray", lty=3)
        lines(x,   pchisqW(x,  df=df, ncp=ncp, var = "f"), col = clr[2], lwd=3)
        lines(x,   pchisqW(x,  df=df, ncp=ncp, var = "s"), col = clr[3], lwd=5)
        points(x0, pchisqW(x0, df=df, ncp=ncp), col = 1, pch=8, cex = 1.5)
        axis(3, at=x0, label=quote(x[0] == nu + lambda))
        if(log == 'x') {
            title(nu.lam.expr(df,ncp))
            legend(x[1], par('usr')[4],c("true","'f'","'s'"),
                   col=clr, pch= c(1,NA,NA), lty=c(NA,1,1), lwd=c(NA,3,5),
                   xjust=0, yjust=1.4, horiz=TRUE)
            op <- par(new=TRUE)
            sx <- sW(x, df=df, ncp=ncp)$s
            plot(x, sx, type = 'l', col= s.col, log=log, axes=FALSE, ylab='')
            x.mx <- 10^par("usr")[2]
            segments(x0, 1, x.mx, 1, col=s.col, lty=2)
            axis(4, col.axis = s.col, col = s.col)
            mtext("s(x)",side=4, col = s.col, adj = 1.02)
        } else {
            title(sprintf('log = "%s"', log))
        }
    }

}# p.W.pchisq()

p.W.pchisq(2 , 1)
p.W.pchisq(1 , 2)
p.W.pchisq(1 , 1)
p.W.pchisq(1 , 0.5)# here 'f' is quite different from 's'
p.W.pchisq(1 , 0.1)
## df < 1 : becomes very bad (for s <= 1); not really good even further up
p.W.pchisq(.5, 0.1)
p.W.pchisq(.1, 0.2, x=2^seq(-50,3, length=2001))
## for larger x : "first" is clearly better than "second"!
p.W.pchisq(.1, 0.2, x=2^seq(-4,20,length=3001))

## df larger --> (always ?)  problem smaller
p.W.pchisq(4 , 1)
p.W.pchisq(3 , 2)

p.W.pchisq(4 , 5)
p.W.pchisq(7 , 3)
p.W.pchisq(7 , 3, x = seq(8,100, length=501))

## ncp >= 80 --- the now (2010++) important region:
p.W.pchisq(7 , 80, x = seq(8,100, length=501))
curve(pchisqW(x, 7, 80), 84,89, n=1024); abline(h = 0.5, lty=2, col="gray50")

## how can I find the region where z < 0 ?
## it's boundaries are where the final sign change should happen, not s <= 1!

p.z.s <- function(df, ncp, n = 513, x = NULL, p.pchi = FALSE)
{
  ## Purpose: Plot "Wiener germ" approx. pnchisq()  problem
  ## ----------------------------------------------------------------------
  ## Author: Martin Mächler, Date: 26 Jan 2004, 21:44

  x0 <- df + ncp
  if(is.null(x))
    x <- seq(0.6*mean(df,ncp), 1.25*x0, length = n)
  else if(is.unsorted(x)) x <- sort(x)

  zfs <- cbind(z.f(x,df,ncp),
               z.s(x,df,ncp))
  zmat <- cbind(z0 (x,df,ncp), zfs)

  ## Select interesting y-range
  yl <- range(zmat)
  ym <- min(zfs)
  yl[1] <- pmax(5*ym,  -0.2,   yl[1])
  yl[2] <- pmin(0.5, 10*(-ym), yl[2])

  ## select x-range from y-range
  xl <- range(x[y.in <- apply(yl[1] <= zmat & zmat <= yl[2], 1, any)])
  ix <- xl[1] <= x & x <= xl[2]
  ## and sub-set
  x <- x[ix]
  zmat <- zmat[ix, ,drop=FALSE]
  zfs  <- zfs [ix, ,drop=FALSE]
  s <- sW(x,df,ncp)$s

  ## really interested in negative z's :
  ineg <- range(which(apply(zfs < 0,  1, any)))# << first and last x-index

  if(p.pchi) {
      op <- mult.fig(2, marP = c(0,0,0,1.5), quiet = TRUE)$old.par
      on.exit(par(op))
      p.W.pchisq(df, ncp, do.log=FALSE)
  }
  matplot(s, zmat, type='l', ylim = yl, xlab = 's(.)', ylab = 'z(.)',
          col = 2:4, lwd = c(1,2,1), lty = 1, main = nu.lam.expr(df,ncp))
  abline(h=0,lty=3)

  usr <- par('usr')
  legend(s[1], usr[4],c("z0","z'f'","z's'"),
         col=2:4, lwd= c(1,2,1), xjust=0, yjust=1.4, horiz=TRUE)
  s0 <- s[ineg]
  x0 <- x[ineg]
  lab <- paste(formatC(s0), paste("x=",formatC(x0)), sep="\n")
  op <- par(mgp=3:1)
  axis(1, at= s0, labels=lab, pos= 0, lwd=2, col="midnight blue")
  par(op)
  list(x.range = x0, s.range = s0)
} ## p.z.s()


p.z.s(3,1, p.p= TRUE)# now shows 'z.s(s ~= 1)' problem, z.f() is fine
p.z.s(1,3, p.p= TRUE)
p.z.s(1,1, p.p= TRUE)
p.z.s(6,1, p.p= TRUE)
## all these have   df+ncp = 16.1 :
p.z.s(16, 0.1, p.p= TRUE)
p.z.s(14, 2.1, p.p= TRUE)
p.z.s(10, 6.1, p.p= TRUE)
p.z.s( 6,10.1, p.p= TRUE)

## (6, 10)  for z.s(), (s ~= 1) problem:
curve(z.s(x,6,10), 16-.05, 16+ .05, ylim = c(0.014,0.018),n=2000)
## "first" doesn't have the problem:
curve(z.f(x,6,10), add = TRUE, n = 2000, col=2)

z.s(16.001, 6,10, verbose=TRUE)
## c(f= 4.33346, s= 1.00004, h= -6.41003e-06, qs=4.33331, z0=3.55211e-05,
##   eta=-3.55025e-05, g=-5538.39, t1=16615.2, t2=-16615.3)
##[1] 0.01601603

p.zs1 <- function(df, ncp,
                  r.ex = c(-8, -5),# depends on (df, ncp) [not much though]
                  n = 513, signs = c(-1,1), type = 'l',
                  xlab = '', ylab = paste("z.s(x0  +/- exp(x)"),
                  main = nu.lam.expr(df,ncp)
                  )
{
  ## Purpose: problematic of z.s() when s --> 1 :
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 27 Jan 2004, 13:33

  x0 <- df + ncp
  cx0 <- formatC(x0)

  xO <- pretty(exp(r.ex))
  xO <- xO[xO > 0]
  xO <- c(xO[1]/2, xO)
  ## x-limits : use "correct" on "+", artificially shifted left by Dx on "-"
  dx <- diff(xl <- rx <- range(r.ex))
  e <- dx/16
  ## final size:   dx + e + dx  ;   e ~= dx/10
  xl[1] <- xl[1] - (e+dx)
  Dx <- dx+e ## = correction for left part

  at <- log(xO)
  at <- list(2*rx[1]-e - at, at)# for (-1 , 1)

  xo <- seq(r.ex[1], r.ex[2], length = n)

  exo <- exp(xo)

  ## Compute everything -> y-layout
  z <- as.list(signs)
  for(is in seq(along=signs)) {
      sig <- signs[is]
      z[[is]] <- z.s(x0 + sig*exo, df, ncp)
  }
  yl <- sfsmisc::rrange(unlist(z), r = 2) # y-limits
  plot(NA,NA, xlim = xl, ylim= yl, type = "n", xaxt = "n",
       xlab=xlab, ylab=ylab, main=main)
  uy <- par('usr')[3:4]
  polygon(x= rx[1]+c(-e,-e,0,0),
          y= uy[c(1:2,2:1)], col = 'gray70', xpd = FALSE)
  mtext(paste(cx0, "+"), side = 1, cex = 1.5, adj = 0.5, line = 0.1)
  op <- par(mgp=c(3,0.5,0))
  axis(3, at = pretty(rx), cex.axis = 0.8)
  par(op)
  mtext(paste("log( x -", cx0,")"), line= 1.2, adj = 0.75)
  for(is in seq(along=signs)) {
      sig <- signs[is]
      lines(if(sig > 0) xo else rev(xo) - Dx,
            y = z[[is]], type = type)
      axis(1, at=at[[is]], labels = formatC(sig*xO))
  }
}

p.zs1(3,1)
p.zs1(6,2)
p.zs1(6,10)
p.zs1(10,10)
p.zs1(10,12)
p.zs1(10,20)## oops! doesn't happen at "30" but slightly afterwards!

p.zs1(10,200, c(-10,-3))# (two regions)!


###-- Find improved z.s() --- need g(1-s)  accurately ----------

tit.zs <- "pnchisq() Wiener germ z.'s' around s ~=1"

curve(gnt(x, 0),-2, 1, n =1001, xlab = "u = 1 - s", main=tit.zs)
curve(gnt(x,10),-2, 1, n =1001, add=TRUE,col=2)
curve(gnt(x, 4),-2, 1, n =1001, add=TRUE,col=3)
curve(gnt(x, 2),-2, 1, n =1001, add=TRUE,col=4)
lines(0, gnt(0,2), type = "h", col="blue3", lty=3)
legend(-1,1.8, paste("g(u)", c("direct",paste(format(c(10,4,2)),"-term",sep=''))),
       lty=1, col = 1:4)

### Look close to 0 :
y <- 2^seq(-50,-10, length=2001)
for(log in c('', "x")) {
    plot(y,gnt(y,0),  type='l', ylim = c(.499,.501), log=log)
    lines(y,gnt(y, 2), col=3)
    if(log=='')title("gnt[d](y)  &  gnt(y,2) {only terms k=0,1,2}")
    else  title("at   y <- 2^ seq(-50,-10, length=2001)")
    ## lines(y,gnt(y,10), col=2)
}; par(op)

## -- now look at relative *error*

p.g <- function(to = 0.1, from = -to, Ns = c(1:4,10), cutoff = 0.2,
                abs.do = FALSE, do.log = abs.do, leg.do = TRUE,
                n.out = 513, main = tit.zs)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 29 Jan 2004, 13:51
    x <- seq(from,to, length=n.out)
    g0 <- ifelse(abs(x) < cutoff, gnt(x,20), gnt(x,0))
    nN <- length(Ns)
    gmat <- matrix(NA, n.out, nN)
    for(j in seq(along=Ns)) { N <- Ns[j]; gmat[,j] <- gnt(x, N) }
    y <- g0 / gmat - 1
    if(abs.do || do.log)
        y <- abs(y)
    if(do.log) { i0 <- y < 2*.Machine$double.eps; y[i0] <- NA }
    matplot(x, y, type = 'l', ylim = sfsmisc::rrange(y), log = if(do.log)"y" else '',
            xlab = '1-s', main = main,
            ylab = paste("rel.Error g*(u) / g(u, n=",deparse(Ns),") - 1",sep=''))
    yU <- par('usr')[3:4] ; if(do.log) yU <- 10^yU
    if(leg.do)
        legend("topleft",#0, if(abs.do) yU[2] else .7*yU[1]+.3*yU[2],
               paste(format(Ns),"-term",sep=''), bty="n",
               col = 1:nN, lty=1:nN, ncol = ceiling(nN/3))
}

p.g(1)
p.g(.1)
p.g(.001)


{ op <- mult.fig(2)$old.par
  p.g(.1)
  p.g(.1, abs = TRUE, leg.do = FALSE)
  par(op) }

## How bad is "direct", close to 0 : (looks worse than it is ?)
##               ## rrange
p.g(.001, Ns = 0)#  +- 3e-5
p.g(.004, Ns = 0)#  +- 4e-7
p.g(.020, Ns = 0)#  +- 4e-9
p.g(.040, Ns = 0)#  +- 4e-10
p.g(.1,   Ns = 0)#  +- 2e-11
p.g(.2,   Ns = 0)#  +- 3e-12


###
epsC <- .Machine$double.eps

###===== Bound for simple two term (k = 0,1: linear) formula: ===============
(5/2*epsC)  ^(1/2) ## 2.356 e-8
(5/2*epsC/2)^(1/2) ## 1.666 e-8
## slightly inside:
y <- 2.2 * 10^seq(-20,-8, length=1001)

## Absolute error:
summary(ae <- gnt(y,10) - gnt(y,1))
##-       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
##- -1.110e-16  0.000e+00  0.000e+00 -2.118e-17  0.000e+00  1.110e-16

## relative error:
summary(re <- 1 - gnt(y,1) / gnt(y,10))
##-       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
##- -2.220e-16  0.000e+00  0.000e+00 -4.237e-17  0.000e+00  2.220e-16
are <- abs(re)

plot(y,are, log = 'x')
plot(y,are, log = 'x',  xlim=c(1e-10,max(y)))
plot(y,are, log = 'xy', xlim=c(1e-10,max(y)))
plot(y,are, log = '',   xlim=c(1e-10,max(y)))
##--> I should rather use epsC / 2 (makes sense from what we know!)

###===== Bound for three term (k = 0,1,2) formula: ====================
(7/2*epsC)  ^(1/3)## 9.19 e-6
(7/2*epsC/2)^(1/3)## 7.30 e-6

## slightly inside:
y <- 9* 10^seq(-20,-6, length=1001)

## Absolute error:
range(ae <- gnt(y,10) - gnt(y,2))## 0  1.11e-16

## relative error:
range(re <- 1 - gnt(y,2) / gnt(y,10))## 0  2.22e-16

plot(y,re, log = 'x')
plot(y,re, log = 'x', xlim=c(1e-6,max(y)))
plot(y,re, log = 'xy', xlim=c(1e-6,max(y)))
plot(y,re, log = '',   xlim=c(1e-6,max(y)))
###--> I should rather use epsC / 2 (makes sense from what we know!)
### or even a bit less ..:  |y| < 1e-5



###===== Bound for four term (k = 0,1,2,3) formula: ====================
(14/3*epsC)  ^(1/4) ## 0.000179
(14/3*epsC/2)^(1/4) ## 0.000151

## slightly _out_side:
y <- 1.8* 10^seq(-10,-4, length=1001)
## Absolute error:
range(ae <- gnt(y,10) - gnt(y,3))    # 0.  1.11e-16
## relative error (very similar)
range(re <- 1 - gnt(y,3) / gnt(y,10))# 0.  2.22e-16


###===== Bound for five term (k = 0,1,2,3,4) formula: ====================
(6*epsC  )^(1/5) ## 0.00106
(6*epsC/2)^(1/5) ## 0.000922

## slightly inside:
y <- 1* 10^seq(-10,-3, length=1001)
## Absolute error:
range(ae <- gnt(y,10) - gnt(y,4))    # 0  1.11e-16
## relative error (very similar)
range(re <- 1 - gnt(y,4) / gnt(y,10))# 0  2.22e-16


###===== Bound for 6 term (k = 0,1,2,3,4,5) formula: ====================
(15/2*epsC  )^(1/6) ## 0.00344

## slightly inside:
y <- 3.4* 10^seq(-10,-3, length=1001)
## Absolute error:
range(ae <- gnt(y,10) - gnt(y,5))    # 0  1.11e-16
## relative error (very similar)
range(re <- 1 - gnt(y,5) / gnt(y,10))# 0  2.22e-16


###===== Bound for seven term (k = 0,1,2,3,4,5,6) formula: ====================
(55/6*epsC  )^(1/7) ## 0.00797
## slightly outside:
y <- 8* 10^seq(-10,-3, length=1001)
## Absolute error:
range(ae <- gnt(y,10) - gnt(y,6))    # 0  2.22-16
## relative error (very similar)
range(re <- 1 - gnt(y,6) / gnt(y,10))# 0  4.44e-16

##==> MM{2014}: I think I settled for  g2() in 2004,
##    -------- which goes up to 4-term ==== because it is good enough.



###===== former ../tests/wienergerm-accuracy.R ==========================================
### Orig. = R script /u/maechler/R/MM/NUMERICS/dpq-functions/wienergerm-accuracy.R

if(!dev.interactive(orNone=TRUE)) pdf("wienergerm-accuracy.pdf")

nuS <- c(1:3, 10, 20, 30,50, 100, 200, 1000, 5000, 1e4)
ncpS <- c(0, .1, 1, 3, 10, 30, 100, 300, 1000, 10000, 1e5)
facs <- c(outer(c(1,2,5), 10^c(-1:3)))
r <- array(NA, dim = c(length(facs), 3, length(nuS), length(ncpS)),
           dimnames= list(x=NULL, c("p.","pW","pP"),
                          df = formatC(nuS), ncp= formatC(ncpS)))
for(df in nuS) {
    for(ncp in ncpS) {
        m <- df + ncp
        s <- sqrt(2*(df + 2*ncp))
        x <- m + s * facs
        r[,, formatC(df), formatC(ncp)] <-
            cbind(pchisq     (x, df=df, ncp=ncp),
                  pchisqW    (x, df=df, ncp=ncp),
                  pnchisqPearson(x, df=df, ncp=ncp))
    }
    cat("done df=",formatC(df),"\n")
}

p.pchisq <- function(x, df, ncp, log = "x",
                     relErr = FALSE, diff = FALSE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 18 Mar 2004, 17:04

    Px <- pchisq(x, df=df, ncp=ncp)
    pW  <- pchisqW(x, df=df, ncp=ncp) # "Wiener2"
    pW1 <- pchisqW(x, df=df, ncp=ncp, variant = 'f')
    pP1 <- pnchisqPearson(x, df=df, ncp=ncp)
    pP2 <- pnchisqPatnaik(x, df=df, ncp=ncp)
    Pmat <- cbind(pW,pW1,pP1,pP2)

    tit <- paste("pchisq[Appr](x, df=",formatC(df),", ncp=",formatC(ncp),")")
    if(diff)
        tit <- paste(tit," -  pchisq[R](.)")
    else if( relErr)
        tit <- paste("rel.Error of", tit," to pchisq[R](.)")
    log.y <- length(grep("y", log))
    col <- 2:5
    lty <- 1:4
    noR <- diff || relErr # subtract from R's pchisq() -- TODO: rather use Rmpfr-ified pnchisq()!
    if(noR) {
         Pmat <- Pmat - Px
         if(relErr)
             Pmat <- abs(Pmat / Px)
         else if(log.y) {
             Pmat <- abs(Pmat)
             if(log.y) tit <- paste("|", tit, "|")
         }
    }
    else {
        Pmat <- cbind(Px, Pmat)
        col <- c(1,col)
        lty <- c(1,lty)
    }
    if(log.y) op <- options(warn=-1)# no "log warning"
    matplot (x, Pmat, log=log, type='l', main = tit,
             ylab = '', ## ylim = rrange(Pmat, range = 0.1),
             lty = lty, col = col)
    if(log.y) op <- options(op)
    u <- par("usr")
    log.x <- length(grep("x", log))
    if(log.x) u[1:2] <- 10^u[1:2]
    if(log.y) u[3:4] <- 10^u[3:4]
    legend(u[1],u[4],
           c(if(!noR)"R", "Wiener2", "Wiener1", "Pearson","Patnaik"),
           lty = lty, col = col)
}

### all these for (0.1, 0.1):
x <- 2^seq(-10,3, length=501)
p.pchisq(x, df=0.1, ncp=0.1)
p.pchisq(2^seq(-12,6, length=2001), df=0.1, ncp=0.1)
p.pchisq(2^seq(-1, 6, length=2001), df=0.1, ncp=0.1)
p.pchisq(2^seq(-1, 6, length=2001), df=0.1, ncp=0.1, diff = TRUE)

p.pchisq(2^seq(-1, 10, length=2001), df=0.1, ncp=0.1, diff = TRUE)
p.pchisq(2^seq(-1, 10, length=2001), df=0.1, ncp=0.1, diff = TRUE,log="xy")
p.pchisq(2^seq(-5,  8, length=2001), df=0.1, ncp=0.1, diff = TRUE,log="xy")
p.pchisq(2^seq(-5,  8, length=2001), df=0.1, ncp=0.1, relE = TRUE,log="xy")

## (1, 1)
p.pchisq(2^seq(-14,8, length=2001), df=1, ncp= 1)
p.pchisq(2^seq(-14,8, length=2001), df=1, ncp= 1,diff=TRUE)
p.pchisq(2^seq(-20,9, length=2001), df=1, ncp= 1,diff=TRUE,log="xy")
p.pchisq(2^seq(-20,9, length=2001), df=1, ncp= 1,relE=TRUE,log="xy")

## (1, 10)
p.pchisq(2^seq(-14,8, length=2001), df=1, ncp= 10,diff=TRUE,log="xy")
## (10, 1)
p.pchisq(2^seq(-14,8, length=2001), df=10, ncp= 1,diff=TRUE,log="xy")

## (10, 10)
p.pchisq(2^seq(-14,8, length=2001), df=10, ncp= 10,diff=TRUE,log="xy")
p.pchisq(2^seq(-1,10, length=2001), df=10, ncp= 10,diff=TRUE,log="xy")

## (10, 100)
p.pchisq(2^seq(-14,8, length=2001), df=10, ncp= 100,diff=TRUE,log="xy")
p.pchisq(2^seq(-1,10, length=2001), df=10, ncp= 100,diff=TRUE,log="xy")
## "good" news: Pearson is "nice" when Wiener fails

## (10, 1e4)
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, log='')
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, diff=TRUE)
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, diff=TRUE, log='y')
p.pchisq(seq(9000, 11000, length=2001), df=10, ncp= 1e4, relE=TRUE, log='y')

## (1000, 1)
p.pchisq(seq(800, 1200, length=2001), df=1000, ncp= 1, log='')

p.pchisq(seq(800, 1200, length=2001), df=1000, ncp= 1, diff=TRUE, log='y')
p.pchisq(seq(800, 1200, length=2001), df=1000, ncp= 1, relE=TRUE, log='y')
## interesting: Pearson really better than Wiener in a whole region
## x in [900, 1100] ...
