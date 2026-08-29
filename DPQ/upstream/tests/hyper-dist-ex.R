#### Used to be part of /u/maechler/R/MM/NUMERICS/hyper-dist.R -- till Jan.24 2020

library(DPQ)

#### First version committed is the code "as in Apr 1999":
####  file | 11782 | Apr 22 1999 | hyper-dist.R

options(rErr.eps = 1e-30)
rErr <- function(approx, true, eps = getOption("rErr.eps", 1e-30))
{
    ifelse(Mod(true) >= eps,
	   1 - approx / true, # relative error
	   true - approx)     # absolute error (e.g. when true=0)
}

if(!dev.interactive(orNone=TRUE)) pdf("hyper-dist-ex.pdf")


### ----------- -----------

rErr.phypDiv <- function(q, m,n,k) {
    ph  <- phyper(q, m=m, n=n, k=k)
    cbind(
      rERR.Ibeta= rErr(phyperIbeta     (q, m=m, n=n, k=k) , ph),
      rERR.as151= rErr(phyperApprAS152 (q, m=m, n=n, k=k) , ph),
      rERR.1mol = rErr(phyper1molenaar (q, m=m, n=n, k=k) , ph),
      rERR.2mol = rErr(phyper2molenaar (q, m=m, n=n, k=k) , ph),
      rERR.Peizer=rErr(phyperPeizer    (q, m=m, n=n, k=k) , ph))
}

## Plotting of the rel.error:
p.rErr.phypDiv <- function(m,n,k, q = .suppHyper(m,n,k), abslog = FALSE, logx= "", type = "b")
{
    rE <- rErr.phypDiv(m=m, n=n, k=k, q=q)
    nc <- ncol(rE)
    cc <- adjustcolor(1:nc, 0.5)
    ppch <- as.character(1:nc)
    cl <- sys.call()
    cl[[1]] <- quote(phyper)
    matplot(q, if(abslog) abs(rE) else rE, type = type, lwd=2, col=cc, pch=ppch,
            ylab = if(abslog) quote(abs(rE)) else quote(rE),
            log = paste0(logx, if(abslog) "y" else ""),
            main = sprintf("%s", deparse1(cl)))#
    mtext(paste(if(abslog) "|rel.Err|" else "rel.Err.",
                "of diverse phyper() pbinom* approximations"))
    if(!abslog) abline(h=0, col=adjustcolor("gray", 0.5), lty=2)
    legend("topright", legend = paste(1:ncol(rE), colnames(rE), sep=": "),
           lwd=2, col=cc, lty=1:4, pch=ppch, bty="n")
    invisible(rE)
}


k <- c(10*(1:9),100*(1:9),1000*(1:9))
k1 <- k
options(digits=4)

 rErr12 <- p.rErr.phypDiv(q=k, 2*k,2*k,2*k)
           p.rErr.phypDiv(q=k, 2*k,2*k,2*k, abslog=TRUE)

k <- c(10*(7:9),100*(1:9),1000*(1:9), 1e4*(1:9))
k2 <- k

p.rErr.phypDiv(q=k, 2*k,2*k,2*k, abslog=TRUE, logx="x")

k <-  c(10*(1:9),100*(1:3)) # small k only (afterwards --> all phyper*() == 1 !)(revert)
## Here, the Ibeta  fails (NaN) ; *also* Molenaar's (???? contradiction to below)
(rErr1215 <- p.rErr.phypDiv(q=k, 1.2*k, 2*k, 1.5*k))
             p.rErr.phypDiv(q=k, 1.2*k, 2*k, 1.5*k, abslog=TRUE, logx="x")

(rErr1215 <- p.rErr.phypDiv(q=k, 1.2*k, 2*k, 1.5*k))
             p.rErr.phypDiv(q=k, 1.2*k, 2*k, 1.5*k, abslog=TRUE, logx="x")


k <- c(10*(7:9), t(outer(10^(2:5), 1:9)))
x <- round(.8*k); ph.k2k <- phyper(x,1.6*k, 2*k,1.8*k)
(rErr1618 <- p.rErr.phypDiv(q=x, 1.6*k, 2*k, 1.8*k, logx="x")) ## AS151 is really unusable here !
             p.rErr.phypDiv(q=x, 1.6*k, 2*k, 1.8*k, abslog=TRUE, logx="x")
## interestingly, Peizer suffers from some erratic behavior for large q
## --> the L has  4  log(.) terms, each of which has argument ~= 1 <--> could use log1p(.) or ???

k <- k1 # the old set
par(mfrow=c(2,1))
for(Reps in c(0,1)) {
    options(rErr.eps = Reps) ## 0: always RELATIVE error; 1: always ABSOLUTE
    x <- round(.6*k); ph.k2k <- phyper(x,1.6*k, 2*k,1.8*k)
    ##         ^^^
    print(ph.mat <-
          cbind(k=k, phyper= ph.k2k,
           rERR.Ibeta= rErr(phyperIbeta     (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.as151= rErr(phyperApprAS152 (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.1mol = rErr(phyper1molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.2mol = rErr(phyper2molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.Peiz = rErr(phyperPeizer    (x,1.6*k,2*k,1.8*k) , ph.k2k))
          )
    ## The two molenaar's  ``break down''; Peizer remains decent!
    Err.phyp <- abs(ph.mat[,-(1:2)])
    matplot(k, Err.phyp,
            ylim= if((ee <- .Options$rErr.eps) > 0) c(1e-200,1),
            main="|relE.{phyper(x = 0.6k, 1.6k, 2k, 1.8k)}|",
            type='o', log='xy')
}


hp1 <- list(n = 10, n1 = 12, n2 = 2)

for(h.nr in 1:1) {
  attach(hplist <- get(paste("hp",h.nr, sep='')) ) # defining  n, n1, n2
  cat("\n");  print(unlist(hplist))
  N <- n1 + n2
  x <- 0:max(0,min(n1, n2 - n))
  d <- dhyper(x, n1, n2, n)
  names(d) <- as.character(x)
  print(d)
}


### Binomial Approximation(s) =================================================

ph.m <- cbind(ph  = phyper             (0:5, 5,15, 7),
              pM1 = phyperBinMolenaar.1(0:5, 5,15, 7),
              pM2 = phyperBinMolenaar.2(0:5, 5,15, 7),
              pM3 = phyperBinMolenaar.3(0:5, 5,15, 7),
              pM4 = phyperBinMolenaar.4(0:5, 5,15, 7))

cc <- adjustcolor(1:(1+ncol(ph.m)), 0.5); ppch <- as.character(0:4)
matplot(0:5, ph.m, type = "b", lwd=2, col=cc, pch=ppch,
        main = "all 4 phyper() binomial approximations via Molenaar's")
legend("right", legend = c("phyper", paste("ph.appr.Mol", 1:4)),
       lwd=2, col=cc, lty=1:5, pch=ppch, bty="n")


## relative errors :
ph.err <- 1 - ph.m[,-1] / ph.m[,"ph"]

matplot(0:5, ph.err, type = "b", lwd=2, col=cc, pch=ppch,
        main = "rel.Err. of the 4 phyper() binom.Molenaar's approximations")
legend("topright", legend = paste("ph.appr.Mol", 1:4),
       lwd=2, col=cc[-1], lty=2:5, pch=ppch[-1], bty="n")
## quite interesting: if we take the *mean* of approx. 1 and 2 -- should become good?

## have rErr() above {but does not work with *matrix* ?}

rErr.Mol.bin <- function(m,n,k, q = .suppHyper(m,n,k), lower.tail=TRUE, log.p=FALSE) {
    ph  <- phyper      (q, m=m, n=n, k=k,      lower.tail=lower.tail, log.p=log.p)
    phM <- phyperAllBinM(m=m, n=n, k=k, q=q, lower.tail=lower.tail, log.p=log.p)
    apply(phM, 2, rErr, true = ph)
}
## Plotting of the rel.error:
p.rErr.Mol.bin <- function(m,n,k, q = .suppHyper(m,n,k), abslog = FALSE,
                           lower.tail=TRUE, log.p=FALSE)
{
    rE <- rErr.Mol.bin(m=m, n=n, k=k, q=q, lower.tail=lower.tail, log.p=log.p)
    cc <- adjustcolor(1:4, 0.5)
    ppch <- as.character(1:4)
    matplot(q, if(abslog) abs(rE) else rE, type = "b", lwd=2, col=cc, pch=ppch,
            log = if(abslog) "y" else "",
            main = sprintf("phyper(*, m = %d, n = %d, k = %d)", m,n,k))
    mtext(paste(if(abslog) "|rel.Err|" else "rel.Err.",
                "of the 4 phyper() binom.Molenaar's approximations"))
    if(!abslog) abline(h=0, col=adjustcolor("gray", 0.5), lty=2)
    legend("topright", legend = paste("phyp.Mol", 1:4),
           lwd=2, col=cc, lty=1:4, pch=ppch, bty="n")
    invisible(rE)
}

p.rErr.Mol.bin (5,15, 7) # same plot, as the "manual" one above
p.rErr.Mol.bin (70,100, 7)# here, 1 and 4 are good
p.rErr.Mol.bin (70,100, 20)#  4 (and 1) ((maybe their *mean* !)
p.rErr.Mol.bin (70,100, 20, abslog=TRUE)#  4 (and 1)
p.rErr.Mol.bin (70,100, 50)#  4
p.rErr.Mol.bin (70,100, 50, abslog=TRUE) -> rE.71c50
(p.rErr.Mol.bin (70,100, 70,abslog=TRUE) -> rE.71c70)#  1 & 2; but *none* is good for small q
par(new=TRUE)
plot(0:70, phyper(0:70, 70,100, 70), type = "l", col="blue", ann=FALSE,axes=FALSE)

## rel error on log-scale ... really nonsense
(p.rErr.Mol.bin (70,100, 70, log.p=TRUE))#  1 & 2; but *none* is good for small q
## however, *again* the *mean* of  1|2  with 3|4 seems fine


## Total Variation error -- what Kuensch(1998) used for binom() -- hyper() approx.
TVerr <- function(approx, true) {
    d <- true - approx
    max(sum(  d[d > 0]),
        sum(- d[d < 0]))
}
supErr <- function(approx, true) max(abs(true - approx))

TVerr.bin <- function(m,n,k, q = .suppHyper(m,n,k), lower.tail=TRUE, log.p=FALSE) {
    ph  <- phyper     (q, m=m, n=n, k=k,    lower.tail=lower.tail, log.p=log.p)
    phM <- phyperAllBin(m=m, n=n, k=k, q=q, lower.tail=lower.tail, log.p=log.p)
    apply(phM, 2, TVerr, true = ph)
}
supErr.binM <- function(m,n,k, q = .suppHyper(m,n,k), lower.tail=TRUE, log.p=FALSE) {
    ph  <- phyper      (q, m=m, n=n, k=k,    lower.tail=lower.tail, log.p=log.p)
    phM <- phyperAllBinM(m=m, n=n, k=k, q=q, lower.tail=lower.tail, log.p=log.p)
    apply(phM, 2, supErr, true = ph)
}

(ph.TV  <- apply(ph.m[,-1], 2, TVerr, true = ph.m[,1]))
## in all three cases,  Molenaar is clearly better than "simple"
TVerr.bin(5,15,7)
TVerr.bin(1000,2000,7) # 1 & 4 are very accurate
TVerr.bin(1000,2000, 1500)# not very good

TVerr.bin(100000, 200000, 150000)
## but TV error really *grows* with the support of the distribution
supErr.binM(100000, 200000, 150000)
supErr.binM(1000000, 2000000, 150000)
supErr.binM(1e6, 1e9, 150000)
round(-log10(supErr.binM(1e6, 1e9, 1e7)), 2) #=> correct #{digits}
##  pM1  pM2  pM3  pM4
## 5.98 7.98 7.99 5.99  -- m = 1e6 < k = 1e7 ==> approx. with size = m are better
system.time(sE <- supErr.binM(1e7, 1e9, 1e7)) ## 20 sec. elapsed! ==> very slow !!
round(-log10(sE), 2)
##  pM1  pM2  pM3  pM4
## 6.00 6.00 6.01 6.01 --- m = k ==> all are the same

system.time(sE <- supErr.binM(1e7, 1e9, 100))# fast
round(-log10(sE), 2)
##   pM1   pM2   pM3   pM4
## 15.35  5.22  5.21 14.65  k << m ==> 1 & 4 are much better {why '1' better than '4' ?}




### Here, we currently only explore normal approx:

k <- c(10*(1:9),100*(1:9),1000*(1:9))

options(digits=4)
ph.k2k <- phyper(k,2*k,2*k,2*k) # rel.err ~ 10^-7
cbind(k=k, phyper= ph.k2k,
      rERR.Ibeta= rErr(phyperIbeta     (k,2*k,2*k,2*k) , ph.k2k),
      rERR.as151= rErr(phyperApprAS152 (k,2*k,2*k,2*k) , ph.k2k),
      rERR.1mol = rErr(phyper1molenaar (k,2*k,2*k,2*k) , ph.k2k),
      rERR.2mol = rErr(phyper2molenaar (k,2*k,2*k,2*k) , ph.k2k),
      rERR.Peizer=rErr(phyperPeizer    (k,2*k,2*k,2*k) , ph.k2k))

## Here, the Ibeta  fails (NaN) ; moleaar's are both very good :
ph.k2k <- phyper(k, 1.2*k, 2*k,1.5*k)
cbind(k=k, phyper= ph.k2k,
      rERR.Ibeta= rErr(phyperIbeta     (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.as151= rErr(phyperApprAS152 (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.1mol = rErr(phyper1molenaar (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.2mol = rErr(phyper2molenaar (k,1.2*k,2*k,1.5*k) , ph.k2k),
      rERR.Peizer=rErr(phyperPeizer    (k,1.2*k,2*k,1.5*k) , ph.k2k))

x <- round(.8*k); ph.k2k <- phyper(x,1.6*k, 2*k,1.8*k)
(ph.mat <-
cbind(k=k, phyper= ph.k2k,
      rERR.Ibeta= rErr(phyperIbeta     (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.as151= rErr(phyperApprAS152 (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.1mol = rErr(phyper1molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.2mol = rErr(phyper2molenaar (x,1.6*k,2*k,1.8*k) , ph.k2k),
      rERR.Peiz = rErr(phyperPeizer    (x,1.6*k,2*k,1.8*k) , ph.k2k))
)
matplot(k, abs(ph.mat[,-(1:2)]), type='o', log='xy')

op <- if(require("sfsmisc")) mult.fig(2)$old.par else par(mfrow=c(2,1))
for(Reps in c(0,1)) {
    options(rErr.eps = Reps) ## 0: always RELATIVE error; 1: always ABSOLUTE
    tit <- paste("phyper() approximations:",
                 if(Reps == 0) "relative" else "absolute", "error")
    cat(tit,":\n")
    x <- round(.6*k); ph.k2k <- phyper(x,1.6*k, 2*k,1.8*k)
    print(ph.mat <-
          cbind(k=k, phyper= ph.k2k,
           rERR.Ibeta= rErr(phyperIbeta        (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.as151= rErr(phyperApprAS152    (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.1mol = rErr(phyper1molenaar    (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.2mol = rErr(phyper2molenaar    (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rERR.Peiz = rErr(phyperPeizer       (x,1.6*k,2*k,1.8*k) , ph.k2k),
           rErr.binMol=rErr(phyperBinMolenaar.1(x,1.6*k,2*k,1.8*k) , ph.k2k))
          )
    ## The two molenaar's  ``break down''; Peizer remains decent!
    Err.phyp <- abs(ph.mat[,-(1:2)])
    matplot(k, Err.phyp,
            ylim= if((ee <- .Options$rErr.eps) > 0) c(1e-200,1),
            main = tit, type='o', log='xy',
            lty=1:6, col=1:6, pch=as.character(1:6))
    if(Reps == 1)
        legend("bottomleft",
               c("Ibeta", "AS 152", "1_Molenaar", "2_Molenaar", "Peizer",
                 "binom_Mol."),
               bty = "n", lty=1:6, col=1:6, pch=as.character(1:6))
}
par(op)

hp1 <- list(n = 10, n1 = 12, n2 = 2)

for(h.nr in 1:1) {
  attach(hplist <- get(paste("hp",h.nr, sep='')) ) # defining  n, n1, n2
  cat("\n");  print(unlist(hplist))
  N <- n1 + n2
  x <- 0:max(0,min(n1, n2 - n))
  d <- dhyper(x, n1, n2, n)
  names(d) <- as.character(x)
  print(d)
}

### ----------- ----------- lfastchoose() etc -----------------------------


sapply(0:10, function(n) lfastchoose(n, c(0,n)))

## System  lchoose  gives non-sense:
sapply(0:10, function(n)    lchoose(n, c(0,n)))
## This is ok:
sapply(0:10, function(n) f05lchoose(n, c(0,n)))


###---------- some  gamma testing: -------------

pi - gamma(1/2)^2
for(n in 1:20) cat(n,":",formatC(log(prod(1:n)) - lgamma(n+1)),"\n")
for(n in 1:20) cat(n,":",formatC(prod(1:n) - gamma(n+1)),"\n")

## in  math/gamma.c :
p1 <- c(0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2)

### Bernoulli numbers  Bern() :

n <- 0:12
Bn <- n; for(i in n) Bn[i+1] <- Bern(i); cbind(n, Bn)
system.time(Bern(30))
system.time(Bern(30))
if(FALSE){ ## rm(.Bernoulli)
system.time(Bern(30))
system.time(Bern(80))
}

### lgammaAsymp() --- Asymptotic log gamma function :

lgammaAsymp(3,0)
for(n in 0:10) print(log(2) - lgammaAsymp(3,n))
for(n in 0:12) print((log(2*3*4*5) - lgammaAsymp(5+1,n)) / .Machine$double.eps)
for(n in 0:12) print((log(720)     - lgammaAsymp(6+1,n)) / .Machine$double.eps)
for(n in 0:12) print((log(720)     - lgamma     (6+1)  ) / .Machine$double.eps)
for(n in 0:12) print((log(5040)    - lgammaAsymp(7+1,n)) / .Machine$double.eps)


## look ok:
curve(lgamma(x),-70,10, n= 1001) # ok
curve(lgamma(x),-70,10, n=20001) #  slowness of x11 !!

curve(abs(gamma(x)),-70,10, n=  201, log = 'y')
curve(abs(gamma(x)),-70,70, n=10001, log = 'y')
par(new=T)
curve(lgamma(x),    -70,70, n= 1001, col = 'red')

## .., we should NOT  use log - scale with negative values... [graph ok, now]
curve(gamma(x),-     40,10, n=  2001, log = 'y')
curve(abs(gamma(x)),-40,10, n=  2001, log = 'y')
##- Warning: NAs produced in function "gamma"

x <- seq(-40,10, length = 201)
plot(x, gamma(x), log = 'y', type='h')


###------------------- early dhyper() problem ---------------------------------

##- Date: 30 Apr 97 11:05:00 EST
##- From: "Ennapadam Venkatraman" <VENKAT@biosta.mskcc.org>
##- Subject: R-beta: dhyper bug??
##- To: "r-testers" <r-testers@stat.math.ethz.ch>
##- Sender: owner-r-help@stat.math.ethz.ch

## I get the following incorrect answer  ---- FIXED

dhyper(34,410,312,49)
#   [1] 2.244973e-118
##- when the coorect value is
##-
dhyper(34,410,312,49)#          (from S-Plus)
##-    [1] 0.0218911

##================== R is now correct !
stopifnot(all.equal(dhyper(34,410,312,49),
                    0.021891095726, tol=1e-11))
##E.S. Venkatraman (venkat@biosta.mskcc.org)


### phyperR() === R version of pre-2004 C version in <R>/src/nmath/phyper.c :

## This takes long, currently (1999??)  {18 sec, on sophie [Ultra 1]}  -- now all fast!
k <- (1:100)*1000
system.time(phyper.k <- phyper(k, 2*k, 2*k, 2*k))

k <- 1000; phyperR(k, 2*k, 2*k, 2*k) - phyper.k[1]

## Debug: ------------
if(FALSE) # for now, when I source this
if(interactive())
  debug(phyperR)
phyperR(k, 2*k, 2*k, 2*k)
## before while():
xb <- 2000
xr <- 0
ltrm <- -2768.21584356709
NR <- 2000
NB <- 0# new value

