#### History(RCS): /u/maechler/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/qbeta-dist.R,v

#### Testing  qbeta(.),   pbeta(.),   qt(.), .....
#### ----------------

source(system.file(package="DPQ", "test-tools.R",
                   mustWork=TRUE))# ../inst/test-tools.R
## => showProc.time(), .. {from Matrix},
##    list_() , loadList() ,  readRDS_() , save2RDS()

(doExtras <- DPQ:::doExtras())
## save directory (to read from):
(sdir <- system.file("safe", package="DPQ"))

if(!dev.interactive(orNone=TRUE)) pdf("qbeta-dist.pdf")
.O.P. <- par(no.readonly=TRUE)
op <- options(nwarnings = 1e5)

fcat <- function(..., f.dig= 4, f.wid = f.dig +5, f.flag = ' ', nl = TRUE,
                 file = "", sep = " ",
                 fill = FALSE, labels = NULL, append = FALSE)
{
  ## Purpose: Formatted CAT -- for printing 'tables'
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 12 May 97
  l <- unlist(lapply(list(...), formatC,
                     wid= f.wid, digits= f.dig, flag= f.flag))
  cat(l, if(nl)"\n", file=file, sep=sep, fill=fill,labels=labels, append=append)
}

format01prec <- function(x, digits = getOption("digits"), width = digits + 2,
                        eps = 1e-6, ...,
                        FUN = function(x,...) formatC(x, flag='-',...))
{
  ## Purpose: format numbers in [0,1] with "precise" result,
  ##          using "1-.." if necessary.
  ## -------------------------------------------------------------------------
  ## Arguments: x:      numbers in [0,1]; (still works if not)
  ##            digits, width: number of digits, width to use with 'FUN'
  ##            eps:    Use '1-' iff  x in  (1-eps, 1] -- 1e-6 is OPTIMAL
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 May 97, 18:07
  if(as.integer(digits) < 4) stop('digits must be >= 4')
  if(eps < 0 || eps > .1) stop('eps must be in [0, .1]')
  i.swap <- is.na(x) | (1-eps < x  &  x <= 1) #-- Use "1- ." if <~ 1,  normal 'FUN' otherwise
  r <- character(length(x))
  if(any(i.swap))
    r[i.swap] <-
      paste("1-", FUN(1-x[i.swap], digits=digits - 5, width=width-2, ...),
            sep='')# -5: '1-' + 4 for exponent -1 for '0' (in other case)
  if(any(!i.swap))
    r[!i.swap] <- FUN(x[!i.swap], digits=digits, width=width,...)
  attributes(r) <- attributes(x)
  r
}


###---------------------- short pre-millennium examples: ----------------

## = ./beta-test-0.R, Jul 9, 1997 (!)
showProc.time()

qbeta(1e-3, 1,10^(1:10))
qbeta(1e-5, 1,10^(1:10))
qbeta(3e-2, 1,10^(1:10))
qbeta(2e-2, 1,10^(1:10))

##-- quite problematic ones :
cbind(qbeta(2e-2, 1,    10^(0:32)))
cbind(qbeta(.2,   1,    10^(0:32)))
cbind(qbeta(.2,  .1,    10^(0:32)))
cbind(qbeta(.1,   1e-1, 10^(0:32)))# "good"
cbind(qbeta(.1,   1e-2, 10^(0:32)))# look good
cbind(qbeta(.1,   1e-2, 10^(0:64)))# look good

plot(qbeta(1/8, 2^-7,  10^(0:64)), log="y")# look good
plot(qbeta(1/8, 2^-8,  10^(0:64)), log="y")# look good (down to 10^-298)
 all(qbeta(1/8, 2^-9,  10^(0:64)) == 0)# TRUE, but ...
 all(qbeta(1/8, 2^-10, 10^(0:64)) == 0)# TRUE, but should be positive (in principle)

stopifnot(qbeta(.1, 1e-3, 10^(0:32)) == 0) # and we cannot do better:

curve(pbeta(x, 1e-3,  1  ), 1e-310, log="x", col=2, xaxt="n")# to give nicer axis:
sfsmisc::eaxis(1); axis(1, at=1); abline(h=1,v=1, lty=2)
curve(pbeta(x, 1e-3, 1000), add=TRUE, col=3)
curve(pbeta(x, 1e-3, 1e10), add=TRUE, col=4)
curve(pbeta(x, 1e-3, 1e20), add=TRUE, col=5)
curve(pbeta(x, 1e-3, 1e50), add=TRUE, col=6)
curve(pbeta(x, 1e-3,1e100), add=TRUE, col=5)
curve(pbeta(x, 1e-3,1e200), add=TRUE, col=4)## ==> Warnings already
curve(pbeta(x, 1e-3,1e300), add=TRUE, col=2)## ==> WARNINGS:
## 1: In pbeta(x, ...) : bgrat(a=1e+300, b=0.001, x=1) *no* convergence: NOTIFY R-core!
summary(warnings()) #                                  ----------------  ^^^^^^^^^^^^^^

## = ./beta-qbeta-mintst.R, July 31, 1997 (!)
1 - qbeta(.05, 5e10, 1/2) #-- infinite loop in  0.49; 0.50-a1
##-- 3.8416e-11 with MM's qbeta(.)
qbeta(.95, 1/2, 5e10) # the "same" (with swap_tail)
showProc.time()

###-----------------------  qt(p,  df --> Inf)  -->  qnorm(p) -----------------

p.<- 5* 10^-c(6:8,10,15,20,100,300)#- well ...
p <- p. #- try things below with this 'EXTREME' p
p <- c(.975, .995, .9995, 1 - 5e-6)
df1 <- c(1:5,10^c(1:21,25,10*(3:10)))
df0 <- df1[df1< 10^30] ## df0  works for  R 0.49 -patch2-
df0f<- unique(sort(c((1:30)/10, df0)))


## system.time() gives 0.0 nowadays ...
sys.time <- function(exp) system.time(for(i in 1:5000) exp)
dig <- 18
dig <- 12
dig <-  8

cat("    df : CPU[s] | p=", formatC(  p ,  form='e', dig=5,  wid=dig+6),"\n")
for(df in df0)
  cat(formatC(df,w=6),":", formatC(sys.time(qq <- qt(p, df))[1],form='f'), "|",
      formatC(qq, form='f', dig=dig,wid=dig+6),"\n")
  cat(" _Inf_ :  ---   |", formatC(qnorm(p), form='f', dig=dig,wid=dig+6),"\n")


### Use  F-distribution :   `` F_{1,nu} ==  (t_{nu}) ^2 ''
### Therefore:  qt(p,n) == sqrt(qf(2*p-1, 1,n))  for p >= 1/2
###	        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -- was true for ORIGINAL in R 0.49
##for(df in df0f) { ##-- ~ TRUE for (0.49) qt/qf which worked with qbeta
if(any(p < 1/2)) stop("must have p >= 1/2 for these examples!")

for(df in df0f) {
  cat(formatC(df,w = 6),":", all.equal( sqrt(qf(2*p-1, 1,df)), qt(p,df)),"\n")
}
## R 3.0.3 : all TRUE, but
## 1e+06 : Mean relative difference: 3.240032e-06
## 1e+07 : Mean relative difference: 3.240022e-07
## 1e+08 : Mean relative difference: 3.240021e-08
## ==> Applying pt(*) ==> (R 3.5.1) shows that qt() is *better* than qf(.) equivalence

## It seems, qt(.) is perfect
p <- 1 - 1/64 # exact rational
df <- 2^seq(10,20, by=1/32)
plot(df,              qt(p, df=df)- qnorm(p), log="xy", type="l")
plot(df[-1],    -diff(qt(p, df=df)),          log="xy", type="l")
plot(df[-(1:2)], diff(qt(p, df=df), diff=2),  log="xy", type="l")
all.equal( sqrt(qf(2*p-1, 1,df)), qt(p,df)) ## -> 3.15636e-07 (R 3.0.3, 64bit)
## looking closely:
plot(df, qf(2*p-1, 1,df) - qt(p,df)^2, type="b")
## shows clear jump at df = 400'000
df <- 1000*seq(390, 410, by=1/16)
plot(df, qf(2*p-1, 1,df) - qt(p,df)^2, type="l")
## ~/R/D/r-devel/R/src/nmath/qf.c : it uses qchisq() for df1 <= df2, df2 > 400000
## ==> FIXME: qbeta() is *clearly* better here (than "df = Inf" !)
##     =====  qf() is  *clearly* improvable

showProc.time()

### Really, the problem is  qbeta (.):
##                          ---------
## qt (for x > 1/2)  __used__ to be defined as qt(x,n) below, but not longer !!
## qt(x,n)  := sqrt((1  / qbeta( 2*(1-x), n/2, 1/2) - 1) * n)
##--- with DEBUG, shows already BIG problems (needs *much* CPU time for the large  df's
if(FALSE) # is too large and slow
for(p in 2^-(30:2)) {
  cat("\n\np = ", formatC(p),"\n===========\n")
  cat(sprintf("%6s : %5.3s | %11s %11s | %s\n",
              "df","cpu","qbeta","qt", "all.equal(sqrt((1/qb - 1)*.), qt(.))"))
  for(df in c(1:5,10^c(1:21,25,10*(3:10)))) {
    cpu <- system.time(qb <- qbeta(2 * p, df/2, 1/2))[1]
    qt. <- qt(p,df, lower.tail=FALSE)
    cat(sprintf("%6g : %5.3g | %11g %11g | %s\n",
                df, cpu, qb, qt., all.equal(sqrt(df*(1/qb - 1)), qt.)))
  }
}



library(DPQ)
showProc.time()

##----------- Test the  'qnormAppr' function as used in  qbeta(.) :

p <- (0:256)/256
matplot(p, cbind(qnorm(p), qnormAppr(p)), type="l", ylab="")
legend("topleft", c("qnorm(p)", "qnormAppr(p)"), col=1:2, lty=1:2, bty="n")
abline(h=0,v=(0:2)/2, lty=2, col="gray60")
## absolute error:
curve(qnorm(x) - qnormAppr(x), col=2, lwd=2, n=1001, 0, 1)
abline(h=0,v=(0:2)/2, lty=2, col="gray60")

##==> Really make use of symmetry --> use qnormUappr():
matplot(p, cbind(qnorm(p), qnormUappr(p)), type="l", ylab="")
abline(h=0,v=(0:2)/2, lty=2, col="gray60")

## absolute error:
c.1 <- rgb(1,0,0); c.2 <- adjustcolor("black",.5)
curve(qnorm(x) - qnormUappr(x), col=c.1, lwd=2, n=1001, 0, 1)
curve(qnorm(x) - qnormAppr (x), col=c.2, lwd=4, n=1001,
      add=TRUE, lty=3)
legend(.1, .0025, paste("qnorm(.) -",c("qnormUappr(.)", " qnormAppr (.)")),
       col=c(c.1,c.2), lwd=c(2,3), lty=c(1,3), bty="n")
abline(h=0,v=(0:2)/2, lty=2, col="gray60")
showProc.time()


## From R's source  <R>/tests/d-p-q-r-tests.R :
options(rErr.eps = 1e-30) # {*not* to be reset via options(op)}
rErr <- function(approx, true, eps = .Options$rErr.eps)
{
    if(is.null(eps)) { eps <- 1e-30; options(rErr.eps = eps) }
    ifelse(Mod(true) >= eps,
	   1 - approx / true, # relative error
	   true - approx)     # absolute error (e.g. when true=0)
}

## The relative error of this approximation is quite Asymmetric: mainly < 0
x <- c(10^(-9:-4),seq(1e-3, .6, len = 1+2^11))
plot(1-x, rErr(qnormAppr(1-x), qnorm(1-x)), type = 'l', col = 'red',
     main = "Rel. Error of  'qnormAppr(1-x)'")
abline(h = 0, col = 'gray')
## Note: qnorm(.) used to use  AS 111 (26)1977, but has been changed to
## ----  AS 241 (37)1988  which is more accurate {and a bit slower}
##
## Much more sensical in modern R {but looks identical, for 1-x > 1/2 }:
lines(1-x, rErr(qnormUappr(x, lower=FALSE), qnorm(x, lower=FALSE)),
      lwd=3, col = adjustcolor("skyblue", 1/2))

curve(qnorm(x) - qnormUappr(x), col=2, lwd=2, n=1001, .44, 1)
abline(h=0, col=adjustcolor(1,.5)); axis(1, at=.44)
## Warning message:
## In qnormUappr(x) : p[.] < 1/2 & lower tail  is inaccurate
curve(qnorm(x) - qnormAppr(x), lwd=5, n=1001, col=adjustcolor(3,1/3), add=TRUE)
##par(new=TRUE)
## 1/10 * rel.error():
cblu <- adjustcolor("blue4",.5)
curve(rErr(qnormUappr(x), qnorm(x)) / 10, n=1001,
      col = cblu, lwd=3, lty=2, add=TRUE)
## axis(4, col=cblu, col.axis=cblu)
showProc.time()


###-------- Testing  my own   qbeta.R(.)  [[and  qbetaAppr(.) used there]

qbeta  (.05 ,  10/2, 1/2) #-==  0.6682436
qbeta  (.05 , 100/2, 1/2) #-==  0.9621292
qbeta  (.05 ,1000/2, 1/2) #-==  0.996164

###=========== FIXME: More qbetaAppr*() versions
###=============================================
##  and use  lower.tail=TRUE/FALSE   *and*  log.p=FALSE/TRUE
##  and  determine regions where  qbetaAppr*() are very accurate ((so, say 1 Newton step then is perfect!))

qb.names <- paste0("qbeta", c('', '.R', 'Appr', 'Appr.1'))
qf <- lapply(setNames(,qb.names), get)
str(qf) ## for (qn in qb.names) cat(qn,":", deparse(args(get(qn))),"\n\n")

p.set2 <- c(1:5,10,20,50,100,500,1000, 10^c(4:20,25,10*3:10))
p.set2 <- c(1:5,10,20,50,100,500,1000, 10^c(4:18,20,30,100))
## TODO: still too expensive
p.set2 <- c(1e-10,1e-5,1e-3,.1,.3,.5, .8, 1,10, 10^c(6:18,20,30,100))
q.set2 <- c(1e-4, 1/2, 1, 3/2, 2, 5, 1e4)

pn2 <- paste("p=",formatC(p.set2,wid = 1),sep = '')
qn2 <- paste("q=",formatC(q.set2,wid = 1),sep = '')

sfil1 <- file.path(sdir, "tests_qbeta-d-ssR.rds")
if(!doExtras && file.exists(sfil1)) {
  ssR_l <- readRDS_(sfil1)
  str(ssR_l)
  loadList(ssR_l)

} else { ## do run the simulation [always if(doExtras)] :

ra <- array(dim = c(length(q.set2), length(qb.names), length(p.set2)),
	    dimnames = list(qn2, qb.names, pn2))
ind_ct <- names(system.time(1))[ c(1,3) ] # "user" and "elapsed"
cta <- array(dim = c(length(ind_ct), dim(ra)),
             dimnames = c(list(ind_ct), dimnames(ra)))
for(qq in q.set2) {
  cat("\nq=",qq,"\n======\n")
  for (p in p.set2) {
    cat(sprintf("p=%-7g : >> ",p))
    for (qn in qb.names) {
      cat(" ",qn,": ", sep = '')
      st <- system.time(
          r <- qf[[qn]](.05, p, qq))[ind_ct]
      cta[, qn2[qq == q.set2], qn, pn2[p == p.set2]] <- st
      ra [  qn2[qq == q.set2], qn, pn2[p == p.set2]] <- r
    };  cat("\n")
  }
  showProc.time() ## -- similar for all q
}

    print(summary(warnings())) # many warnings from qbeta() inaccuracies

    save2RDS(list_(ra, cta), file = sfil1)

}## {run simulations} -------------------------------------------------------

## Timings :
1e4*apply(cta, 2:3,    mean)
1e4*apply(cta, c(2,4), mean)

noquote(format01prec(aperm(ra)))

qbeta.R(.05 ,  10/2, 1/2, verbose=2)
qbeta.R(.05 , 100/2, 1/2, verbose=2)
qbeta.R(.05 ,1000/2, 1/2, verbose=2)

showProc.time()#-----------------------------------------------------------------------------

DEBUG <- TRUE ## FIXME

## This is a "version" of qnormUappr() --
## but that has "regular"
##
##	lower.tail=TRUE, log.p=FALSE
##                                        --- MM(2019) to MM(2008): I think this is a bad idea
## ---> ../TODO:  Use *argument*  lower.tail=FALSE for all qnormU*() :  "U" means "Upper tail" !!!!
qnormU <- function(u, lu) {
    if(missing(u))
        qnorm(lu, lower.tail=FALSE,log.p=TRUE)
    else
        qnorm(u,  lower.tail=FALSE)
}

nln <- "\n----------------\n"
for(a in (10^(2:8))/2)
    cat("p=", a, qbeta.R(.05, a, 1/2, low.bnd = 1e-9, qnormU = qnormU),
        nln)
showProc.time()


for(a in (10^(2:16))/2) cat("p=",a, qb <- qbeta.R(.05, a, 1/2, low.bnd = 0), 1-qb,nln)
##            ----  -- SENSATIONAL !                           --------
a <- 5e13; cat("p=",a, qb <- qbeta.R(.05, a, 1/2, low.bnd = 0), 1-qb,nln) #-- problem: 0 outer it.
a <- 5e14; cat("p=",a, qb <- qbeta.R(.05, a, 1/2, low.bnd = 0), 1-qb,nln) #-- problem: 19 outer it
a <- 5e15; cat("p=",a, qb <- qbeta.R(.05, a, 1/2, low.bnd = 0), 1-qb,nln) #-- the last one "working"
a <- 5e16; cat("p=",a, qb <- qbeta.R(.05, a, 1/2, low.bnd = 0), 1-qb,nln) # still fails: init xinbta =1

for(a in (10^(2:16))/2) cat("p=",a, qb <- qbeta.R(.001, a, 1/2, low.bnd = 0), 1-qb,nln)
##                                                ----

## BOTH p & q are BIG: --- sometimes  LONG loop !
for(a in (10^(2:15))/2) cat("p=",a, qb <- qbeta.R(.05, a, a/2, low.bnd = 0), 1-qb,nln)
##            ---- !WOW!
for(a in (10^(2:15))/2) cat("p=",a, qb <- qbeta.R(.05, a, a/2, low.bnd = 0, qnormU = qnormU,
                                                  f.acu = function(...)1e-26),1-qb,nln)

options(op)# revert to usual
showProc.time()


## When using a new  qbeta(.) C-function  with  low.bnd=0.0, up.bnd = 1.0
##--- STILL FAILS at 10^10 ... why ? ...
op <- options(warn=1) # immediate {there are  only 2, now}
for(a in (10^(2:16))/2)
    cat("p=",a, qb <- qbeta(.05, a, 1/2), 1-qb, " relE{pb, .05): ", formatC(pbeta(qb, a, 1/2)/.05 - 1),
        nln)
options(op)

showProc.time()


## Had tests for qbeta() built on AS 109 --- but now R *does* use  AS 109 !
## This will probably all be irrelevant now, see also
## ==> ~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/qbeta109.R



pr.val <- c(.2,.4, 10^(-1:-5)); pr.val <- sort(c(pr.val, 1/2, 1-pr.val))

p.seq <- sort(outer(c(1,3),10^(3:16)))
p.seq  <- c(2^(-15:6),10^(2:18))
p.seq <- c(1e5^c(-4:-1), 10^c(-2:2,10:18))

## q = 1  ---------- know exact formula ...
for(x in pr.val) {
  cat("-----------\nx=", format01prec(x, digits = 6),":\n~~~~~~~~~~\n")
  for(p in p.seq) {
    t0 <- proc.time()[1]
    qb <- qbeta(x, p, 1)  # takes forever
    C1 <- (t1 <- proc.time()[1]) - t0
    pb <- pbeta(x^(1/p), p, 1)
    C2 <- (t2 <- proc.time()[1]) - t1
    fcat(
        ## NO Debug:
        paste0(formatC(p,w = 9), "|",formatC(C1, digits = 2, wid = 5),"|"),
        ## In any case:
        format01prec(qb, digits = 10),
        qb - x^(1/p), qb /( x^(1/p))-1, pb - x
	 ##, f.dig=3, f.wid=8 # << Debug
    )
  }
}

summary(warnings())
showProc.time()

## q = 1/2  (T-distrib.)
## qt is defined as
## qt(x,n) := - sqrt((1  / qbeta( 2*x    , n/2, 1/2) - 1) * n)  ## for x <= 1/2
## qt(x,n) :=   sqrt((1  / qbeta( 2*(1-x), n/2, 1/2) - 1) * n)  ## for x >  1/2
for(x in pr.val) {
  cat("-----------\nx=", format01prec(x, digits = 6),":\n~~~~~~~~~~\n")
  for(p in p.seq[p.seq >= 1/2]) {
    qb0 <- 1/(1+ qt(x/2, df = 2*p)^2/(2*p))
    t0 <- proc.time()[1]
    qb <- qbeta(x, p, 1/2)  # takes forever
    C1 <- (t1 <- proc.time()[1]) - t0
    pb <- pbeta(qb, p, 1/2)
    C2 <- (t2 <- proc.time()[1]) - t1
    ## NO Debug:
    fcat(p,paste("|",t2-t0,"|",sep = ''),format01prec(c(qb,qb0),digits = 10),
    ## Debug: fcat(format01prec(qb,digits=10),
	 paste("|",formatC(qb0-qb,dig = 8,w = 12)),
	 qb/qb0 - 1,  pb/x - 1, f.wid = 12)
  }
}

summary(warnings())
showProc.time()

###  q = p  !!
for(x in pr.val) {
  cat("-----------\nx=", format01prec(x, digits = 6),":\n~~~~~~~~~~\n")
  for(p in p.seq[p.seq >= 1/2]) { #-- otherwise, qt(.) is not defined
    y <- 1/(1+ qt(1-x, 2*p)^2/(2*p))
    qb0 <- 1/2 + ifelse(x < 1/2, -1, 1)* sqrt(1-y)/2  #-- should be = qbeta(..)
    t0 <- proc.time()[1]

    qb <- qbeta(x, p,p)
    C1 <- (t1 <- proc.time()[1]) - t0
    pb <- pbeta(qb, p,p)
    C2 <- (t2 <- proc.time()[1]) - t1
    ## NO Debug:
    fcat(p,paste("|",t2-t0,"|",sep = ''),format01prec(qb,digits = 10),
	 paste("|",formatC(qb0-qb,dig = 8,w = 12)), qb/qb0 - 1,  pb/x - 1, f.wid = 12)
    ## Debug:
    ##fcat(format01prec(qb,digits=10,flag='-'),paste("|",formatC(qb0-qb,dig=8,w=12)),
    ##     qb/qb0 - 1,  pb/x - 1, f.wid=12)
  }
}

summary(warnings())
showProc.time()


###---- 2008-11-12 : another example, not sure if covered above:
###  Letting  a --> 0   --- Distribution only defined for  a > 0 ---
##
curve(qbeta(.95, x, 20), 1e-300, 0.1,           n=1001)# a -> 0 : qbeta() -> 0

## now, "zoom in logarithmically":
curve(qbeta(.95, x, 20), 1e-7, 0.1, log="xy", n=1001, yaxt="n")
if(require("sfsmisc")) { eaxis(2, at= c(10^-c(314,309), axTicks(2, nintLog=15)))
} else { axis(2)
    warning("could not use the cool eaxis() function from CRAN package 'sfsmisc' as that is not yet installed.")
}

## -> clearly should go to 0;
##  now does ... not continually, as with pbeta() warning about bgrat() underflow
##  now without any warnings and continually: not 0 but staying at
unique(qbeta(.95, 10^-(100:5), 20))
##  5.5626846462680034577e-309

## limit :
qbeta(.95, 0, 20)# gave '1' instead of 0 --> now 0  {as for dbeta(*, a=0, *)

## Variation on the above  (a --> 0) :
##                         ---------
## this seems ok; is "the" behavior at large"
x. <- seq(0,1,len=1001); plot(x., pbeta(x., 1e-20,20), type="l")
## but if you log-zoom in ,
## the precision problem is already in pbeta():  ----- fixed on 2009-10-16  (in Skafidia, Greece)
##-- TOMS 708: 'a0 small .. -> apser()
curve(pbeta(exp(x), 1e-16, 10, log.p=TRUE), -10, 0, n=1000)
curve(pbeta(exp(x), 1e-16, 10, log.p=TRUE), -900, 0, n=1000)
                             ## _jumps_ to 0 - because exp() *must* [no bug in pbeta !]
## but this is somewhat similar, *not* using  apser(), but L20, then
## a) (x < 0.29): bup() + L131 bgrat() & w = 1-w1 ~= 1 -- now fixed
## b)  x >=0.29:  bpser()
xpb <- curve(pbeta(exp(x), 2e-15, 10, log.p=TRUE), -21, 0, n=1000)
xpb <- within(as.data.frame(xpb), {exp.x <- exp(x); bp <- exp(x) >= 0.29})
## unfinished
bps <- bpser(2e-15, 10, x = with(xpb, exp.x[bp]), log.p=TRUE, verbose=TRUE)


u <- seq(-16, 1.25, length=501)
op <- mult.fig(mfrow=c(3,5), marP = c(0,-1.5,-1,-.8),
               main = expression(pbeta(e^u, alpha, beta, ~~ log=='TRUE')))$old.par
for(b in c(.5, 2, 5))
    for(a in c(.1, .2, .5, 1, 2)/100) {
        plot(u, (qp <- pbeta(exp(u), a,b, log.p=TRUE)), ylab = NA,
             type="l", main = substitute(list(alpha == a, beta == b), list(a=a,b=b)))
    }
par(op)


summary(warnings())
showProc.time()

###---- Another (or the same?) set of problems:
## From: Martin Maechler <maechler@stat.math.ethz.ch>
## Sender: r-devel-bounces@r-project.org
## To: leydold@statistik.wu-wien.ac.at
## Cc: r-devel@r-project.org
## Subject: Re: [Rd] Buglet in qbeta?
## Date: Wed, 7 Oct 2009 16:04:47 +0200

## >>>>> "JL" == Josef Leydold <leydold@statmath.wu.ac.at>
## >>>>>     on Wed, 7 Oct 2009 14:48:26 +0200 writes:

##     JL> Hi,
##     JL> I sometimes play around with extreme parameters for distributions and
##     JL> found that qbeta is not always monotone as the following example shows.
##     JL> I don't know whether this is serious enough to submit a bug report (as
##     JL> this example is near to the limitations of floating point arithmetic).

## well, it's not near enough the limits to *not* warrant a bug
## report!

## I "know" (I have saved a collection of ..) about several
## accuracy problems of pbeta() / qbeta(), most cases indeed
## "borderline" (at the limits of FP arithmetic) but this one may
## well be a new one.

## A bit more succint than the dump of your numbers is a simple
## graphic

p <- (1:100)/100; qp <- qbeta(p, 0.01, 5)
plot(p,qp, type="o", log = "y") # now fine

## or, already suggesting a fix to the bug:

plot(p,qp, type="o", log = "xy")## now fine
## which is definitely a bug...  though not a terrible one.
## ...

## more experiments -- adapt to log-scale:
if(!exists("lseq", mode="function"))# package "sfsmisc"
   lseq <- function (from, to, length) exp(seq(log(from), log(to), length.out = length))
p <- lseq(.005,1, len = 500) ; qp <- qbeta(p, 0.01,5)
## (no warnings)
plot(p,qp, ylab = "qbeta(p, .01, 5)", type="l", log = "xy")

## p --> 0 even closer -- next problem {but that's easier forgivable ..}
p <- lseq(.0005,1, len = 500)
a <- .01 ; b <- 5
## nomore (gave 36 warnings "qbeta: ..precision not reached ..")
plot(p, (qp <- qbeta(p, a,b)), ylab = expression(qbeta(p, alpha, beta)),
     type="l", log = "xy", main = substitute(list(alpha == a, beta == b), list(a=a,b=b)))

## Seems almost independent of  beta =: b
op <- mult.fig(mfrow = c(3,5), marP = c(0,-.8,-1,-.8),
               main = quote('qbeta(p, <small>, .) vs p  for  '~ {p %->% 0}))$old.par
for(b in c(.5, 2, 5))
    for(a in c(.1, .2, .5, 1, 2)/100) { # last one, a = 0.02 does not underflow
        plot(p, (qp <- qbeta(p, a,b)), ylab = expression(qbeta(p, alpha, beta)), col=2,
             type="l", log = "xy", main = substitute(list(alpha == a, beta == b), list(a=a,b=b)))
    }
par(op)

summary(warnings())
showProc.time()

## x <- qbeta((0:100)/100,0.01,5)
## x
## ....
## order(x)
## 1 2 3 .... 50 51 55 68 72 56 59 ...
## pbeta(x,0.01,5)
## ...
## >> version
## JL> _
## JL> platform       x86_64-unknown-linux-gnu
## JL> arch           x86_64
## JL> os             linux-gnu
## JL> system         x86_64, linux-gnu
## JL> status         Under development (unstable)
## JL> major          2
## JL> minor          11.0
## JL> year           2009
## JL> month          10
## JL> day            07
## JL> svn rev        49963
## JL> language       R
## JL> version.string R version 2.11.0 Under development (unstable) (2009-10-07
## JL> r49963)

## JL> p.s. there are similar results for R-2.9.2 in Windows (with
## JL> different round-off errors).


## --> only vary a == alpha == shape1:

p.qbeta.p0 <- function(p1 = .005, p2 = 1, beta = 1, mark.p = 0.5,
                       alphas = c(.1, .15, .2, .3, .4, .5, .8, 1, 2)/100,
                       n = 500, verbose = getOption("verbose"))
{
    ## Purpose: Plot qbeta() investigations for  small alpha & p --> 0
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  7 Oct 2009, 22:46
    stopifnot(exprs = {
        require("sfsmisc") # lseq(), mult.fig()
        require("DPQ")     # qbetaAppr() ..
        p1 < p2
        beta > 0
        alphas > 0 ; diff(alphas) > 0
        n == round(n) ; n > 2
    })
    p <- lseq(p1, p2, len = n)
    tit <- expression(qbeta(p, alpha, beta) * "  for   small"~alpha~
                          " and " * {p %down% 0})
    cols <- adjustcolor(1:3, 1/2); lwds <- c(2,4,2); ltys <- c(2,1,4)
    op <- mult.fig(length(alphas), main = tit, marP = c(-1,-1,-1,-.6))$old.par
    on.exit(par(op))
    for(a in alphas) {
        matplot(p, cbind(qbetaAppr  (p, a,beta, verbose=verbose),
                         qbeta      (p, a,beta),
                         qbetaAppr.3(p, a,beta)),
                col = cols, lty = ltys, lwd = lwds,
                ylab = "qbeta", type= "l", log = "xy",
                main = substitute(list(alpha == a, beta == b),
                list(a=a, b=beta)))
        if(match(a, alphas) == length(alphas) %/% 2)
           legend("left", c("qbetaAppr", "qbeta", "qbeta.a..3"),
                  col = cols, lty = ltys, lwd = lwds, inset = .02, bty="n")
        if(is.numeric(mark.p) && length(mark.p)) {
            abline(  v = mark.p, lty = 3, lwd = 2, col = "light blue")
            axis(1, at = mark.p, line=-1, lwd = 2, col = "light blue")
        }
    }
    invisible(p)
}

p.qbeta.p0(.05, mark=NULL, verbose=2)
p.qbeta.p0(.0005)
p.qbeta.p0(.00002)
## but for beta <~ 0.5, the approximation seems problematic
p.qbeta.p0(beta = 0.8)
p.qbeta.p0(beta = 0.4)# approx. does not go up to 1
## but then, the "doc" says to use it  for x < 1/2
p.qbeta.p0(beta = 0.2)
p.qbeta.p0(beta = 0.1)## -> first plot, alpha = .001 is "catastrophe" << still not ok
p.qbeta.p0(.00001, beta = .0001)# here, the approx. is completely hopeless
## larger alphas : another problem: the *approximation* is bad!
p.qbeta.p0(.001, 0.5, n=2000, alpha= c(.2, .5, 1, 2))

showProc.time()

p <- lseq(.0001, 1, len = 2000)
b <- 1
tit <- expression(qb == qbeta(p, alpha, beta) * "  for   small"~alpha~" and " * {p %down% 0})
aset <- c(.1, .15, .2, .3, .4, .5, .8, 1, 2)/100
op <- mult.fig(length(aset), main = tit, marP = c(-1,-1,-1,-.6))$old.par
c2 <- adjustcolor(2, 1/2)
for(a in aset) {
    qpA <- qbetaAppr(p, a,b)
     qp <- qbeta    (p, a,b)
    plot (p, qpA, ylab = "qb", type="l", log = "xy", ylim=pmax(1e-20, range(qp, qpA)),
          main = substitute(list(alpha == a, beta == b), list(a=a,b=b)))
    lines(p, qp, col = c2, lwd=3)
    legend("bottom", c("qbetaAppr()", "qbeta()"),
           lty=1, lwd=c(1,3), col=c(palette()[1], c2), bty="n")
}
par(op)


## Ok, let's look at pbeta() "here"
curve(pbeta(x, 0.01, 5, log.p=TRUE), 0, 1,  n=1000) # fine (but not close enough to x = 0)
## log-zoom:
curve(pbeta(x, 0.01, 5, log.p=TRUE), 1e-100, 1, n=10000, log="x") # fine
curve(pbeta(x, 0.01, 5, log.p=TRUE), 1e-250, 1, n=10000, log="x") # fine
## it does *not* seem to be the problem of pbeta
## Ok, the "best" qbeta() is the one matching  pbeta:
p. <- lseq(1e-3, 1, length=200)
qb  <- qbeta       (p., 0.02, 0.01)
qb3 <- qbetaAppr.3(p., 0.02, 0.01)
## qbetaAppr.3() is good for *small* p., but otherwise breaks down *before*
## ------------   *qbeta()*
cbind(p., pbeta(qb,0.02, 0.01), pbeta(qb3,0.02, 0.01))
showProc.time()


## Excursion:  In order to check if the extreme-left tail approximation is ok,
##             we need to compute  log(p*beta(p,q)) accurately for small p
## This suffers from "cancellation" and is for *EXPERIMENTATION*
## --> c12pBeta() further below is for "production"
c12.pBeta <- function(q, c2.only=TRUE) {
    ## Purpose: Compute 1st and 2nd coefficient of p*Beta(p,q)= 1 + c_1*p + c_2*p^2 + ..
    ## ----------------------------------------------------------------------
    ## Arguments: q -- can be a vector
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Oct 2009
    stopifnot(q > 0)
    psi1 <- digamma(1)# == - (Euler's)gamma = -0.5772157
    psiD1 <- trigamma(1)

    c1 <- psi1 - digamma(q)
    ## Note:  psi(1)^2 - psi(q)^2 = (psi(1) - psi(q))*(psi(1)+psi(q)) =
    ## and (trigamma(1) - trigamma(q) + c1*(psi1+psi(q)))/2 - psi(q) * c1
    ##....
    c2 <- (psiD1 - trigamma(q) + c1^2)/2
    ## Note: cancellation for  q -> 0 !! ---> c12pBeta()  below
    ##  digamma(q) = psi (q) ~= -1/q + psi(1) + q*psi'(1) + q^2/2*psi''(1) +...
    ## trigamma(q) = psi'(q) ~= 1/q^2 + psi'(1) + q*psi''(1) ...
    ## ==> 2*c_2 = psi'(1) - psi'(q) + (psi(q) - psi(1))^2  =
    ##  = - 1/q^2 -q*psi''(1) -.. + (-1/q + q*psi'(1) + q^2/2*psi''(1)+..)^2 =
    ##  = -q*psi''(1) -2*psi'(1) - q*psi''(1) = 2*(-psi'(1)-q*psi''(1))+ O(q^2)
    ## ==>  c_2(q) = -(psi'(1) + q*psi''(1)) + O(q^2)
    ##      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(c2.only) c2 else cbind(c1 = c1, c2 = c2)
}
c12.pBeta(1e-5,FALSE)
summary(warnings())

showProc.time()

##         c1        c2
## [1,] 1e+05 -1.644912
trigamma(1)
## [1] 1.644934
## which confirms the expansion of c_2(q) above

curve(c12.pBeta, .001, 100) ##->
curve(c12.pBeta, 1, 1000, log="xy") ## dominated by  -digamma(q)^2 / 2
## Cancellation problem when q -> 0 as well ... probably not really relevant
## *but* fixed in c12pBeta() below:
curve(c12.pBeta, 5e-9, 100, log="x")

tit.pBeta <- function() {
    title("c12.pBeta(x) and its (log-)linear approximation")
    mtext(substitute(c[1] + c[2]*x == C1 + C2*x,
                     setNames(as.list(cq12), c("C1", "C2"))), col=2)
}
cq12 <- -psigamma(1, 1:2) ## == c( -digamma(1), -trigamma(1) )
curve(cq12[1] + cq12[2]*x, add=TRUE, col=2, n=1000)
tit.pBeta()

## a bit zoomed:
curve(c12.pBeta, 1e-8, .2, log="x")
curve(cq12[1] + cq12[2]*x, add=TRUE, col="red", n=1000); tit.pBeta()
## difference only:
curve(c12.pBeta(x) - (cq12[1] + cq12[2]*x),
      1e-8, .2, log="x", ylim = .01*c(-1,1)); abline(h=0, lty=3, col="gray60")
## | . - . | log-scaled:
curve(abs(c12.pBeta(x) - (cq12[1] + cq12[2]*x)),
      1e-8, .2, log="xy")
abline(v= 4e-4, col="lightblue")
##==> use  q = 4e-4 as cutoff

qq <- lseq(1e-6, 1e-2, length= 64)
1 - (cq12[1] + cq12[2]*qq) / c12.pBeta(qq)
plot(qq, 1 - c12.pBeta(qq)/(cq12[1] + cq12[2]*qq), type="o", log="x")
plot(qq, abs(1 - c12.pBeta(qq)/(cq12[1] + cq12[2]*qq)), type="o", log="xy")

showProc.time()

##  Compute  log(p*beta(p,q)) accurately for small p
c12pBeta <- function(q, c2.only=FALSE, cutOff.q = 4e-4) {
    ## Purpose: Compute 1st and 2nd coefficient  of   p*Beta(p,q) = 1 + c_1*p + c_2*p^2 + ..
    ##  NB: The Taylor expansion of   log(p*Beta(p,q))    =     c_1*p + d_2*p^2 + ..
    ##      is simpler: c_1=d_1, and d_j = \psi^(j)(1) - \psi^(j)(q) = psigamma(1,j) - psigamma(q,j)
    ## ----------------------------------------------------------------------
    ## Arguments: q -- can be a vector
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Oct 2009
    stopifnot(q > 0, cutOff.q >= 0)
    psi1 <- digamma(1)  # == - (Euler's)gamma = -0.5772157
    psiD1 <- trigamma(1)# ==  pi^2 / 6

    ## vectorize in q:
    r <- q
    if(!c2.only)
        c1 <- psi1 - digamma(q)
    if(any(sml <- q < cutOff.q)) {
        ## To avoid cancellation, use Taylor expansion for small q:
        ##  digamma(q)= psi (q) ~= -1/q + psi(1) + q*psi'(1)+q^2/2*psi''(1)+..
        ## trigamma(q)= psi'(q) ~= 1/q^2+ psi'(1)+ q*psi''(1) ...
        ## ==> 2*c_2 = psi'(1) - psi'(q) + (psi(q) - psi(1))^2  =
        ##  = - 1/q^2 -q*psi''(1) -.. + (-1/q + q*psi'(1) + q^2/2*psi''(1)+..)^2=
        ##  = -q*psi''(1) -2*psi'(1) - q*psi''(1) = 2*(-psi'(1)-q*psi''(1))+ ..
        ## ==>  c_2(q) = -(psi'(1) + q*psi''(1)) + O(q^2)
        ##      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        r[sml] <- -(psiD1 + psigamma(1,2) * q[sml])
    }
    if(!all(sml)) { ## regular case:
        i.ok <- !sml
        q <- q[i.ok]
        c1. <- if(c2.only) psi1 - digamma(q) else c1[i.ok]

        ## Note:  psi(1)^2 - psi(q)^2 = (psi(1) - psi(q))*(psi(1)+psi(q)) =
        ## and (trigamma(1) - trigamma(q) + c1*(psi1+psi(q)))/2 - psi(q) * c1
        ##....
        r[i.ok] <- (psiD1 - trigamma(q) + c1.^2)/2
    }
    if(c2.only) r else cbind(c1 = c1, c2 = r)
}

tit <- quote(c[2](q) ~~~ "from Beta expansion" ~~ p*Beta(p,q) %~~% 1 + c[1]*p + c[2]*p^2)
doX <- function(cutOff.q = eval(formals(c12pBeta)$cutOff.q)) {
    eaxis(1); py <- par("usr")[3:4]
    abline(h=0, v=1, lty=3)
    text(1, py[1]+0.95*diff(py), quote(q==1), adj=-.005)
    abline(v = cutOff.q, col=2, lty=2, lwd=2)
    text(cutOff.q, py[1]+0.9*diff(py),
         paste0("cutOff.q=",formatC(cutOff.q)), col=2, adj=-0.05)
}
curve(c12pBeta(x,TRUE), 5e-9, 100, log="x", n=1025, xaxt='n',
      xlab=quote(q), main = tit); doX()#--> no cancellation problems
showProc.time()

## zoom in
curve(c12pBeta(x,TRUE), 3e-4, 5e-4, log="x", n=1025, xaxt='n',
      xlab=quote(q), main = tit); doX()#--> no break..
showProc.time()

## zoom out
curve(c12pBeta(x,TRUE), 1e-100, 1e110, log="x", n=1025, xaxt='n',
      xlab=quote(q), main = tit); doX()# .. neither problems ..

curve(c12pBeta(x,TRUE), 1e-100, 1e110, log="x", n=1025, xaxt='n', ylim = c(-2,1),
      xlab=quote(q), main = tit); doX()# .. neither problems ..

summary(warnings())
showProc.time()

## In order for the first oder approximation  p*beta(p,q) = 1 + c_1*p to be ok,
## that  c_2*p^2 << 1, i.e.  p^2 < eps.c / c_2  <=>  p < sqrt(eps.c / c_2)
##                                                   =====================
## --> draw the cutoff line below ( --> 'p.cut')

p.pBeta <- function(q, p.min=1e-20, p.max=1e-6, n=1000, do.ord2=TRUE, log = "xy", ylim=NULL)
{
    ## Purpose: Plot   log(p * beta(p,q))   for small p  -- to visualize cancellation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Oct 2009, 09:16
    stopifnot(0 < p.min, p.min < p.max, n > 100, n == round(n),
              is.character(log), length(log) == 1,
              length(q) == 1, 0 < q)
    sgn <- if(q > 1) -1 else 1
    curve(sgn*log(x*beta(x,q)), p.min, p.max, n=n, log=log, ylim=ylim,
          col="gray", lwd=3, ylab="", xlab = expression(p),
          axes = FALSE, frame.plot = TRUE, # rather use eaxis() below
          main = substitute("Accurately computing "~ log(p %.% B(p,q)) ~
          "   for  "~ {p %->% 0} *",     "~ (q == Q), list(Q=format(q,digits=15))))
    eaxis(1); eaxis(2)
    curve(sgn*(log(x)+lbeta(x,q)), add = TRUE, col="red3", n=n)
    ## Here comes the 1st order cancellation solution
    c1 <- sgn*( digamma(1) -  digamma(q))
    c2 <- sgn*(trigamma(1) - trigamma(q))/2
    colc1 <- adjustcolor("skyblue", 1/2)
    c1y <- curve(c1 * x, add = TRUE, col=colc1, lwd=2, n=n)$y
    ## the 2nd order approx.
    if(do.ord2) {
        colc2 <- "#60C07064" # opaque (!)
        c12y <- curve(x*(c1 + c2* x), add = TRUE, col=colc2, lwd=3, n=n)$y
    }
    ## rather:
    collog1 <- adjustcolor("forestgreen", 1/2)
    lc1y <- curve(sgn*log1p(sgn*c1*x), add = TRUE, col=collog1, lty=5, lwd=3, n=n)$y
    ## the cutoff, from where 1st order should be perfect:
    p.cut <- sqrt(.Machine$double.eps / abs(c2))
    txt.cut <- "cutoff(1st order approx.)"
    cat(txt.cut,": ", formatC(p.cut),"\n", sep="")
    ## FIXME: These *badly* fail, if we do NOT have  log = "xy" :
    ux <- par("xaxp")[1:2] # e.g. [1e-20, 1e-6]
    uy <- par("yaxp")[1:2]
    in.x <- function(f) { stopifnot(0 <= f, f <= 1); 10^c(c(1-f, f) %*% log10(ux)) }
    in.y <- function(f) { stopifnot(0 <= f, f <= 1); 10^c(c(1-f, f) %*% log10(uy)) }
    if(ux[1] <= p.cut & p.cut <= ux[2]) {
        abline(v = p.cut, col = "blue", lwd=2, lty=2)
        text(p.cut, in.y(.96), txt.cut, col="blue")
    }
    else if(p.cut > ux[2]) text(in.x(.98), in.y(0.5), paste("-->", txt.cut), col = "blue", lwd=2, lty=2)
    else if(p.cut < ux[1]) text(in.x(.02), in.y(0.5), paste(txt.cut, "<--"), col = "blue", lwd=2, lty=2)
    leg <-
        if(sgn == -1) expression(-log(p*B(p,q)), -(log(p) + log(B(p,q))), -(psi(1)-psi(q))*p, -log1p(c[1]*p))
        else          expression(log(p*B(p,q)), log(p) + log(B(p,q)), (psi(1)-psi(q))*p, log1p(c[1](q)*p))
    if(do.ord2) { leg <- leg[c(1:4,4)]; leg[[4]] <- quote(p*(c[1](q) + c[2](q)* p)) }
    legend("topleft", legend= leg, inset=.01,
           lwd= c(3,1,2,if(do.ord2)3, 3),
           col= c("gray","red3", colc1, if(do.ord2)colc2, collog1))
    mtext(R.version.string, side=4, cex = 3/5, adj=0)
}

showProc.time()

p.pBeta(q = 1.1) ##q should not matter much; q very close to 1 will be another challenge
## zoom in:
p.pBeta(1.1, p.max = 1e-12)
p.pBeta( 11, p.max = 1e-11)## log(..) is not much worse here
p.pBeta(101, p.max = 1e-11)## !! interestingly, here log(..) is not worse
p.pBeta(1001, p.max = 1e-11)## both versions are identical

showProc.time()
## q -> 1
p.pBeta(1.001, p.max = 1e-10)
p.pBeta(1+1e-5, 1e-15,  1e-8)# problems start earlier!
p.pBeta(1+1e-7, 1e-11,  1e-6)# problems start even earlier!
p.pBeta(1+1e-9, 1e-9,  1e-4)# problems start even earlier!
##
showProc.time()
## q < 1
p.pBeta(0.9)
p.pBeta(0.1, p.max = 1e-8)
p.pBeta(0.01)
p.pBeta(0.0001)
p.pBeta(0.0001,1e-23, 1e-4)
p.pBeta(0.0001,1e-12, 1e-2)
p.pBeta(0.0001,1e-6, 1e-3)
p.pBeta(0.0001,1e-6, 1e-1)

summary(warnings())
showProc.time()

if(!require("Rmpfr")) quit("no")
##--===============---====------------ Rmpfr needed from here -----------------------------


p.err.pBeta <- function(q, p.min=1e-12, p.max=1e-6, n=1000,
                        kind = c("relErr", "absErr", "error"))
{
    ## Purpose: Plot   log(p * beta(p,q))   for small p  -- to visualize cancellation
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 30 Oct 2009, 09:16
    stopifnot(0 < p.min, p.min < p.max, n > 100, n == round(n),
              length(q) == 1, 0 < q)

    c12 <- c12pBeta(q)
    d2 <- (trigamma(1) - trigamma(q))/2
    c1 <- c12[1] ## == digamma(1) - digamma(q)
    c2 <- c12[2]

    p <- lseq(p.min, p.max, length=n)

    T0 <- log(p*beta(p,q))# which we know has its problems, too
    T1 <- c1*p
    T2 <- p*(c1 + d2*p)
    T3 <- log1p(T1)
    T4 <- log1p(p*(c1 + c2*p))

    ## Better: compute "true value" very accurately:
    stopifnot(require("Rmpfr"))
    p. <- mpfr(p, prec=200)
    tt <- log(p.*beta(p.,q))

    p.cut <- sqrt(.Machine$double.eps / abs(c2))
    cat("cutoff :", formatC(p.cut),"\n")

    err.mat <- as(tt - cbind(T0, T1, T3, T2, T4),
                  "matrix")## classic numeric matrix
    leg <- expression(T - log(p*beta(p,q)) ~ " (double)",
                      T - c[1]*p,		#  T1
                      T - log1p(c[1]*p),	#  T3
                      T - (c[1]*p+d[2]*p^2),	#  T2
                      T - log1p(c[1]*p+c[2]*p^2))# T4
    stopifnot(length(leg) == ncol(err.mat))
    main <- paste("log(p*beta(p,q))  approx. errors,   q = ", formatC(q))
    ## `col` below is from
    ##     dput(RColorBrewer::brewer.pal(length(leg), "Dark2"))
    col <- c("#1B9E77", "#D95F02", "#E7298A", "#7570B3", "#66A61E")
    ## put the color which is "most gray" (i.e equal R,G,B values) first:
    nc <- length(col)
    cm <- col2rgb(col)
    i1 <- which.min(colSums(sweep(cm, 2, colMeans(cm))^2))
    if(i1 > 1) col <- col[c(i1:nc, 1:(i1-1))]
    m.matplot <- function(x,y, log, main, ..., LEG, xy.leg) {
        log.. <- strsplit(log, "")[[1]]
        log.x <- any("x" == log..)
        log.y <- any("y" == log..)

        matplot(x,y, log=log, ..., xlab = quote(p), type = "l", col=col, main=main,
                xaxt = if(log.x) "n" else "s",
                yaxt = if(log.y) "n" else "s")
        if(log.x) eaxis(1, nintLog = 15)
        if(log.y) eaxis(2, nintLog = 20) # "FIXME" eaxis(): default nintLog=10 does not cover well
        legend(xy.leg, LEG, lty=1:5, col=col, inset=.01)
    }
    kind <- match.arg(kind)
    switch(kind,
           "error" =
               m.matplot(p, err.mat, log="x", main=main,
                         ylab = "log(p*beta(p,q)) - 'approx.'",
                         LEG = leg, xy.leg = "left"),
           "absErr" = {

               m.matplot(p, abs(err.mat), ## <- absolute -- log scale
                         log="xy", main=main,
                         ylab = "|log(p*beta(p,q)) - 'approx.'|",
                         LEG = as.expression(lapply(leg,
                             function(.) substitute(abs(EXPR), list(EXPR= .)))),
                         xy.leg = "topleft")
           },
           "relErr" = {
               m.matplot(p, abs(err.mat/as(tt,"numeric")), ## <- RELATIVE errors
                         log="xy", main = paste(
                                   "log(p*beta(p,q))  RELATIVE approx. errors,   q = ",
                                   formatC(q)),
                         ylab = "|1 - 'approx'/log(p*beta(p,q))|",
                         LEG = as.expression(lapply(leg,
                             function(.) substitute(abs(EXPR), list(EXPR= .)))),
                         xy.leg = "topleft")
           },
           stop("invalid 'kind' : ", kind)
           )# end{switch}

    abline(v = p.cut, col = "blue", lwd=2, lty=2)
    mtext(R.version.string, side=4, cex = 3/5, adj=0)
}

## with q >> 1,  the  log1p() approx. are *worse* (see more below) !
p.err.pBeta(21, 1e-18) # cutoff : 5.527e-09
p.err.pBeta(21, , 1e-3)
summary(warnings())
showProc.time()

if(doExtras) { # each plot is somewhat large & expensive  ..-------------------
## hmm, with q >> 1,  the  log1p() approx. are *worse* ..
p.err.pBeta(q = 2.1, kind="absErr")
p.err.pBeta(q = 2.1, 1e-13, kind="absErr")
p.err.pBeta(1.1, 1e-13, 0.1)##  log1p() still very slightly *worse*
## here, q < 1 and log1p() are already (slightly) better
p.err.pBeta(0.8)
p.err.pBeta(0.8, 1e-10, 1e-2)
p.err.pBeta(0.1)  #  log1p(<2 trms>) is best;  zoom out (to guess "cutoff_2")
p.err.pBeta(0.1, 1e-15, 1e-3)
p.err.pBeta(0.001,, 1e-2)## Here, log1p(<2 terms>) is clearly best
## reminding of the test  qbeta(1e-300, 2^-12, 2^-10):
p.err.pBeta(2^-10, p.max = 2^-10); abline(v=2^-12, lty=3, col="gray20")
## --> ok, the (*, 2^-12, 2^-10)  is *not* yet critical

##--->  q < 1 ==>  the log1p() versions are superior
##--->  q > 1 ==>  the direct  versions are superior
##  in both cases:  The quadratic term is an order of magnitude better
##  the above  p.cut is "order-of-magnitude" correct, but not really numerically ok!

print(summary(warnings()))
showProc.time()
}# only if(doExtras) ----------------------------------------------------------------------

pBeta <- function(p, q)
{
    ## Purpose: Compute   log(p * beta(p,q))  accurately also for small p
    ## ----------------------------------------------------------------------
    ## Arguments: p: (vector)  q: scalar
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 31 Oct 2009, 18:38

    stopifnot(0 < p, 0 < q, length(q) == 1, (n <- length(p)) > 0)

    if(q == 1) return(0*p)

    c12 <- c12pBeta(q)
    d2 <- (trigamma(1) - trigamma(q))/2
    c1 <- c12[1] ## == digamma(1) - digamma(q)
    c2 <- c12[2]

    ## the cutoff  from where the approximation is "perfect"
    p.cut <- sqrt(.Machine$double.eps / abs(d2))
    r <- p
    if(any(sml <- p <= p.cut)) {
        p. <- p[sml]
        r[sml] <- if(q < 1) log1p(p.*(c1 + c2*p.)) else p.*(c1 + d2*p.)
    }
    if(any(ok <- !sml)) {
        p <- p[ok]
        r[ok] <- log(p * beta(p, q))
    }
    r
}

###  More on  log(p * Beta(p,q))  for  p << 1
###           ==================


## 1. Slightly nicer   Gamma(.)  plot than in  example(Gamma) : ----------
op <- options(warn = -1)
x <- sort(c(seq(-3, 4, length.out = 201), outer(0:-3, (-1:1)*1e-6, "+")))
plot(x, gamma(x), ylim = c(-20,20), col = "red", type = "l", lwd = 2,
     main = expression(Gamma(x)))
abline(h = 0, v = -3:0, lty = 3, col = "midnightblue")
## Warning .. in gamma(x) : NaNs produced
axis(4, at=factorial(1:3), las=2, cex.axis=3/4)
abline(h=1, lty=4, col="thistle")
for(t in c("p","h")) points(1:4, gamma(1:4), pch=4, type=t, lty=2)

## Now try, comparing with Rmpfr exact values:
if(!require("Rmpfr")) { message("Cannot attach 'Rmpfr' package"); q("no") }

a. <- 1/256 # is exact!
b. <- 1/2   # "   "
a <- mpfr(a., 512)
b <- mpfr(b., 512)
(aBa1 <- log(a*beta(a,b)))
## 0.0053902550661545715269612410527830537772316302960547912531493485401255362864762019....
(aBa2 <- log(a)+lbeta(a,b))#
## 0.0053902550661545715269612410527830537772316302960547912531493485401255362864762019....
## how close?
asNumeric((aBa1-aBa2)*2/(aBa1+aBa2)) # -3.426748e-152; (NB 512 bits ~= 154.1 digits)

## Double precison accuracy
dput(log(a. * beta(a.,b.)))
## 0.00539025506615481
dput(log(a.)+lbeta(a.,b.))# (to my surprise!) slightly worse
## 0.00539025506615509
asNumeric( log(a. * beta(a.,b.))/(log(a)+lbeta(a,b)) - 1) # 4.504144e-14
asNumeric((log(a.)+lbeta(a.,b.))/(log(a)+lbeta(a,b)) - 1) # 9.653358e-14

psigamma(1) # -0.5772157 = Euler's gamma

## Use my approx:
a.* (psigamma(1) - psigamma(b.)) ## 0.005415212
a.* (psigamma(1) - psigamma(b.) + a./2*(psigamma(1, d=1) - psigamma(b., d=1)))
##  0.005390113 {clearly better, still only ~ 4 correct digits

showProc.time()
