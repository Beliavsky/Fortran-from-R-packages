
options(digits=16)
Aeq <- function(x,y) all.equal(x,y, tol = 1e-14)
stopifnot(exprs = {
    identical(ppois(1.5, 2), ppois(1, 2)) # by internal code definition
    Aeq(ppois(1,2), sum(dpois(0:1, 2)))
})

## From: Matthias Kohl <Matthias.Kohl@uni-bayreuth.de>
## To: r-help <r-help@stat.math.ethz.ch>
## Subject: [R] Exactness of ppois
## Date: Thu, 15 Jan 2004 13:55:22 +0000

##- by checking the precision of a convolution algorithm, we found the
##- following "inexactness":
##- We work with R Version 1.8.1 (2003-11-21) on Windows systems (NT, 2000, XP)

## Try the code:

## Kolmogorov distance between two methods to
## determine P(Poisson(lambda)<=x)
Kolm.dist <- function(lam, eps) {
    x <- 0:qpois(eps, lambda = lam, lower.tail=FALSE)
    max(abs(ppois(x, lambda = lam) - cumsum(dpois(x, lambda = lam))))
}
str(erg <- optimize(Kolm.dist, lower = 800, upper = 2000,
                    maximum = TRUE, eps = 1e-15), digits=12)
## max at  lambda = 1262.6, but value is simply 2.22e-16 ...
## i.e. very small! ==> no problem anymore
Kolm.dist(lam = erg$max, eps = 1e-14)

Kolm1.dist <- function(lam, eps) {
    x <- 0:qpois(eps, lambda = lam, lower.tail=FALSE)
    which.max(abs(ppois(x, lambda = lam) - cumsum(dpois(x, lambda = lam))))
}
(x. <- Kolm1.dist(lam = erg$max, eps = 1e-15))
## 1360,  was  1422 # was '1001' ## or '1301' or .. :  in any case ==  "alphlimit + 1" !!



##- So for lambda=977.8 and x=1001 we get a distance of about 5.2e-06.
##- (This inexactness seems to hold for all lambda values greater than
##- about 900.)

## MM: --- now all (we see here) is fine : Morten Welinder's pgamma() since ~ 2005:
lam <- 977.8
(p1 <- ppois(1001, lam))
(p2 <- sum(dpois(0:1001, lam)))
1 - p1/p2# -2.13e-14 , was -6.7e-6, then 0
stopifnot(abs(1 - p1/p2) < 1e-12)

## new pgamma(alphalimit = 1300):
lam <- 1274.487
(p1 <- ppois(1301, lam))
(p2 <- sum(dpois(0:1301, lam)))
1 - p1/p2# -1.942e014,  was -5.13e-6, then -2.22e-16,
stopifnot(abs(1 - p1/p2) < 1e-12)

## Also:
lam <- 1274.487
x <- 1301
## This is by R's definition of  ppois() :
stopifnot(
    pgamma(lam, x+1, lower.tail=FALSE) == print( ppois(x,lam) )
)


##--- More generally:
library(DPQ)
(doExtras <- DPQ:::doExtras())

## dyn.load("/u/maechler/R/MM/NUMERICS/dpq-functions/ppois-direct.so")

## "TODO": currently quite a few checks already happen in the examples in
## >> ../man/ppoisson.Rd <<
ppD0 <- function(x, lambda, all.from.0 = TRUE)
    cumsum(dpois(if(all.from.0) 0:x else x,
                 lambda=lambda))

(pE  <- ppoisErr(900))
(pEd <- ppoisErr(900, ppFUN = ppD0))
options(digits = 7)# back to normal
all.equal(pE, pEd) # not quite: 'x0' differs slightly

## Large Lam  where  exp(-lambda) underflows to 0 (even in 'long double'):

## "the minimal" large lambda is
okLD <- okLongDouble(999)
notXorLD <- !okLD || !doExtras
Lam  <- 11400 ## but we choose a slightly larger, where one sees more
set.seed(12)
for(Lam in c(11400 + c(0, sort(round(rlnorm(10, 4, 2), 1))))) {
    M    <- ceiling((Lam + 4*sqrt(Lam))/50)*50
    tit <- sprintf("Lam = %g; M = %7.0f", Lam, M)
    cat(tit, ":\n")
    kM <- 0:M
    pp. <- ppoisD(M, lambda=Lam) ## with current DEBUG output vives
    ## ppoisD(*, lambda=11600): expl(-ldlam)=0= 0 ==> llam=9.35876, exp_arg=-11600
    ##  .. i=30, finally new f = expl(exp_arg = -11393.9) = 4.95747e-4949 > 0
    ## ____________or________________
    ## ppoisD(lambda=11444, expl(-ldlam)=0= 0 ==> llam=9.34522, exp_arg=-11444
    ##  .. i=6, finally new f = expl(exp_arg = -11394.5) = 2.69745e-4949 > 0
    pp <-  ppois(kM, lambda=Lam)
    pp0 <- cumsum(dp <- dpois(kM, lambda=Lam))
    if(doExtras)
        print(system.time(
            ppSL <- ppoisD(kM, lambda=Lam, all.from.0=FALSE)
        ))
    ip <- dp > 0
    plot(kM[ip], pp[ip], type="l", main = tit)
    plot(kM[ip], pp[ip], type="l", main = tit, log = "y")# LOG: here *no* warnings
    lines(kM[ip], pp0[ip], col=2)
    lines(kM[ip], pp.[ip], col=adjustcolor("blue", 1/2), lwd=4)
    abline(v = Lam, lty=2, col="gray")
    ##
    plot (kM[ip], pp.[ip] - pp[ip], xlab = "k", type="l", main = tit)
    lines(kM[ip], pp0[ip] - pp[ip], col = adjustcolor("red", 1/2), lwd = 2)
    if(doExtras) lines(kM[ip], ppSL[ip] - pp[ip], col = adjustcolor(3, 1/2), lwd=3)
    abline(v = Lam, lty=2, col="gray")
    stopifnot(exprs =
      if(doExtras) {
        !okLD    || all.equal(pp,  pp. , tol = 1e-12)
        notXorLD || all.equal(pp,  ppSL, tol = 1e-12)
        notXorLD || all.equal(pp., ppSL, tol = 1e-14)# both ppoisD()  "fast" or "slow" should be very close
      } else {
        !okLD    || all.equal(pp,  pp. , tol = 1e-12)
      })
}



## MM(2018-08):  'alphLim' below must have been a way to set 'alphlimit'  in pgamma()'s C code.
## ----------
##
## alphLim <- 1000 # old
alphLim <- 100000  # new --- in R's C code for pgamma() since R 1.8.0

## when using ppoisErr <- ppoisErr1 {much slower ==> save here:
(sdir <- system.file("safe", package="DPQ"))
sfil <- file.path(sdir, "tests_ppoisErr.rda")
if(!doExtras && file.exists(sfil)) {
    cat(sprintf("load(%s) creates ", sfil), load(sfil), "\n")
} else {
    ##             2^20 is much too large and the last few take much time!
    l2ex <- if(doExtras) seq(1, 15, by=1/8) else seq(1, 12, by=1/4)
    errL <- lams <- 2^l2ex
    ## Use for() loop instead of *apply() to see progress and time of each:
    cat(head <- paste(
            "  lambda |     errLam |  user  system elapsed [x 1000, i.e. ms]",
            "---------+------------+---------------------\n", sep="\n"))
    for(i in seq_along(lams)) {
        st <- system.time(errL[i] <- ppoisErr(lams[i]))
        cat(sprintf("%8.2f | %10.4g |", lams[i], errL[i]),
            paste(sprintf("%5.0f", 1000*as.vector(st)[1:3]), collapse="   "),
            "\n")
    } ; cat(head)
    save(lams,errL,alphLim,  file = sfil)
}

plot(lams, abs(errL), log = "xy", xaxt="n", xlab = "lambda",
     main = paste("|rel.Error { ppois(x,lambda) } |    (alphlimit=",
                  formatC(alphLim)," in pgamma)"))
abline(v= alphLim, lty=3)
at.l <- c(outer(c(1,2,5),10^pretty(log10(lams))))
at.l <- at.l[10^par("usr")[1] < at.l & at.l < 10^par("usr")[2]]
## Improve  axis(1, at=at.l, las=2) :
axis(1, at=at.l, labels=NA)
u.y <- par("usr")[3:4]
yA <- function(e) 10^c(c(1+e, -e) %*% u.y)
for(ll in at.l) {
    if(ll < 1e4) { # normal
        yt <- yA(.03); srt <- 0 ;  adj <- 0.5
    } else { # sloping
        yt <- yA(.03); srt <- -30 ; adj <- c(0.25, 0.5)
    }
    text(ll, yt, formatC(ll), srt = srt, adj = adj, xpd = NA)
}

lm1 <- lm(log10(abs(errL)) ~ log10(lams), subset = errL != 0 & lams < alphLim - 150)
##                                       # 855.1, 975.5 already bad  ^^^^^
abline(lm1, col = "blue")
## -----------------------
cat("#{lams > alphLim} : ", sum(lams > alphLim),"\n")
##
if(sum(lams > alphLim) >= 2) { ## always FALSE nowadays
lm2 <- lm(log10(abs(errL)) ~ log10(lams), subset = lams > alphLim)
##                          # alphlimit :                 ^^^^
abline(lm2, col = "purple")
## where do the lines cross:
c1 <- coef(lm1)
c2 <- coef(lm2)
x. <- 10^print(log10x <- (c1[1]-c2[1])/(c2[2]-c1[2]))
x. # 684'300
## draw segment to cross point:
points(x., 10^predict(lm1, new=data.frame(lams=x.)), type = 'h', col=2)
text(x., 1e-16, paste("lambda=",formatC(x.)), srt=90, adj = c(0,0), col = 2)
} # cannot draw these here-----------------------------

