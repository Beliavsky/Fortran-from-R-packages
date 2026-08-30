library("mvtnorm")

chk <- function(...) isTRUE(all.equal(...))

## Showing the TVPACK() gives *NON*-random results:
(cor1 <- toeplitz(c(1, 1/4, -1/8)))
(up1  <- c(1/4, 7/4, 5/8))
d <- length(up1) # = 3
pmvt.. <- function(df, algorithm)
    vapply(df, function(df) pmvt(upper=up1, corr=cor1, df=df, algorithm=algorithm),
           numeric(1))

dfs <- 1:9
pmvt_TV.7 <- replicate(7, pmvt..(dfs, TVPACK()))

stopifnot(identical(unique(c(pmvt_TV.7 - pmvt_TV.7[,1])), 0))
(pmvt.TV. <- pmvt_TV.7[,1])
(pmvt.TV  <- pmvt..(dfs, TVPACK(1e-14)))# has no effect here
chk(max(abs(pmvt.TV - pmvt.TV.)), 0) ## all 0 {unexpectedly ??}


set.seed(47) ## and default algorithm: -> *random* result
pmvt_7 <- replicate(7, vapply(dfs, function(df) pmvt(df=df, upper=up1, corr=cor1), numeric(1)))
## relative errors
relE <- 1 - pmvt_7 / pmvt.TV
rng.rE <- range(abs(relE))
stopifnot(1e-6 < rng.rE[1], rng.rE[2] < 7e-4)
stopifnot(chk(
    colMeans(abs(relE)),
    c(88, 64, 105, 73, 52, 90, 87)*1e-6, tol= 1e-3))


set.seed(29)

########################################################################
## 3 dim example
corr <- cov2cor(crossprod(matrix(runif(9,-1,1),3,3))+diag(3))
df <- rpois(1,3)+1

## central t distribution (-Inf,upper)
ctrl <- GenzBretz(maxpts = 2500000, abseps = 0.000001, releps = 0)
upper <- rexp(3,1)
pmvt(upper=upper, corr=corr, df = df, algorithm = ctrl)
pmvt(upper=upper, corr=corr, df = df, algorithm = TVPACK())

## central t distribution (lower,Inf)
lower <- -rexp(3,1)
pmvt(lower=lower, upper=rep(Inf,3), corr=corr, df = df, algorithm = ctrl)
pmvt(lower=lower, upper=rep(Inf,3), corr=corr, df = df, algorithm = TVPACK())

## non-central t (not possible for TVPACK)
delt <- rexp(3,1/10)
upper <- delt+runif(3)
ctrl <- GenzBretz(maxpts = 2500000, abseps = 0.000001, releps = 0)
pmvt(upper=upper, corr=corr, df = df, algorithm = ctrl, delta = delt)
tools::assertError(pmvt(upper=upper, corr=corr, df = df, algorithm = TVPACK(), delta = delt))

## central mvn (-Inf, upper)
upper <- rexp(3,1)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK())

## central mvn (lower, Inf)
lower <- rexp(3,5)
pmvnorm(lower=lower,upper=rep(Inf, 3), corr=corr, algorithm = ctrl)
pmvnorm(lower=lower,upper=rep(Inf, 3), corr=corr, algorithm = TVPACK())

## non-central mvn
delt <- rexp(3,1/10)
upper <- delt+rexp(3,1)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl,     mean = delt)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK(), mean = delt) # should not error

########################################################################
## 2 dim example
corr <- cov2cor(crossprod(matrix(runif(4,-1,1),2,2))+diag(2))
upper <- rexp(2,1)
df <- rpois(1, runif(1, 0, 20))

## central t (-Inf, upper)
pmvt(upper=upper, corr=corr, df = df, algorithm = ctrl)
pmvt(upper=upper, corr=corr, df = df, algorithm = TVPACK())

## central t (lower, Inf)
pmvt(lower=-upper, upper=rep(Inf, 2), corr=corr, df = df, algorithm = ctrl)
pmvt(lower=-upper, upper=rep(Inf, 2), corr=corr, df = df, algorithm = TVPACK())

## non-central t
delt <- rexp(2,1/5)
upper <- delt+rexp(2,1)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl, mean = delt)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK(), mean = delt)

########################################################################
## comparison with Miwa
## 2d
corr <- cov2cor(crossprod(matrix(runif(4,-1,1),2,2))+diag(2))
upper <- rexp(2, 1)

pmvnorm(upper=upper, corr=corr, algorithm = Miwa(steps=128))
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK())

## 3d
corr <- cov2cor(crossprod(matrix(runif(9,-1,1),3,3))+diag(3))
upper <- rexp(3, 1)

ctrl <- Miwa(steps=128)
pmvnorm(upper=upper, corr=corr, algorithm = ctrl)
pmvnorm(upper=upper, corr=corr, algorithm = TVPACK())

##==== Cases where some  (lower[j], upper[j]) == (-Inf, Inf) :
S <- toeplitz(c(1, 1/2, 1/4))

set.seed(11)
P0 <- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S)
P1 <- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = TVPACK()) # had failed
P2 <- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = Miwa())
P2a<- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = Miwa(512))
P2.<- pmvnorm(lower=c(-Inf, 0, 0), upper=Inf, corr=S, algorithm = Miwa(2048))

stopifnot(chk(1/3, c(P0), tol=1e-14),
          chk(1/3, c(P1), tol=1e-14),
          chk(1/3, c(P2), tol=1e-9 ), # 3.765e-10
          chk(1/3, c(P2a),tol=4e-12), # 8.32 e-13
          chk(1/3, c(P2.),tol=2e-12) # 5.28 e-13
)

## t-dist [TVPACK() had failed] :
set.seed(11)
Ptdef <- replicate(20, c(pmvt(lower=c(-Inf, 1, 2), upper=Inf, df=2, corr=S)))
unique(Ptdef)# see length 1; i.e., same result [even though default is Monte-Carlo ??]
Pt1 <- pmvt(lower=c(-Inf, 1, 2), upper=Inf, df=2, corr=S, algorithm = TVPACK())
P. <- 0.0570404044526986
stopifnot(exprs = {
    chk(P., c(Pt1), tol = 1e-14)# seen 3.65 e-16
    abs(P. - Ptdef) < 1e-15 # seen 1.39 e-17
})

##-------- Fix dimension reduction to dimension 1  for algo  Miwa() and TVPACK(): ---------

### Default algorithm  Gentz  works fine :
## 3-D
r1 <- pmvnorm(lower= rep(-Inf,3),
              upper= c(-1, rep(Inf,2)), sigma=diag(3))
str(r1)
stopifnot(all.equal(c(r1), pnorm(-1), tolerance = 4e-16))

## TVPACK()
r2 <- pmvnorm(lower= rep(-Inf,3),
             upper= c(-1, rep(Inf,2)), sigma=diag(3), algorithm = TVPACK())
r2
stopifnot(all.equal(c(r2), pnorm(-1), tolerance = 0),
          identical("Normal Complettion (dim reduced to 1)", attr(r2, "msg")))
## Previously gave error: need n = 2 or 3 for TVPACK() algorithm

## Miwa()
r <- pmvnorm(lower= rep(-Inf,3),
             upper= c(-1, rep(Inf,2)), sigma=diag(3), algorithm = Miwa())
## gave  *** caught segfault *** ... address 0x7fff76ac80a0, cause 'memory not mapped'
r
stopifnot(all.equal(r, r2, tolerance=0))

##-------- simple 2D :  works fine with default algo ---------------
str(t1 <- pmvnorm(lower= c(-Inf,-Inf),
                  upper= c(- 1,  Inf), sigma=diag(2)))
stopifnot(all.equal(c(t1), pnorm(-1), tolerance = 4e-16))

(tT <- pmvnorm(lower= c(-Inf,-Inf),
               upper= c(- 1,  Inf), sigma=diag(2), algorithm = TVPACK()))
## gave Error   either needs all(lower == -Inf) or all(upper == Inf).

tM <- pmvnorm(lower= c(-Inf,-Inf),
               upper= c(- 1,  Inf), sigma=diag(2), algorithm = Miwa())
## -- gave --  Warning message:
## In probval.Miwa(algorithm, n, df, lower, upper, infin, corr, delta) :
##   Approximating +/-Inf by +/-1000
stopifnot(all.equal(c(tT), pnorm(-1), tolerance = 0),
          all.equal(c(tM), pnorm(-1), tolerance = 0))
