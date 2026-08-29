#### qpois(), qbinom() and qnbinom() overflow much too early in the upper tail
#### ~~~~~    ~~~~~~       ~~~~~~~

if(!dev.interactive(orNone=TRUE)) pdf("qPoisBinom-ex.pdf")
.O.P. <- par(no.readonly=TRUE)

## NB:  No bug here anymore  --- fixed by this :
## ------
## r86504 | maechler | 2024-05-01 10:43:29 +0200 (Mi, 01 Mai 2024) | 1 line
##
## fix discrete inversion p |--> q; for qbinom() and more: PR#18711

## source  ~/R/D/r-devel/R/src/nmath/qDiscrete_search.h

relErrV <- sfsmisc::relErrV

e <- c(-1000, -500, -200, -100, seq(-70,-1/4, by=1/4))
lambda <- 10000
qp <- qpois(2^e, lambda=lambda, lower.tail=FALSE)

## 'cbind_no_rownames' :
cbNoRN <- function(...) {
    r <- cbind(...)
    dimnames(r)[[1]] <- rep.int("", nrow(r))
    r
}

cbNoRN(e, p=2^e, qpois=qp)
## is all good now

plot(qp ~ e, type = "b", subset = -(1:5),
     main = paste0("qpois(2^e, lambda=",lambda,
                   ") - looks smooth now"))

## p(q(.)) ~= identity (but have *discrete*ness errors) :
l2p <- log2(ppois(qp, lambda=lambda, lower.tail=FALSE))
relE <- relErrV(e, l2p)

    all.equal(e, l2p, tolerance = 0) # 0.001125
stopifnot(exprs = {
    !is.unsorted(e)
    9800 < qp; qp < 14000
    diff(qp) < 0
    all.equal(e, l2p, tolerance = 0.002)
    0 < relE ; relE < 9e-3
})


### __________ NB   no bug anymore  !!! ____________________________

## qbinom() "same" problem -- ... solved now
qBin <- qbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE)
cbNoRN(e, p=2^e, qBin)[c(1, 75:85),]

plot(qBin ~ e, type = "b", subset = -(1:5),
     main = paste0("qbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE) -- ok now"))
abline(h=100, lty=3)
d.qB <- diff(qBin)

pB <- pbinom(qBin, size=100, prob = 0.4, lower.tail=FALSE)
l2pB <- log2(pB[pB != 0])
stopifnot(exprs = {
    100 >= qBin; qBin >= 35
    -8 <= d.qB; d.qB <= 0
    max(relErrV(e[pB != 0], l2pB)) < 0.20
     all.equal( e[pB != 0], l2pB, tolerance = 0.04) # 0.0281
})

## qnbinom() "same"  --

qNB <- qnbinom(2^e, size = 100, prob = 0.4, lower.tail=FALSE)
cbNoRN(e, p=2^e, qNB)
## now also fine:
range(dqNB <- diff(qNB))
stopifnot(1949 >= qNB, qNB >= 131,
          -773 <= dqNB, dqNB <= 0)

(pNB <- pnbinom(qNB, size=100, prob = 0.4, lower.tail=FALSE))
l2pN <- log2(pNB[pNB != 0])
     all.equal(e[pNB != 0], l2pN, tolerance = 0.) # 0.00377
stopifnot(exprs = {
     all.equal( e[pNB != 0], l2pN, tolerance = 0.01)
    max(relErrV(e[pNB != 0], l2pN)) <= 0.1 # 0.0756
})



require(DPQ)
## Embarrassing forgotten FIXME (for boring boundary cases only):
M <- 2^31; pr <- 1e-9
stopifnot(exprs = {
    qbinomR (0:1, size=M, prob=pr) == c(0, M) # was 0 1
    qnbinomR(0:1, size=M, prob=pr) == c(0, Inf) # " "
    qpoisR  (0:1,      M)          == c(0, Inf)
    qbinomR (c(-Inf,0), size=M, prob=pr, log.p=TRUE) == c(0, M)
    qnbinomR(c(-Inf,0), size=M, prob=pr, log.p=TRUE) == c(0, Inf)
    qpoisR  (c(-Inf,0),      M,          log.p=TRUE) == c(0, Inf)

    ## The same for the plain R versions
    qbinom (0:1, size=M, prob=pr) == c(0, M)
    qnbinom(0:1, size=M, prob=pr) == c(0, Inf)
    qpois  (0:1,      M)          == c(0, Inf)
    qbinom (c(-Inf,0), size=M, prob=pr, log.p=TRUE) == c(0, M)
    qnbinom(c(-Inf,0), size=M, prob=pr, log.p=TRUE) == c(0, Inf)
    qpois  (c(-Inf,0),      M,          log.p=TRUE) == c(0, Inf)
})
