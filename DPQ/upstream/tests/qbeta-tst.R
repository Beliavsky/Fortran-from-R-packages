library(DPQ)

### Log-scale  {was part of my function code -- should work in R after 2014 :}
lp <- -10*(150:20)
qb <- qbeta(lp, 2,3, log.p=TRUE)
pb <- pbeta(qb, 2,3, log.p=TRUE)
all.equal(lp, pb, countEQ=TRUE) # 0.003525 ... not particularly good; now even *worse*: 0.01925502
## here's the reason : qbeta() stays at minimum > 0, but it should rather underflow
cbind(lp, qb, pb, D = lp-pb, relD = 1 - pb/lp)

require("graphics")

plot(lp, qb, xlab = quote(log(p)),
     main = "qbeta(lp, 2,3, log=TRUE)", log="y", type="l")
## in the limit, it's fine:
stopifnot(identical(0, qbeta(-Inf, 2,3, log.p=TRUE)))

## starting from subnormal
qq <- 2^(-1074:0)
pb. <- pbeta(qq, 2, 3, log.p=TRUE)
## like a straight line
plot(qq, pb., ylab="pbeta(*, 2, 3, log.p=TRUE)",
     log="x", type="l", xaxt="n")
sfsmisc::eaxis(1)
