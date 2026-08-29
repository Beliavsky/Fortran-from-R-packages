library(DPQ)

### Some "manual" tests for ../R/pnbeta-approx2.R
###                         ~~~~~~~~~~~~~~~~~~~~~

## Some values from  Table 2 and 3 (p.152) of the paper:
## T.1
pnbetaAppr2(0.868,  10, 20, ncp=150)# 0.937569774351
## T.2
pnbetaAppr2(0.644,  10,  5, ncp= 20)# 0.0497384
pnbetaAppr2(0.80,   30, 30, ncp=140)# 0.7816822

## more extreme cases
x <- c(0.9, 0.95, 0.96, 0.97, 0.99)
pnbetaAppr2(x, 10, 20, 1000)# 7.953858e-08 8.276896e-02 3.709090e-01 8.216320e-01 9.999989e-01
pbeta       (x, 10, 20, 1000)# 6.959528e-08 8.289800e-02 3.709609e-01 8.214000e-01 9.999992e-01
## in the tail, they differ very considerably; in the center (P in (1e-7, 1-1e-7), things are ok:
pnbetaAppr2(x, 10, 30, 5000)
pbeta       (x, 10, 30, 5000)

## really large lambda: now "traditional" pbeta() seems hopeless
x <- 1 - 2^(-5:-15)
px <- cbind(appr2 = pnbetaAppr2(x, 10, 30, 1e6, log.p=TRUE),
            R     = pbeta      (x, 10, 30, 1e6, log.p=TRUE))
cbind(x, px)
# either '-Inf'; or full warnings about precision
## now (TOMS algo): much better and basically correct (but still -Inf !!)
matplot(x, px, type ="b")
