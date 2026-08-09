# BellX.R
library(nlraa)
require(ggplot2)
set.seed(1234)
x <- 1:20
y <- bell(x, 8, -0.0314, 0.000317, 13) + rnorm(length(x), 0, 0.5)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSbell(x, ymax, a, b, xc), data = dat)
fit
gfrm<- y ~ ymax * exp(a *(x - xc)^2)
gfrmx<- y ~ ymax * exp(a *(x - xc)^2 + b*(x - xc)^3)
stgauss<-c(ymax=max(y), a=-0.5/(sd(x)^2), xc=mean(x))
cat("stgauss:"); print(stgauss)
st2<-c(ymax=8, a= 0.03, xc= 13)
st3<-c(ymax=8, a= 0.03, b=0, xc= 13)
fit2 <- nls(gfrm, start=st2, data=dat, trace=TRUE)
summary(fit2)
# Use the selfStart model
xfrm <- y ~ SSbell(x, ymax, a, b, xc)
library(nlsr)
fx2mx <- nlxb(xfrm, start=st3, data=dat, lower=c(0, -1e5, 0, 0),
             upper=c(1e5, 1e5, 0, 1e5), trace=TRUE, control=list(japprox="jacentral"))
fx2mx
# Or the formula for the same model with Jacobian approximation
fx2mgc <- nlxb(gfrmx, start=st3, data=dat, lower=c(0, -1e5, 0, 0),
              upper=c(1e5, 1e5, 0, 1e5), trace=TRUE, control=list(japprox="jacentral"))
fx2mgc
# Or the formula and analytic derivatives
fx2mg <- nlxb(gfrmx, start=st3, data=dat, lower=c(0, -1e5, 0, 0),
               upper=c(1e5, 1e5, 0, 1e5), trace=TRUE)
fx2mg
# Display results together
pshort(fx2mx)
pshort(fx2mgc)
pshort(fx2mg)
