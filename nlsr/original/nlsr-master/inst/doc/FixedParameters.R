## ----opx----------------------------------------------------------------------
require(optimx)
sq<-function(x){
   nn<-length(x)
   yy<-1:nn
   f<-sum((yy-x)^2)
   f
}
sq.g <- function(x){
   nn<-length(x)
   yy<-1:nn
   gg<- 2*(x - yy)
}
xx <- c(.3, 4)
uncans <- Rvmmin(xx, sq, sq.g)
proptimr(uncans)
mybm <- c(0,1) # fix parameter 1
cans <- Rvmmin(xx, sq, sq.g, bdmsk=mybm)
proptimr(cans)
require(nlsr)
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
   38.558, 50.156, 62.948, 75.995, 91.972)
ii <- 1:12
wdf <- data.frame(weed, ii)
weedux <- nlxb(weed~b1/(1+b2*exp(-b3*ii)), start=c(b1=200, b2=50, b3=0.3)) 
weedux
weedcx <- nlxb(weed~b1/(1+b2*exp(-b3*ii)), start=c(b1=200, b2=50, b3=0.3), masked=c("b1")) 
weedcx
rfn <- function(bvec, weed=weed, ii=ii){
  res <- rep(NA, length(ii))
  for (i in ii){
    res[i]<- bvec[1]/(1+bvec[2]*exp(-bvec[3]*i))-weed[i]
  }
  res
}
weeduf <- nlfb(start=c(200, 50, 0.3),resfn=rfn,weed=weed, ii=ii,
               control=list(japprox="jacentral"))
weeduf
weedcf <- nlfb(start=c(200, 50, 0.3),resfn=rfn,weed=weed, ii=ii, lower=c(200, 0, 0),
               upper=c(200, 100,100), control=list(japprox="jacentral"))
weedcf

## ----bellx, code=readLines("./BellX.R"), echo=TRUE----------------------------
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

