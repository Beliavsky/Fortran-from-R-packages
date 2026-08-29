library(statmod)
options(warnPartialMatchArgs=TRUE,warnPartialMatchAttr=TRUE,warnPartialMatchDollar=TRUE)

set.seed(0); u <- runif(100)

### fitNBP

y <- matrix(rnbinom(2*4,mu=4,size=1.5),2,4)
lib.size <- rep(50000,4)
group <- c(1,1,2,2)
fitNBP(y,group=group,lib.size=lib.size)

### glmgam.fit

glmgam.fit(1,1)
glmgam.fit(c(1,1),c(0,4))
glmgam.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=rchisq(5,df=1))

### glmnb.fit

y <- rnbinom(5,mu=10,size=10)
glmnb.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=y,dispersion=0.1)
glmnb.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=y,dispersion=runif(6))
glmnb.fit(X=cbind(1,c(1,1,0,0,0)),y=c(0,0,6,2,9),dispersion=0.1)
fit <- glmnb.fit(X=cbind(1,c(1,1,0,0,0)),y=c(0,0,0,0,0),dispersion=0.1)
fit$coefficients <- zapsmall(fit$coefficients,digits=15)
fit
X <- matrix(rnorm(10),5,2)
glmnb.fit(X,y=c(0,0,0,0,0),offset=rnorm(5),dispersion=0.05)

### mixedModel2

y <- rnorm(6)
x <- rnorm(6)
z <- c(1,1,2,2,3,3)
mixedModel2(y~x,random=z)

### mixedModel2Fit

y <- c(-1,1,-2,2,0.5,1.7,-0.1)
X <- matrix(1,7,1)
Z <- model.matrix(~0+factor(c(1,1,2,2,3,3,4)))
mixedModel2Fit(y,X,Z)

### qresiduals

y <- rnorm(6)
fit <- glm(y~1)
residuals(fit)
qresiduals(fit)
qresiduals(fit,dispersion=1)

if(require("MASS")) {
  fit <- glm(Days~Age,family=negative.binomial(2),data=quine)
  print(summary(qresiduals(fit)))
  options(warnPartialMatchArgs=FALSE)
  fit <- glm.nb(Days~Age,link=log,data = quine)
  options(warnPartialMatchArgs=TRUE)
  print(summary(qresiduals(fit)))
}

### gauss.quad

options(digits=10)
g <- gauss.quad(5,"legendre")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad(5,"chebyshev1")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad(5,"chebyshev2")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad(5,"hermite")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad(5,"laguerre",alpha=5)
zapsmall(data.frame(g),digits=15)
g <- gauss.quad(5,"jacobi",alpha=5,beta=1.1)
zapsmall(data.frame(g),digits=15)
g <- gauss.quad.prob(5,dist="uniform")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad.prob(5,dist="normal")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad.prob(5,dist="beta")
zapsmall(data.frame(g),digits=15)
g <- gauss.quad.prob(5,dist="gamma")
zapsmall(data.frame(g),digits=15)

### invgauss

pinvgauss(c(0,0.1,1,2.3,3.1,NA),mean=c(1,2,3,0,1,2),dispersion=0.5)
pinvgauss(c(0,0.1,1,2.3,3.1,NA),mean=c(1,2,3,0,1,2),dispersion=0.5,log.p=TRUE)
pinvgauss(c(0,0.1,1,2.3,3.1,NA),mean=c(1,2,3,0,1,2),dispersion=0.5,lower.tail=FALSE,log.p=TRUE)
pinvgauss(1,mean=c(1,2,NA))
p <- c(0,0.001,0.5,0.999,1)
qinvgauss(p,mean=1.3,dispersion=0.6)
qinvgauss(p,mean=1.3,dispersion=0.6,lower.tail=FALSE)
qinvgauss(0.5,mean=c(1,2,NA))
qinvgauss(log(p),mean=1.3,dispersion=0.6,log.p=TRUE)
qinvgauss(log(p),mean=1.3,dispersion=0.6,lower.tail=FALSE,log.p=TRUE)
rinvgauss(5,mean=c(1,NA,3,Inf,1e10),dispersion=c(2,3,NA,Inf,4))

### tweedie

tw <- tweedie(var.power=1.25, link.power=0)
tw$linkinv( matrix(u[1:10],5,2,dimnames=list(R=LETTERS[1:5],C=letters[1:2])) )

### expectedDeviance
expectedDeviance(c(0,0.4,1),family="binomial",binom.size=2)
expectedDeviance(matrix(c(0,NA,1,Inf),2,2),family="gaussian")
expectedDeviance(c(0,1,Inf),family="Gamma",gamma.shape=2)
expectedDeviance(c(1,2),family="inverse.gaussian")
expectedDeviance(c(0,1,2),family="negative.binomial",nbinom.size=2)
expectedDeviance(c(0,2,Inf),family="poisson")

### extra tests done only locally

#GKSTest <- Sys.getenv("GKSTest")
#if(GKSTest=="on") {
#print("hello")
#}
