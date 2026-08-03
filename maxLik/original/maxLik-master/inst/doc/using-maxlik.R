### R code from vignette source 'using-maxlik.Rnw'

###################################################
### code chunk number 1: using-maxlik.Rnw:46-48
###################################################
library(maxLik)
set.seed(6)


###################################################
### code chunk number 2: using-maxlik.Rnw:98-107
###################################################
x <- rnorm(100)  # data.  true mu = 0, sigma = 1
loglik <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   sum(dnorm(x, mean=mu, sd=sigma, log=TRUE))
}
m <- maxLik(loglik, start=c(mu=1, sigma=2))
                           # give start value somewhat off
summary(m)


###################################################
### code chunk number 3: using-maxlik.Rnw:147-148
###################################################
coef(m)


###################################################
### code chunk number 4: using-maxlik.Rnw:151-152
###################################################
stdEr(m)


###################################################
### code chunk number 5: using-maxlik.Rnw:208-212
###################################################
## create 3 variables with very different scale
X <- cbind(rnorm(100), rnorm(100, sd=1e3), rnorm(100, sd=1e7))
## note: correct coefficients are 1, 1, 1
y <- X %*% c(1,1,1) + rnorm(100)


###################################################
### code chunk number 6: using-maxlik.Rnw:224-232
###################################################
negSSE <- function(beta) {
   e <- y - X %*% beta
   -crossprod(e)
                           # note '-': we are maximizing
}
m <- maxLik(negSSE, start=c(0,0,0))
                           # give start values a bit off
summary(m, eigentol=1e-15)


###################################################
### code chunk number 7: using-maxlik.Rnw:256-259
###################################################
grad <- function(beta) {
   2*t(y - X %*% beta) %*% X
}


###################################################
### code chunk number 8: using-maxlik.Rnw:263-265
###################################################
m <- maxLik(negSSE, grad=grad, start=c(0,0,0))
summary(m, eigentol=1e-15)


###################################################
### code chunk number 9: using-maxlik.Rnw:278-281
###################################################
hess <- function(beta) {
   -2*crossprod(X)
}


###################################################
### code chunk number 10: hessianExample
###################################################
m <- maxLik(negSSE, grad=grad, hess=hess, start=c(0,0,0))
summary(m, eigentol=1e-15)


###################################################
### code chunk number 11: SSEA
###################################################
negSSEA <- function(beta) {
   ## negative SSE with attributes
   e <- y - X %*% beta  # we will re-use 'e'
   sse <- -crossprod(e)
                           # note '-': we are maximizing
   attr(sse, "gradient") <- 2*t(e) %*% X
   attr(sse, "Hessian") <- -2*crossprod(X)
   sse
}
m <- maxLik(negSSEA, start=c(0,0,0))
summary(m, eigentol=1e-15)


###################################################
### code chunk number 12: using-maxlik.Rnw:338-340
###################################################
compareDerivatives(negSSE, grad, t0=c(0,0,0))
                           # 't0' is the parameter value


###################################################
### code chunk number 13: BFGS
###################################################
m <- maxLik(loglik, start=c(mu=1, sigma=2),
            method="BFGS")
summary(m)


###################################################
### code chunk number 14: using-maxlik.Rnw:473-493
###################################################
loglik <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   N <- length(x)
   -N*log(sqrt(2*pi)) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
                           # sum over observations
}
gradlikB <- function(theta) {
   ## BHHH-compatible gradient
   mu <- theta[1]
   sigma <- theta[2]
   N <- length(x)  # number of observations
   gradient <- matrix(0, N, 2)  # gradient is matrix:
                           # N datapoints (rows), 2 components
   gradient[, 1] <- (x - mu)/sigma^2
                           # first column: derivative wrt mu
   gradient[, 2] <- -1/sigma + (x - mu)^2/sigma^3
                           # second column: derivative wrt sigma
   gradient
}


###################################################
### code chunk number 15: using-maxlik.Rnw:503-506
###################################################
m <- maxLik(loglik, gradlikB, start=c(mu=1, sigma=2),
            method="BHHH")
summary(m)


###################################################
### code chunk number 16: using-maxlik.Rnw:514-525
###################################################
loglikB <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   -log(sqrt(2*pi)) - log(sigma) - 0.5*(x - mu)^2/sigma^2
                           # no summing here
                           # also no 'N*' terms as we work by
                           # individual observations
}
m <- maxLik(loglikB, start=c(mu=1, sigma=2),
            method="BHHH")
summary(m)


###################################################
### code chunk number 17: using-maxlik.Rnw:557-561
###################################################
m <- maxLik(loglikB, start=c(mu=1, sigma=2),
            method="BHHH",
            control=list(printLevel=3, iterlim=2))
summary(m)


###################################################
### code chunk number 18: using-maxlik.Rnw:601-605
###################################################
m <- maxLik(loglikB, start=c(mu=1, sigma=2),
            method="BHHH",
            control=list(reltol=0, gradtol=0))
summary(m)


###################################################
### code chunk number 19: using-maxlik.Rnw:639-648
###################################################
loglik <- function(theta, x) {
   mu <- theta[1]
   sigma <- theta[2]
   sum(dnorm(x, mean=mu, sd=sigma, log=TRUE))
}
m <- maxLik(loglik, start=c(mu=1, sigma=2), x=x)
                           # named argument 'x' will be passed
                           # to loglik
summary(m)


###################################################
### code chunk number 20: using-maxlik.Rnw:680-689
###################################################
f <- function(theta) {
   x <- theta[1]
   y <- theta[2]
   exp(-x^2 - y^2)
                           # optimum at (0, 0)
}
m <- maxBFGS(f, start=c(1,1))
                           # give start value a bit off
summary(m)


###################################################
### code chunk number 21: using-maxlik.Rnw:710-720
###################################################
## create 3 variables, two independent, third collinear
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- x1 + x2 + rnorm(100, sd=1e-6)  # highly correlated w/x1, x2
X <- cbind(x1, x2, x3)
y <- X %*% c(1, 1, 1) + rnorm(100)
m <- maxLik(negSSEA, start=c(x1=0, x2=0, x3=0))
                           # negSSEA: negative sum of squared errors
                           # with gradient, hessian attribute
summary(m)


###################################################
### code chunk number 22: using-maxlik.Rnw:733-734
###################################################
condiNumber(X)


###################################################
### code chunk number 23: using-maxlik.Rnw:767-780
###################################################
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
X <- cbind(x1, x2, x3)
y <- X %*% c(1, 1, 1) > 0
                           # y values 1/0 linearly separated
loglik <- function(beta) {
   link <- X %*% beta
   sum(ifelse(y > 0, plogis(link, log=TRUE),
              plogis(-link, log=TRUE)))
}
m <- maxLik(loglik, start=c(x1=0, x2=0, x3=0))
summary(m)


###################################################
### code chunk number 24: using-maxlik.Rnw:784-785
###################################################
condiNumber(X)


###################################################
### code chunk number 25: using-maxlik.Rnw:789-790
###################################################
condiNumber(hessian(m))


###################################################
### code chunk number 26: using-maxlik.Rnw:815-825
###################################################
x <- rnorm(100)
loglik <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   sum(dnorm(x, mean=mu, sd=sigma, log=TRUE))
}
m <- maxLik(loglik, start=c(mu=1, sigma=1),
            fixed="sigma")
                           # fix the component named 'sigma'
summary(m)


###################################################
### code chunk number 27: using-maxlik.Rnw:863-874
###################################################
f <- function(theta) {
   x <- theta[1]
   y <- theta[2]
   exp(-x^2 - y^2)
                           # optimum at (0, 0)
}
A <- matrix(c(1, 1), ncol=2)
B <- -1
m <- maxNR(f, start=c(1,1),
           constraints=list(eqA=A, eqB=B))
summary(m)


###################################################
### code chunk number 28: using-maxlik.Rnw:908-913
###################################################
A <- matrix(c(1, 1), ncol=2)
B <- -1
m <- maxBFGS(f, start=c(1,1),
             constraints=list(ineqA=A, ineqB=B))
summary(m)


###################################################
### code chunk number 29: using-maxlik.Rnw:947-952
###################################################
A <- matrix(c(1, 1, 1, -1), ncol=2)
B <- c(-1, -1)
m <- maxBFGS(f, start=c(2, 0),
             constraints=list(ineqA=A, ineqB=B))
summary(m)


