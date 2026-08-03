### Various tests for constrained optimization
###
options(digits=4)

### -------------------- Normal mixture likelihood, no additional parameters --------------------
### param = c(rho, mean1, mean2)
### 
### X = N(mean1) w/Pr rho
### X = N(mean2) w/Pr 1-rho
### 
logLikMix <- function(param) {
   ## a single likelihood value
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
   ll <- sum(ll)
   ll
}

gradLikMix <- function(param) {
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   f1 <- dnorm(x - mu1)
   f2 <- dnorm(x - mu2)
   L <- rho*f1 + (1 - rho)*f2
   g <- matrix(0, length(x), 3)
   g[,1] <- (f1 - f2)/L
   g[,2] <- rho*(x - mu1)*f1/L
   g[,3] <- (1 - rho)*(x - mu2)*f2/L
   colSums(g)
   g
}

hessLikMix <- function(param) {
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   f1 <- dnorm(x - mu1)
   f2 <- dnorm(x - mu2)
   L <- rho*f1 + (1 - rho)*f2
   dldrho <- (f1 - f2)/L
   dldmu1 <- rho*(x - mu1)*f1/L
   dldmu2 <- (1 - rho)*(x - mu2)*f2/L
   h <- matrix(0, 3, 3)
   h[1,1] <- -sum(dldrho*(f1 - f2)/L)
   h[2,1] <- h[1,2] <- sum((x - mu1)*f1/L - dldmu1*dldrho)
   h[3,1] <- h[1,3] <- sum(-(x - mu2)*f2/L - dldmu2*dldrho)
   h[2,2] <- sum(rho*(-f1 + (x - mu1)^2*f1)/L - dldmu1^2)
   h[2,3] <- h[3,2] <- -sum(dldmu1*dldmu2)
   h[3,3] <- sum((1 - rho)*(-f2 + (x - mu2)^2*f2)/L - dldmu2^2)
   h
}

logLikMixInd <- function(param) {
   ## individual obs-wise likelihood values
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
   ll <- sum(ll)
   ll
}

gradLikMixInd <- function(param) {
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   f1 <- dnorm(x - mu1)
   f2 <- dnorm(x - mu2)
   L <- rho*f1 + (1 - rho)*f2
   g <- matrix(0, length(x), 3)
   g[,1] <- (f1 - f2)/L
   g[,2] <- rho*(x - mu1)*f1/L
   g[,3] <- (1 - rho)*(x - mu2)*f2/L
   colSums(g)
   g
}

### --------------------------
library(maxLik)
## mixed normal
set.seed(1)
N <- 100
x <- c(rnorm(N, mean=-1), rnorm(N, mean=1))

## ---------- INEQUALITY CONSTRAINTS -----------
## First test inequality constraints, numeric/analytical gradients
## Inequality constraints: rho < 0.5, mu1 < -0.1, mu2 > 0.1
A <- matrix(c(-1, 0, 0,
              0, -1, 0,
              0, 0, 1), 3, 3, byrow=TRUE)
B <- c(0.5, 0.1, 0.1)
start <- c(0.4, 0, 0.9)
ineqCon <- list(ineqA=A, ineqB=B)
## analytic gradient
cat("Inequality constraints, analytic gradient & Hessian\n")
a <- maxLik(logLikMix, grad=gradLikMix, hess=hessLikMix,
            start=start,
            constraints=ineqCon)
all.equal(coef(a), c(0.5, -1, 1), tolerance=0.1)
                           # TRUE: relative tolerance 0.045
## No analytic gradient
cat("Inequality constraints, numeric gradient & Hessian\n")
a <- maxLik(logLikMix, 
            start=start,
            constraints=ineqCon)
all.equal(coef(a), c(0.5, -1, 1), tolerance=0.1)
                           # should be close to the true values, but N is too small
## NR method with inequality constraints
try( maxLik(logLikMix, start = start, constraints = ineqCon, method = "NR" ) )
                           # Error in maxRoutine(fn = logLik, grad = grad, hess = hess, start = start,  : 
                           #   Inequality constraints not implemented for maxNR

## BHHH method with inequality constraints
try( maxLik(logLikMix, start = start, constraints = ineqCon, method = "BHHH" ) )
                           # Error in maxNR(fn = fn, grad = grad, hess = hess, start = start, finalHessian = finalHessian,  : 
                           # Inequality constraints not implemented for maxNR

## ---------- EQUALITY CONSTRAINTS -----------------
cat("Test for equality constraints mu1 + 2*mu2 = 0\n")
A <- matrix(c(0, 1, 2), 1, 3)
B <- 0
eqCon <- list( eqA = A, eqB = B )
## default, numeric gradient
mlEq <- maxLik(logLikMix, start = start, constraints = eqCon, tol=0)
                           # only rely on gradient stopping condition
all.equal(coef(mlEq), c(0.33, -1.45, 0.72), tolerance=0.01, scale=1)
## default, individual likelihood
mlEqInd <- maxLik(logLikMixInd, start = start, constraints = eqCon, tol=0)
                           # only rely on gradient stopping condition
all.equal(coef(mlEq), coef(mlEqInd), tol=1e-4)
## default, analytic gradient
mlEqG <- maxLik(logLikMix, grad=gradLikMix,
                  start = start, constraints = eqCon )
all.equal(coef(mlEq), coef(mlEqG), tolerance=1e-4)
## default, analytic gradient, individual likelihood
mlEqGInd <- maxLik(logLikMixInd, grad=gradLikMixInd,
                     start = start, constraints = eqCon )
all.equal(coef(mlEqG), coef(mlEqGInd), tolerance=1e-4)
## default, analytic Hessian
mlEqH <- maxLik(logLikMix, grad=gradLikMix, hess=hessLikMix,
                  start=start,
                  constraints=eqCon)
all.equal(coef(mlEqG), coef(mlEqH), tolerance=1e-4)


## BFGS, numeric gradient
eqBFGS <- maxLik(logLikMix, 
                 start=start, method="bfgs",
                 constraints=eqCon,
                 SUMTRho0=1)
all.equal(coef(eqBFGS), c(0.33, -1.45, 0.72), tolerance=0.01, scale=1)

## BHHH, analytic gradient (numeric does not converge?)
eqBHHH <- maxLik(logLikMix, gradLikMix,
                 start=start, method="bhhh",
                 constraints=eqCon,
                 SUMTRho0=1)
all.equal(coef(eqBFGS), coef(eqBHHH), tol=1e-4)


### ------------------ Now test additional parameters for the function ----
### similar mixture as above but rho is give as an extra parameter
### 
logLikMix2 <- function(param, rho) {
   mu1 <- param[1]
   mu2 <- param[2]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
#   ll <- sum(ll)
   ll
}

gradLikMix2 <- function(param, rho) {
   mu1 <- param[1]
   mu2 <- param[2]
   f1 <- dnorm(x - mu1)
   f2 <- dnorm(x - mu2)
   L <- rho*f1 + (1 - rho)*f2
   g <- matrix(0, length(x), 2)
   g[,1] <- rho*(x - mu1)*f1/L
   g[,2] <- (1 - rho)*(x - mu2)*f2/L
#   colSums(g)
   g
}

hessLikMix2 <- function(param, rho) {
   mu1 <- param[1]
   mu2 <- param[2]
   f1 <- dnorm(x - mu1)
   f2 <- dnorm(x - mu2)
   L <- rho*f1 + (1 - rho)*f2
   dldrho <- (f1 - f2)/L
   dldmu1 <- rho*(x - mu1)*f1/L
   dldmu2 <- (1 - rho)*(x - mu2)*f2/L
   h <- matrix(0, 2, 2)
   h[1,1] <- sum(rho*(-f1 + (x - mu1)^2*f1)/L - dldmu1^2)
   h[1,2] <- h[2,1] <- -sum(dldmu1*dldmu2)
   h[2,2] <- sum((1 - rho)*(-f2 + (x - mu2)^2*f2)/L - dldmu2^2)
   h
}

## ---------- Equality constraints & extra parameters ------------
A <- matrix(c(1, 2), 1, 2)
B <- 0
start <- c(0, 1)
## We run only a few iterations as we want to test correct handling
## of parameters, not the final value.  We also avoid any
## debug information
iterlim <- 3
cat("Test for extra parameters for the function\n")
## NR, numeric gradient
cat("Newton-Raphson, numeric gradient\n")
a <- maxLik(logLikMix2,
            start=start, method="nr",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## NR, numeric hessian
a <- maxLik(logLikMix2, gradLikMix2, 
            start=start, method="nr",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## nr, analytic hessian
a <- maxLik(logLikMix2, gradLikMix2, hessLikMix2,
            start=start, method="nr",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## BHHH
cat("BHHH, analytic gradient, numeric Hessian\n")
a <- maxLik(logLikMix2, gradLikMix2, 
            start=start, method="bhhh",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## BHHH, analytic
a <- maxLik(logLikMix2, gradLikMix2, 
            start=start, method="bhhh",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## bfgs, no analytic gradient
a <- maxLik(logLikMix2, 
            start=start, method="bfgs",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## bfgs, analytic gradient
a <- maxLik(logLikMix2, 
            start=start, method="bfgs",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## SANN, analytic gradient
a <- maxLik(logLikMix2, gradLikMix2,
            start=start, method="SANN",
            constraints=list(eqA=A, eqB=B),
            iterlim=iterlim, SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)
## NM, numeric
a <- maxLik(logLikMix2, 
            start=start, method="nm",
            constraints=list(eqA=A, eqB=B),
            iterlim=100,
                           # use more iters for NM
            SUMTRho0=1, rho=0.5)
all.equal(coef(a), c(-1.36, 0.68), tol=0.01)

## -------------------- NR, multiple constraints --------------------
f <- function(theta) exp(-theta %*% theta)
                           # test quadratic function
## constraints:
## theta1 + theta3 = 1
## theta1 + theta2 = 1
A <- matrix(c(1, 0, 1,
              1, 1, 0), 2, 3, byrow=TRUE)
B <- c(-1, -1)
cat("NR, multiple constraints\n")
a <- maxNR(f, start=c(1,1.1,2), constraints=list(eqA=A, eqB=B))
theta <- coef(a)
all.equal(c(theta[1] + theta[3], theta[1] + theta[2]), c(1,1), tolerance=1e-4)
## Error handling for equality constraints
A <- matrix(c(1, 1), 1, 2)
B <- -1
cat("Error handling: ncol(A) != lengths(start)\n")
try(a <- maxNR(f, start=c(1, 2, 3), constraints=list(eqA=A, eqB=B)))
                           # ncol(A) != length(start)
A <- matrix(c(1, 1), 1, 2)
B <- c(-1, 2)
try(a <- maxNR(f, start=c(1, 2), constraints=list(eqA=A, eqB=B)))
                           # nrow(A) != nrow(B)
##                           
## -------------- inequality constraints & extra paramters ----------------
##
## mu1 < 1
## mu2 > -1
A <- matrix(c(-1, 0,
              0,  1), 2,2, byrow=TRUE)
B <- c(1,1)
start <- c(0.8, 0.9)
##
inEGrad <- maxLik(logLikMix2, gradLikMix2,
            start=start, method="bfgs",
            constraints=list(ineqA=A, ineqB=B),
            rho=0.5)
all.equal(coef(inEGrad), c(-0.98, 1.12), tol=0.01)
##
inE <- maxLik(logLikMix2, 
            start=start, method="bfgs",
            constraints=list(ineqA=A, ineqB=B),
            rho=0.5)
all.equal(coef(inEGrad), coef(inE), tol=1e-4)
##
inENM <- maxLik(logLikMix2, gradLikMix2,
                start=start, method="nm",
                constraints=list(ineqA=A, ineqB=B),
                rho=0.5)
all.equal(coef(inEGrad), coef(inENM), tol=1e-3)
                           # this is further off than gradient-based methods
## ---------- test vector B for inequality  --------------
## mu1 < 1
## mu2 > 2
A <- matrix(c(-1, 0,
              0,  1), 2,2, byrow=TRUE)
B1 <- c(1,-2)
a <- maxLik(logLikMix2, gradLikMix2,
            start=c(0.5, 2.5), method="bfgs",
            constraints=list(ineqA=A, ineqB=B1),
            rho=0.5)
theta <- coef(a)
all.equal(c(theta[1] < 1, theta[2] > 2), c(TRUE, TRUE))
                           # components should be larger than
                           # (-1, -2)

##
## ----  ERROR HANDLING: insert wrong A and B forms ----
##
A2 <- c(-1, 0, 0, 1)
try(maxLik(logLikMix2, gradLikMix2,
           start=start, method="bfgs",
           constraints=list(ineqA=A2, ineqB=B),
           print.level=1, rho=0.5)
    )
                           # should explain that matrix needed
A2 <- matrix(c(-1, 0, 0, 1), 1, 4)
try(maxLik(logLikMix2, gradLikMix2,
           start=start, method="bfgs",
           constraints=list(ineqA=A2, ineqB=B),
           print.level=1, rho=0.5)
    )
                           # should explain that wrong matrix
                           # dimension
B2 <- 1:3
try(maxLik(logLikMix2, gradLikMix2,
           start=start, method="bfgs",
           constraints=list(ineqA=A, ineqB=B2),
           print.level=1, rho=0.5)
    )
                           # A & B do not match
cat("A & B do not match\n")
B2 <- matrix(1,2,2)
try(maxLik(logLikMix2, gradLikMix2,
           start=start, method="bfgs",
           constraints=list(ineqA=A, ineqB=B2),
           print.level=1, rho=0.5)
    )
                           # B must be a vector

## ----  fixed parameters with constrained optimization -----
## Thanks to Bob Loos for finding this error.
## Optimize 3D hat with one parameter fixed (== 2D hat).
## Add an equality constraint on that
cat("Constraints + fixed parameters\n")
hat3 <- function(param) {
   ## Hat function.  Hessian negative definite if sqrt(x^2 + y^2) < 0.5
   x <- param[1]
   y <- param[2]
   z <- param[3]
   exp(-x^2-y^2-z^2)
}
sv <- c(1,1,1)
## constraints: x + y + z >= 2.5
A <- matrix(c(x=1,y=1,z=1), 1, 3)
B <- -2.5
constraints <- list(ineqA=A, ineqB=B)
res <- maxBFGS(hat3, start=sv, constraints=constraints, fixed=3,
               iterlim=3)
all.equal(coef(res), c(0.770, 0.770, 1), tol=0.01)
