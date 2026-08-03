### general optimization tests for the functions of various forms
### test for:
### 1. numeric gradient, Hessian
### 2. analytic gradient, numeric Hessian
### 3. analytic gradient, Hessian
###
### a) maxLik(, method="NR")
### c) maxLik(, method="BFGS")
### b) maxLik(, method="BHHH")
###
### i) maxNR()
### ii) maxBFGS()
if(!requireNamespace("tinytest", quietly = TRUE)) {
   cat("These tests require 'tinytest' package\n")
   q("no")
}
library(maxLik)

## ---------- define log-likelihood functions ----------
## log-likelihood function(s)
logLL <- function(x, X)   # per observation for maxLik
   dgamma(x = X, shape = x[1], scale = x[2], log = TRUE)
logLLSum <- function(x, X)
   sum(logLL(x, X))

# gradient of log-likelihood function
d.logLL <- function(x, X){   # analytic 1. derivatives
   shape <- x[1]
   scale <- x[2]
   cbind(shape= log(X) - log(scale) - psigamma(shape, 0),
         scale= (X/scale - shape)/scale
         )
}
d.logLLSum <- function(x, X) {
   ## analytic 1. derivatives, summed
   colSums(d.logLL(x, X))
}

## Hessian of log-likelihood function
dd.logLL <- function(x, X){   # analytic 2. derivatives
   shape <- x[1]
   scale <- x[2]
   hessian <- matrix(0, 2, 2)
   hessian[1,1] <- -psigamma(shape, 1)*length(X)
   hessian[2,2] <- (shape*length(X) - 2*sum(X)/scale)/scale^2
   hessian[cbind(c(2,1), c(1,2))] <- -length(X)/scale
   return(hessian)
}

## ---------- create data ----------
## sample size 1000 should give precision 0.1 or better
param <- c(1.5, 2)
set.seed(100)
testData <- rgamma(1000, shape=param[1], scale=param[2])
start <- c(1,1)
mTol <- .Machine$double.eps^0.25

## estimation with maxLik() / NR
doTests <- function(method="NR") {
   suppressWarnings(rLLSum <- maxLik( logLLSum, start=start, method=method, X=testData ))
   stdDev <- stdEr(rLLSum)
   tol <- 2*max(stdDev)
   expect_equal(coef(rLLSum), param, tolerance=tol,
                info=paste("coefficient values should be close to the true values", paste(param, collapse=", ")))
                           # should equal to param, but as N is small, it may be way off
   ##
   rLL <- suppressWarnings(maxLik( logLL, start = start, method=method, X=testData ))
   expect_equal(coef(rLL), coef(rLLSum), tolerance=mTol)
   ##
   rLLSumGSum <- suppressWarnings(maxLik( logLLSum, grad=d.logLLSum, start = start, method=method, X=testData ))
   expect_equal(coef(rLLSumGSum), coef(rLLSum), tolerance=mTol)
   rLLG <- suppressWarnings(maxLik( logLL, grad=d.logLL, start = start, method=method, X=testData ))
   expect_equal(coef(rLLG), coef(rLLSum), tolerance=mTol)
   rLLGH <- suppressWarnings(maxLik( logLL, grad=d.logLL, hess=dd.logLL, start = start, method=method, X=testData ))
   expect_equal(coef(rLLGH), coef(rLLSum), tolerance=mTol)
}

doTests("NR")
doTests("BFGS")
## maxBHHH: cannot run the same tests
method <- "BHHH"
expect_error(
   maxLik( logLLSum, start=start, method=method, X=testData),
   pattern = "not provided by .* returns a numeric vector"
)
rLL <- suppressWarnings(maxLik( logLL, start = start, method=method, X=testData ))
stdDev <- stdEr(rLL)
tol <- 2*max(stdDev)
expect_equal(coef(rLL), param, tolerance=tol,
             info=paste("coefficient values should be close to the true values", paste(param, collapse=", ")))
                           # should equal to param, but as N is small, it may be way off
##
rLLG <- suppressWarnings(maxLik( logLL, grad=d.logLL, start = start, method=method, X=testData ))
expect_equal(coef(rLLG), coef(rLL), tolerance=mTol)

## Do the other basic functions work?
expect_equal(class(logLik(rLL)), "numeric")
expect_equal(class(gradient(rLL)), "numeric")
expect_true(inherits(hessian(rLL), "matrix"),
            info="Hessian must inherit from matrix class")

## test maxNR with gradient and hessian as attributes
W <- matrix(-c(4,1,2,4), 2, 2)
c <- c(1,2)
start <- c(0,0)
f <- function(x) {
   hess <- 2*W
   grad <- 2*W %*% (x - c)
   val <- t(x - c) %*% W %*% (x - c)
   attr(val, "gradient") <- as.vector(grad)
                           # gradient matrices only work for BHHH-type problems
   attr(val, "hessian") <- hess
   val
}
res <- maxNR(f, start=start)
expect_equal(coef(res), c, tolerance=mTol)
expect_equal(sqrt(sum(gradient(res)^2)), 0, tolerance=mTol)
expect_equal(maxValue(res), 0, tolerance=mTol)
