### BFGSR-related tests

## 1. Test maximization algorithm for convex regions
## 
## Optimize quadratic form t(D) %*% W %*% D with p.d. weight matrix
## (ie unbounded problems).
## All solutions should go to large values with a message about successful convergence
set.seed(0)
options(digits=4)
quadForm <- function(D) {
   C <- seq(1, N)
   return( - t(D - C) %*% W %*% ( D - C) )
}
N <- 3
                           # 3-dimensional case
## a) test quadratic function t(D) %*% D
library(maxLik)
W <- diag(N)
D <- rep(1/N, N)
res <- maxBFGSR(quadForm, start=D)
all.equal(coef(res), 1:3, tolerance=1e-4)
all.equal(gradient(res), rep(0,3), tolerance=1e-3)
all.equal(nIter(res) < 100, TRUE)
all.equal(returnCode(res) < 4, TRUE)

## Next, optimize hat function in non-concave region.  Does not work well.
hat <- function(param) {
   ## Hat function.  Hessian negative definite if sqrt(x^2 + y^2) < 0.5
   x <- param[1]
   y <- param[2]
   exp(-(x-2)^2 - (y-2)^2)
}

hatNC <- maxBFGSR(hat, start=c(1,1), tol=0, reltol=0)
all.equal(coef(hatNC), rep(2,2), tolerance=1e-4)
all.equal(gradient(hatNC), rep(0,2), tolerance=1e-3)
all.equal(nIter(hatNC) < 100, TRUE)
all.equal(returnCode(hatNC) < 4, TRUE)

## Test BFGSR with fixed parameters and equality constraints
## Optimize 3D hat with one parameter fixed (== 2D hat).
## Add an equality constraint on that
hat3 <- function(param) {
   ## Hat function.  Hessian negative definite if sqrt((x-2)^2 + (y-2)^2) < 0.5
   x <- param[1]
   y <- param[2]
   z <- param[3]
   exp(-(x-2)^2-(y-2)^2-(z-2)^2)
}
sv <- c(x=1,y=1,z=1)
## constraints: x + y + z = 8
A <- matrix(c(1,1,1), 1, 3)
B <- -8
constraints <- list(eqA=A, eqB=B)
hat3CF <- maxBFGSR(hat3, start=sv, constraints=constraints, fixed=3)
all.equal(coef(hat3CF), c(x=3.5, y=3.5, z=1), tolerance=1e-4)
all.equal(nIter(hat3CF) < 100, TRUE)
all.equal(returnCode(hat3CF) < 4, TRUE)
all.equal(sum(coef(hat3CF)), 8, tolerance=1e-4)
