## Bottom of page 16 "5. A Numerical Example"
## of the  (2013) paper about the Python implementation

## Edited into R code by Martin Maechler
LBound <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
UBound <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
Mean <- c(1.175, 1.19, 0.396, 1.12, 0.346, 0.679, 0.089, 0.73, 0.481, 1.08)
Cov  <- c(0.4075516, 0.0317584, 0.9063047, 0.0518392, 0.0313639, 0.194909, 0.056639, 0.0268726, 0.0440849, 0.1952847, 0.0330226, 0.0191717, 0.0300677, 0.0277735, 0.3405911, 0.0082778, 0.0093438, 0.0132274, 0.0052667, 0.0077706, 0.1598387, 0.0216594, 0.0249504, 0.0352597, 0.0137581, 0.0206784, 0.0210558, 0.6805671, 0.0133242, 0.0076104, 0.0115493, 0.0078088, 0.0073641, 0.0051869, 0.0137788, 0.9552692, 0.0343476, 0.0287487, 0.0427563, 0.0291418, 0.0254266, 0.0172374, 0.0462703, 0.0106553, 0.3168158, 0.022499, 0.0133687, 0.020573, 0.0164038, 0.0128408, 0.0072378, 0.0192609, 0.0076096, 0.0185432, 0.1107929)

(p <- length(Mean))# 10
stopifnot(p == length(LBound), p == length(UBound),
          p*(p+1)/2 == length(Cov))

require(Matrix)
options(width = 110)# so printing shows the matrix nicely

C <- matrix(0, p,p); C[ lower.tri(C, diag=TRUE) ]  <- Cov
(C <- as(C, "Matrix")) # wrong

## try again
C <- matrix(0, p,p); C[ upper.tri(C, diag=TRUE) ]  <- Cov
(C <- t(as(C, "Matrix")))# yes

##====== get data from
if(FALSE) # << maybe offline
str(a10 <- read.csv("http://www.quantresearch.info/CLA_Data.csv.txt"))
## MM: local copy in package source:
str(a10 <- read.csv("CLA_Data.csv.txt")
a10 <- as.matrix(a10)
str(mu.10 <- a10[1,])
stopifnot(a10[2,] == 0)# == Lower_bound
stopifnot(a10[3,] == 1)# == Upper_bound
str(cc <- a10[-(1:3),]) # 10 x 10
isSymmetric(cc, check.attributes=FALSE)

LT <- lower.tri(cc, diag=TRUE)
all.equal(cc[LT], C[LT])) # ==> mean relative difference: 3.59e-7
muS.10ex <- list(mu = mu.10, # including names "X1" .. "X10"
                 covar = unname(cc))

if(FALSE) ## as maintainer, did once
 save(muS.10ex, file = "~/R/Pkgs/CLA/data/muS.10ex.rda", compress = "xz")

CLA.10ex <- with(muS.10ex, CLA(mu, covar, lB=0, uB=1)) # works after 'drop = FALSE' fix in getMatrices() !!
drop0(zapsmall(CLA.10ex$weights_set))

## or to look similar as in Alexander Norrington's M.thesis, p.33
## (but has 1st row doubled, *and* last row [lambda=0] missing (!?)
t(round(CLA.10ex$weights_set, 3))
##          X1    X2    X3    X4    X5    X6 X7    X8    X9   X10
##  [1,] 0.000 1.000 0.000 0.000 0.000 0.000  0 0.000 0.000 0.000
##  [2,] 0.000 1.000 0.000 0.000 0.000 0.000  0 0.000 0.000 0.000
##  [3,] 0.649 0.351 0.000 0.000 0.000 0.000  0 0.000 0.000 0.000
##  [4,] 0.434 0.231 0.000 0.335 0.000 0.000  0 0.000 0.000 0.000
##  [5,] 0.127 0.072 0.000 0.281 0.000 0.000  0 0.000 0.000 0.520
##  [6,] 0.123 0.070 0.000 0.279 0.000 0.000  0 0.006 0.000 0.521
##  [7,] 0.087 0.050 0.000 0.224 0.000 0.174  0 0.030 0.000 0.435
##  [8,] 0.085 0.049 0.000 0.220 0.000 0.180  0 0.031 0.006 0.429
##  [9,] 0.074 0.044 0.000 0.199 0.026 0.198  0 0.033 0.028 0.398
## [10,] 0.068 0.041 0.015 0.188 0.034 0.202  0 0.034 0.034 0.383
