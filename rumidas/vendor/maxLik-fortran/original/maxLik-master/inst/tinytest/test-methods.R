## Test methods.  Note: only test if methods work in terms of dim, length, etc,
## not in terms of values here
##
## ...
## * printing summary with max.columns, max.rows
## 
if(!requireNamespace("tinytest", quietly = TRUE)) {
   message("These tests require 'tinytest' package\n")
   q("no")
}
require(sandwich)
library(maxLik)
set.seed(0)
compareTolerance = 0.001
                           # tolerance when comparing different optimizers

## Test standard methods for "lm"
x <- runif(20)
y <- x + rnorm(20)
m <- lm(y ~ x)
expect_equal(
   nObs(m), length(y),
   info = "nObs.lm must be correct"
)
expect_equal(
   stdEr(m),
   c(`(Intercept)` = 0.357862322670879, x = 0.568707094458801)
)

## Test maxControl methods:
set.seed(9)
x <- rnorm(20, sd=2)
ll1 <- function(par) dnorm(x, mean=par, sd=1, log=TRUE)
ll2 <- function(par) dnorm(x, mean=par[1], sd=par[2], log=TRUE)
for(method in c("NR", "BFGS", "BFGSR")) {
   m <- maxLik(ll2, start=c(0, 2), method=method, control=list(iterlim=1))
   expect_equal(maxValue(m), -41.35, tolerance=0.01)
   expect_true(is.vector(gradient(m)),
               info=paste0("'gradient' returns a vector for ", method))
   expect_equal(length(gradient(m)), 2, info="'gradient(m)' is of length 2")
   expect_true(is.matrix(estfun(m)), info="'estfun' returns a matrix")
   expect_equal(dim(estfun(m)), c(20,2), info="'estfun(m)' is 20x2 matrix")
   expect_stdout(
      show(maxControl(m)),
      pattern = "Adam_momentum2 = 0\\.999"
   )
}

## Test methods for non-likelihood optimization
hatf <- function(theta) exp(- theta %*% theta)
for(optimizer in c(maxNR, maxBFGSR, maxBFGS, maxNM, maxSANN, maxCG)) {
   name <- as.character(quote(optimizer))
   res <- optimizer(hatf, start=c(1,1))
   if(name %in% c("maxNR", "maxBFGS", "maxNM", "maxCG")) {
      expect_equal(coef(res), c(0,0), tol=1e-5,
                   info=paste0(name, ": result (0,0)"))
   }
   expect_equal(objectiveFn(res), hatf, info=paste0(name, ": objectiveFn correct"))
}

## Test maxLik vcov related methods
set.seed( 15 )
t <- rexp(20, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
a <- maxLik(loglik, start=1)
expect_equal(dim(vcov(a)), c(1,1), info="vcov 1D numeric correct")
expect_equal(length(stdEr(a)), 1, info="stdEr 1D numeric correct")
a <- maxLik(loglik, gradlik, hesslik, start=1)
expect_equal(dim(vcov(a)), c(1,1), info="vcov 1D analytic correct")
expect_equal(length(stdEr(a)), 1, info="stdEr 1D analytic correct")

## ---------- both individual and aggregated likelihood ----------
NOBS <- 100
x <- rnorm(NOBS, 2, 1)
## log likelihood function
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   sum(dnorm(x, mu, sigma, log=TRUE))
}
## log likelihood function (individual observations)
llfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   return( llValues )
}
startVal <- c(mu=2, sigma=1)
ml <- maxLik( llf, start = startVal)
mlInd <- maxLik( llfInd, start = startVal)

## ---------- Various summary methods ----------
## These should work and produce consistent results
expect_stdout(
   show(confint(ml)),
   pattern = "2.5 % +97.5 %\nmu +[[:digit:] .]+\n"
)
expect_stdout(
   show(glance(ml)),
   pattern = "df logLik   AIC +nobs.*1     2  -140.  284. NA"
)
expect_stdout(
   show(glance(mlInd)),
   pattern = "df logLik   AIC  nobs.*1     2  -140.  284.   100"
)
expect_stdout(
   show(tidy(ml)),
   pattern = "term.*estimate std.error statistic.*p.value"
)


### ---------- estfun, bread, sandwich ----------
expect_error( estfun( ml ) )
expect_equal(dim(estfun( mlInd )), c(NOBS, 2))
expect_equal(colnames(estfun( mlInd )), names(startVal))

expect_error(bread( ml ) )
expect_equal(dim(bread( mlInd )), c(2, 2))
expect_equal(colnames(bread( mlInd )), names(startVal))
expect_equal(rownames(bread( mlInd )), names(startVal))

expect_error(sandwich( ml ) )
expect_equal(dim(sandwich( mlInd )), c(2, 2))
expect_equal(colnames(sandwich( mlInd )), names(startVal))
expect_equal(rownames(sandwich( mlInd )), names(startVal))
