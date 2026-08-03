
### Test battery for various optimization parameters for different optimizers.
###
### ...
### 
library(maxLik)
library(tinytest)

tol <- .Machine$double.eps^(0.25)
set.seed( 123 )
# generate a variable from normally distributed random numbers
N <- 50
x <- rnorm(N, 1, 2 )

## log likelihood function
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   return( llValue )
}

# start values
startVal <- c( mu = 0, sigma = 1 )

# 
expect_silent(ml <- maxLik( llf, start = startVal ))
expect_equivalent(coef(ml), c(1.069, 1.833), tolerance=tol)
## tol
expect_silent(mlTol <- maxLik( llf, start = startVal, tol=1))
expect_equal(returnCode(mlTol), 2)
                           # tolerance limit
expect_silent(mlTolC <- maxLik(llf, start=startVal, control=list(tol=1)))
expect_equal(coef(mlTol), coef(mlTolC))
expect_equal(hessian(mlTol), hessian(mlTolC))
expect_equal(returnCode(mlTol), returnCode(mlTolC))
expect_silent(ml <- maxLik( llf, start = startVal, tol=-1))
                           # negative tol switches tol off
expect_silent(ml <- maxLik( llf, start = startVal, control=list(tol=-1)))
expect_false(returnCode(ml) == 2)
                           # should not be w/in tolerance limit
expect_error(ml <- maxLik( llf, start = startVal, tol=c(1,2)),
             pattern="'tol' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, control=list(tol=c(1,2))),
             pattern="'tol' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, tol=TRUE),
             pattern="object of class \"logical\" is not valid for slot 'tol'")
expect_error(ml <- maxLik( llf, start = startVal, control=list(tol=TRUE)),
             pattern="object of class \"logical\" is not valid for slot 'tol'")

## ----- reltol: play w/reltol, leave other tolerances at default value -----
expect_silent(mlRelTol <- maxLik( llf, start = startVal, reltol=1))
expect_equal(returnCode(mlRelTol), 8)
mlRelTolC <- maxLik(llf, start=startVal, control=list(reltol=1))
expect_equal(coef(mlRelTol), coef(mlRelTolC))
expect_silent(ml0 <- maxLik( llf, start = startVal, reltol=0))
expect_true(nIter(ml0) > nIter(mlRelTol))
                           # switching off reltol makes more iterations
expect_silent(ml1 <- maxLik( llf, start = startVal, reltol=-1))
expect_equal(nIter(ml0), nIter(ml1))
expect_error(ml <- maxLik( llf, start = startVal, reltol=c(1,2)),
             pattern="invalid class \"MaxControl\" object: 'reltol' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, control=list(reltol=c(1,2))),
             pattern="invalid class \"MaxControl\" object: 'reltol' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, reltol=TRUE),
             pattern="assignment of an object of class \"logical\" is not valid for slot 'reltol'")
expect_error(ml <- maxLik( llf, start = startVal, control=list(reltol=TRUE)),
             pattern="assignment of an object of class \"logical\" is not valid for slot 'reltol'")
## gradtol
expect_silent(mlGradtol <- maxLik( llf, start = startVal, gradtol=0.1))
expect_equal(returnCode(mlGradtol), 1)
mlGradtolC <- maxLik(llf, start=startVal, control=list(gradtol=0.1))
expect_equal(coef(mlGradtol), coef(mlGradtolC))
expect_silent(ml <- maxLik( llf, start = startVal, gradtol=-1))
expect_true(nIter(ml) > nIter(mlGradtol))
                           # switching off gradtol makes more iterations
expect_error(ml <- maxLik( llf, start = startVal, gradtol=c(1,2)),
             pattern="object: 'gradtol' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, control=list(gradtol=c(1,2))),
             pattern="object: 'gradtol' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, gradtol=TRUE),
             pattern="assignment of an object of class \"logical\" is not valid for slot 'gradtol' ")
expect_error(ml <- maxLik( llf, start = startVal, control=list(gradtol=TRUE)),
             pattern="assignment of an object of class \"logical\" is not valid for slot 'gradtol' ")
## examples with steptol, lambdatol
## qac
expect_silent(mlMarq <- maxLik( llf, start = startVal, qac="marquardt"))
expect_equal(maximType(mlMarq),
             "Newton-Raphson maximisation with Marquardt (1963) Hessian correction")
expect_silent(mlMarqC <- maxLik(llf, start=startVal, control=list(qac="marquardt")))
expect_equal(coef(mlMarq), coef(mlMarqC))
expect_error(ml <- maxLik( llf, start = startVal, qac=-1),
             pattern = "assignment of an object of class \"numeric\" is not valid for slot 'qac'")
                           # qac should be "stephalving" or "marquardt"
expect_error(ml <- maxLik( llf, start = startVal, qac=c("a", "b")),
             pattern = "invalid class \"MaxControl\" object: 'qac' must be of length 1, not 2")
expect_error(ml <- maxLik( llf, start = startVal, qac=TRUE),
             pattern = "assignment of an object of class \"logical\" is not valid for slot 'qac'")
mlMarqCl <- maxLik(llf, start = startVal,
                        control=list(qac="marquardt", lambda0=1000, lambdaStep=4))
expect_equal(coef(mlMarqCl), coef(mlMarq))
## NM: alpha, beta, gamma
expect_silent(mlNMAlpha <- maxLik(llf, start=startVal, method="nm", beta=0.8))
expect_silent(mlNMAlphaC <- maxLik(llf, start=startVal, method="nm", control=list(beta=0.8)))
expect_equal(coef(mlNMAlpha), coef(mlNMAlphaC))

## likelihood function with additional parameter
llf1 <- function( param, sigma ) {
   mu <- param
   N <- length( x )
   ll <- -0.5*N*log( 2 * pi ) - N*log( sigma ) -
      0.5*sum( ( x - mu )^2/sigma^2 )
   ll
}

## log-lik mixture
logLikMix <- function(param) {
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
   ll
}

## loglik mixture with additional parameter
logLikMixA <- function(param, rho) {
   mu1 <- param[1]
   mu2 <- param[2]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
   ll
}

## Test the following with all the main optimizers:
pl2Patterns <- c(NR = "----- Initial parameters: -----\n.*-----Iteration 1 -----",
                 BFGS = "initial  value.*final  value",
                 BFGSR = "-------- Initial parameters: -------\n.*Iteration  1")
for(method in c("NR", "BFGS", "BFGSR")) {
   ## create data in loop, we need to mess with 'x' for constraints
   N <- 100
   x <- rnorm(N, 1, 2 )
   startVal <- c(1,2)
   ## two parameters at the same time
   ## iterlim, printLevel
   expect_stdout(ml2 <- maxLik(llf, start=startVal, method=method,
                               iterlim=1, printLevel=2),
                 pattern = pl2Patterns[method])
   expect_stdout(ml2C <- maxLik(llf, start=startVal, method=method,
                                control=list(iterlim=1, printLevel=2)),
                 pattern = pl2Patterns[method])
   expect_equal(coef(ml2), coef(ml2C))
   ## what about additional parameters for the loglik function?
   expect_silent(mlsM <- maxLik(llf1, start=0, method=method, tol=1, sigma=1))
   expect_silent(mlsCM <- maxLik(llf1, start=0, method=method, control=list(tol=1), sigma=1))
   expect_equal(coef(mlsM), coef(mlsCM))
   ## And what about unused parameters?
   expect_error(maxLik(llf1, start=0, method=method, control=list(tol=1),
                       sigma=1, unusedPar=2),
                pattern = "unused argument")
   N <- 100
   ## Does this work with constraints?
   x <- c(rnorm(N, mean=-1), rnorm(N, mean=1))
   ## First test inequality constraints
   ## Inequality constraints: x + y + z < 0.5
   A <- matrix(c(-1, 0, 0,
                 0, -1, 0,
                 0, 0, 1), 3, 3, byrow=TRUE)
   B <- rep(0.5, 3)
   start <- c(0.4, 0, 0.9)
   ## analytic gradient
   if(!(method %in% c("NR", "BFGSR"))) {
      expect_silent(mix <- maxLik(logLikMix, 
                                  start=start, method=method,
                                  constraints=list(ineqA=A, ineqB=B)))
      expect_silent(mixGT <- try(maxLik(logLikMix, 
                                        start=start, method=method,
                                        constraints=list(ineqA=A, ineqB=B),
                                        tol=1)))
      expect_silent(
         mixGTC <- try(maxLik(logLikMix, 
                              start=start, method=method,
                              constraints=list(ineqA=A, ineqB=B),
                              control=list(tol=1)))
      )
      ## 2d inequality constraints: x + y < 0.5
      A2 <- matrix(c(-1, -1), 1, 2, byrow=TRUE)
      B2 <- 0.5
      start2 <- c(-0.5, 0.5)
      expect_silent(
         mixA <- maxLik(logLikMixA, 
                        start=start2, method=method,
                        constraints=list(ineqA=A2, ineqB=B2),
                        tol=1,
                        rho=0.5)
      )
      expect_silent(
         mixAC <- maxLik(logLikMixA, 
                         start=start2, method=method,
                         constraints=list(ineqA=A2, ineqB=B2),
                         control=list(tol=1),
                         rho=0.5)
      )
      expect_equal(coef(mixA), coef(mixAC))
      expect_equal(hessian(mixA), hessian(mixAC))
   }
}

### Test adding both default and user-specified parameters through control list
estimate <- function(control=NULL, ...) {
   maxLik(llf, start=c(1,1),
          control=c(list(iterlim=100), control),
          ...)
}
expect_silent(m <- estimate(control=list(iterlim=1), fixed=2))
expect_stdout(show(maxControl(m)),
              pattern = "iterlim = 1")
                           # iterlim should be 1
expect_equal(coef(m)[2], 1)
                           # sigma should be 1.000
## Does print.level overwrite 'printLevel'?
expect_silent(m <- estimate(control=list(printLevel=2, print.level=1)))
expect_stdout(show(maxControl(m)),
              pattern = "printLevel = 1")

## Does open parameters override everything?
expect_silent(m <- estimate(control=list(printLevel=2, print.level=1), print.level=0))
expect_stdout(show(maxControl(m)),
              pattern = "printLevel = 0")

### does both printLevel, print.level work for condiNumber?
expect_silent(condiNumber(hessian(m), print.level=0))
expect_silent(condiNumber(hessian(m), printLevel=0))
expect_silent(condiNumber(hessian(m), printLevel=0, print.level=1))
