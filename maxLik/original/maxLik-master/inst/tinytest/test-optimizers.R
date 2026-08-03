### This code tests all the methods and main parameters.  It includes:
### * analytic gradients/Hessian
### * fixed parameters
### * inequality constraints
### * equality constraints

## do not run unless 'NOT_CRAN' explicitly defined
## (Suggested by Sebastian Meyer and others)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    message("skipping slow optimizer tests")
    q("no")
}
if(!requireNamespace("tinytest", quietly = TRUE)) {
   message("These tests require 'tinytest' package\n")
   q("no")
}
library(maxLik)

## data to fit a normal distribution
# set seed for pseudo random numbers
set.seed( 123 )
tol <- .Machine$double.eps^0.25
## generate a variable from normally distributed random numbers
truePar <- c(mu=1, sigma=2)
NOBS <- 100
x <- rnorm(NOBS, truePar[1], truePar[2] )
xSaved <- x

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

## function to calculate analytical gradients
gf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llGrad <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llGrad )
}

## function to calculate analytical gradients (individual observations)
gfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llGrads <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llGrads )
}

## log likelihood function with gradients as attributes
llfGrad <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   attributes( llValue )$gradient <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llValue )
}

## log likelihood function with gradients as attributes (individual observations)
llfGradInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   attributes( llValues )$gradient <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llValues )
}

## function to calculate analytical Hessians
hf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llHess <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llHess )
}

## log likelihood function with gradients and Hessian as attributes
llfGradHess <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   attributes( llValue )$gradient <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   attributes( llValue )$hessian <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llValue )
}

## log likelihood function with gradients as attributes (individual observations)
llfGradHessInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   attributes( llValues )$gradient <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   attributes( llValues )$hessian <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llValues )
}

# start values
startVal <- c( mu = 0, sigma = 1 )

## basic NR: test if all methods work
ml <- maxLik( llf, start = startVal )
expect_equal(
   coef(ml), truePar, tol=2*max(stdEr(ml))
)
expect_stdout(
   print( ml ),
   pattern = "Estimate\\(s\\): 1.18.*1.81"
)
expect_stdout(
   print( summary( ml )),
   pattern = "Estimates:"
)
expect_equal(
   activePar( ml ), c(mu=TRUE, sigma=TRUE)
)
expect_equal(
   AIC( ml ), 407.167892384587,
   tol = 0.1, check.attributes=FALSE
)
expect_equal(
   coef( ml ), c(mu=1.181, sigma=1.816),
   tol = 0.001
)
expect_stdout(
   condiNumber( ml, digits = 3),
   "mu[[:space:]]+1[[:space:]\n]+sigma[[:space:]]+1\\."
)
expect_equal(
   hessian( ml), matrix(c(-30.3, 0, 0, -60.6), 2, 2),
   tol = 0.01, check.attributes = FALSE
)
expect_equal(
   logLik( ml ), -201.583946192294,
   tol = tol, check.attributes = FALSE
)
expect_equal(
   maximType( ml ), "Newton-Raphson maximisation"
)
expect_equal(
   nIter( ml ) > 5, TRUE
)
expect_error(
   nObs( ml ),
   "cannot return the number of observations"
)
expect_equal(
   nParam( ml ), 2
)
expect_equal(
   returnCode( ml ), 1
)
expect_equal(
   returnMessage( ml ), "gradient close to zero (gradtol)"
)
expect_equal(
   vcov( ml ), matrix(c(0.032975, 0, 0, 0.0165), 2, 2),
   tol=0.01, check.attributes = FALSE
)
expect_equal(
   logLik( summary( ml ) ), logLik(ml)
)
mlInd <- maxLik( llfInd, start = startVal )
expect_stdout(
   print( summary( mlInd ), digits = 2 ),
   "mu +1\\.18"
)
expect_equal(
   nObs( mlInd ), length(x)
)
## Marquardt (1963) correction
mlM <- maxLik( llf, start = startVal, qac="marquardt")
expect_equal(
   coef(mlM), coef(ml),
                           # coefficients should be the same as above
   tol=tol
)
expect_equal(
   returnMessage(mlM), returnMessage(ml)
)

## test plain results with analytical gradients
## compare coefficients, Hessian
mlg <- maxLik(llf, gf, start = startVal )
expect_equal(coef(ml), coef(mlg), tol=tol)
expect_equal(hessian(ml), hessian(mlg), tolerance = 1e-2)
## gradient with individual components
mlgInd <- maxLik( llfInd, gfInd, start = startVal )
expect_equal(coef(mlInd), coef(mlgInd), tolerance = 1e-3)
expect_equal(hessian(mlg), hessian(mlgInd), tolerance = 1e-3)

## with analytical gradients as attribute
mlG <- maxLik( llfGrad, start = startVal )
expect_equal(coef(mlG), coef(mlg), tolerance = tol)
expect_equivalent(gradient(mlG), gf( coef( mlG ) ), tolerance = tol)
mlGInd <- maxLik( llfGradInd, start = startVal )
expect_equal(coef(mlGInd), coef(mlgInd), tolerance = tol)
expect_equivalent(gradient(mlGInd), colSums( gfInd( coef( mlGInd ) ) ), tolerance = tol)
expect_equivalent(estfun(mlGInd), gfInd( coef( mlGInd ) ), tolerance=tol)

## with analytical gradients as argument and attribute
expect_warning(mlgG <- maxLik( llfGrad, gf, start = startVal))
expect_equal(coef(mlgG), coef(mlg), tolerance = tol)

## with analytical gradients and Hessians
mlgh <- maxLik( llf, gf, hf, start = startVal )
expect_equal(coef(mlg), coef(mlgh), tolerance = tol)

## with analytical gradients and Hessian as attribute
mlGH <- maxLik( llfGradHess, start = startVal )
expect_equal(coef(mlGH), coef(mlgh), tolerance = tol)

## with analytical gradients and Hessian as argument and attribute
expect_warning(mlgGhH <- maxLik( llfGradHess, gf, hf, start = startVal ))
expect_equal(coef(mlgGhH), coef(mlgh), tolerance = tol)


## ---------- BHHH method ----------
## cannot do BHHH if llf not provided by individual
x <- xSaved[1]
expect_error( maxLik( llfInd, start = startVal, method = "BHHH" ) )
## 2 observations: can do BHHH
x <- xSaved[1:2]
expect_silent( maxLik( llfInd, start = startVal, method = "BHHH" ) )
##
x <- xSaved
mlBHHH <- maxLik( llfInd, start = startVal, method = "BHHH" )
expect_stdout(print( mlBHHH ),
              pattern = "Estimate\\(s\\): 1\\.18.* 1\\.81")
expect_stdout(print(summary( mlBHHH)), pattern = "mu *1.18")
expect_equivalent(activePar( mlBHHH ), c(TRUE, TRUE))
expect_equivalent(AIC( mlBHHH ), 407.168, tolerance=0.01)
expect_equal(coef( mlBHHH ), setNames(c(1.180808, 1.816485), c("mu", "sigma")), tolerance=tol)
expect_equal(condiNumber( mlBHHH, printLevel=0),
             setNames(c(1, 1.72), c("mu", "sigma")), tol=0.01)
expect_equivalent(hessian( mlBHHH ),
                  matrix(c(-30.306411, -1.833632, -1.833632, -55.731646), 2, 2),
                  tolerance=0.01)
expect_equivalent(logLik( mlBHHH ), -201.583946192983, tolerance=tol)
expect_equal(maximType( mlBHHH ), "BHHH maximisation")
expect_equal(nIter(mlBHHH) > 3, TRUE)
                           # here 12 iterations
expect_equal(nParam( mlBHHH ), 2)
expect_equal(returnCode( mlBHHH ), 8)
expect_equal(returnMessage( mlBHHH ),
             "successive function values within relative tolerance limit (reltol)")
expect_equivalent(vcov( mlBHHH ),
                  matrix(c(0.03306213, -0.00108778, -0.00108778, 0.01797892), 2, 2),
                  tol=0.001)
expect_equivalent(logLik(summary(mlBHHH)), -201.583946192983, tolerance=tol)
expect_equal(coef(ml), coef(mlBHHH), tol=tol)
expect_equal(stdEr(ml), stdEr(mlBHHH), tol=0.1)
expect_equal(nObs( mlBHHH ), length(x))
# final Hessian = usual Hessian
expect_silent(mlBhhhH <- maxLik( llfInd, start = startVal, method = "BHHH", 
                                finalHessian = TRUE )
              )
                           # do not test Hessian equality--BHHH may be imprecise, at least
                           # for diagonal elements
expect_stdout(print(hessian( mlBhhhH )),
              pattern="mu.*\nsigma.+")
## Marquardt (1963) correction
expect_silent(mlBHHHM <- maxLik( llfInd, start = startVal, method = "BHHH", qac="marquardt"))
expect_equal(coef(mlBHHHM), coef(mlBHHH), tolerance=tol)
expect_equal(returnMessage(mlBHHHM), "successive function values within relative tolerance limit (reltol)")

## BHHH with analytical gradients
expect_error( maxLik( llf, gf, start = startVal, method = "BHHH" ) )
                           # need individual log-likelihood
expect_error( maxLik( llfInd, gf, start = startVal, method = "BHHH" ) )
                           # need individual gradient
x <- xSaved[1]  # test with a single observation
expect_error(maxLik( llf, gfInd, start = startVal, method = "BHHH" ))
                           # gradient must have >= 2 rows
expect_error( maxLik( llfInd, gfInd, start = startVal, method = "BHHH" ) )
                           # ditto even if individual likelihood components
x <- xSaved[1:2]  # test with 2 observations
expect_silent(maxLik( llf, gfInd, start = startVal, method = "BHHH",
                     iterlim=1))
                           # should work with 2 obs
expect_silent( maxLik( llfInd, gfInd, start = startVal, method = "BHHH",
                      iterlim=1) )
                           # should work with 2 obs
x <- xSaved
expect_silent(mlgBHHH <- maxLik( llfInd, gfInd, start = startVal, method = "BHHH" ))
                           # individual log-likelihood, gradient
expect_equal(coef(mlBHHH), coef(mlgBHHH), tolerance = tol)
expect_equal(coef(mlg), coef(mlgBHHH), tolerance = tol)
expect_silent(mlgBHHH2 <- maxLik( llf, gfInd, start = startVal, method = "BHHH" ))
                           # aggregated log-likelihood, individual gradient
expect_equal(coef(mlgBHHH), coef(mlgBHHH2), tolerance=tol)
                           # final Hessian = usual Hessian
expect_silent(
   mlgBhhhH <- maxLik( llf, gfInd, start = startVal, method = "BHHH", 
                      finalHessian = TRUE )
)
expect_equal(hessian(mlgBhhhH), hessian(mlBhhhH), tolerance = 1e-2)

## with analytical gradients as attribute
expect_error( maxLik( llfGrad, start = startVal, method = "BHHH" ) )
                           # no individual gradients provided
x <- xSaved[1]
expect_error( maxLik( llfGrad, start = startVal, method = "BHHH" ),
             pattern = "gradient is not a matrix")
                           # get an error about need a matrix
expect_error( maxLik( llfGradInd, start = startVal, method = "BHHH" ),
             pattern = "at least as many rows")
                           # need at least two obs
x <- xSaved[1:2]
expect_error( maxLik( llfGrad, start = startVal, method = "BHHH" ),
             pattern = "gradient is not a matrix")
                           # enough obs but no individual grad
x <- xSaved
expect_silent(mlGBHHH <- maxLik( llfGradInd, start = startVal, method = "BHHH" ))
expect_equal(coef(mlGBHHH), coef(mlgBHHH), tolerance = tol)
                           # final Hessian = usual Hessian
expect_silent(mlGBhhhH <- maxLik( llfGradInd, start = startVal, method = "BHHH", 
                                 finalHessian = TRUE ))
expect_equal(hessian(mlGBhhhH), hessian(mlgBhhhH), tolerance = tol)

## with analytical gradients as argument and attribute
expect_warning(mlgGBHHH <- maxLik( llfGradInd, gfInd, start = startVal, method = "BHHH" ),
               pattern = "both as attribute 'gradient' and as argument 'grad'")
                           # warn about double gradient
expect_equal(coef(mlgGBHHH), coef(mlgBHHH), tolerance = tol)
## with unused Hessian
expect_silent(mlghBHHH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH" ))
expect_equal(coef(mlgBHHH), coef(mlghBHHH), tolerance = tol)
## final Hessian = usual Hessian
expect_silent(
   mlghBhhhH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH", 
                       finalHessian = TRUE )
)
expect_equivalent(hessian(mlghBhhhH), hessian(mlghBHHH), tolerance = 0.2)
                           # BHHH and ordinary hessian differ quite a bit
## with unused Hessian as attribute
expect_silent(mlGHBHHH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH" ))
expect_equal(coef(mlGHBHHH), coef(mlghBHHH), tolerance = tol)
## final Hessian = usual Hessian
expect_silent(mlGHBhhhH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH", 
                                  finalHessian = TRUE ))
expect_equal(hessian(mlGHBhhhH), hessian(mlghBhhhH), tolerance = tol)
## with analytical gradients and Hessian as argument and attribute
expect_warning(
   mlgGhHBHHH <- maxLik( llfGradHessInd, gfInd, hf, start = startVal, method = "BHHH" ),
   pattern = "both as attribute 'gradient' and as argument 'grad': ignoring"
)
expect_equal(coef(mlgGhHBHHH), coef(mlghBHHH), tolerance = tol)
expect_equal(hessian(mlgGhHBHHH), hessian(mlGHBHHH), tolerance = tol)

## ---------- Test BFGS methods ----------
optimizerNames <- c(bfgsr = "BFGSR", bfgs = "BFGS", nm = "Nelder-Mead",
                    sann = "SANN", cg = "CG")
successCodes <- list(bfgsr = 1:4, bfgs = 0, nm = 0, sann = 0, cg = 0)
successMsgs <- list(bfgsr = c("successive function values within tolerance limit (tol)"),
                    bfgs = c("successful convergence "),
                           # includes space at end...
                    nm = c("successful convergence "),
                    sann = c("successful convergence "),
                    cg = c("successful convergence ")
                    )
for(optimizer in c("bfgsr", "bfgs", "nm", "sann", "cg")) {
   expect_silent(mlResult <- maxLik( llf, start = startVal, method = optimizer ))
   expect_stdout(print( mlResult ),
                 pattern = paste0(optimizerNames[optimizer], " maximization")
                 )
   expect_stdout(print( summary( mlResult )),
                 pattern = paste0(optimizerNames[optimizer], " maximization,.*Estimates:")
                 )
   expect_equal(coef(ml), coef(mlResult), tolerance=0.001)
   expect_equal(stdEr(ml), stdEr(mlResult), tolerance=0.01)
   expect_equal(activePar( mlResult ), c(mu=TRUE, sigma=TRUE))
   expect_equivalent(AIC( mlResult ), 407.167893392749, tolerance=tol)
   expect_equivalent( hessian( mlResult ),
                     matrix(c(-30.32596, 0.00000, 0.00000, -60.59508), 2, 2),
                     tolerance = 0.01)
   expect_equivalent(logLik( mlResult ), -201.5839, tolerance = 0.01)
   expect_equal(maximType( mlResult ),
                paste0(optimizerNames[optimizer], " maximization")
                )
   expect_true(nIter( mlResult ) > 1 & is.integer(nIter(mlResult)))
   expect_error( nObs( mlResult ),
                pattern = "cannot return the number of observations")
   expect_equal(nParam( mlResult ), 2)
   expect_true(returnCode( mlResult ) %in% successCodes[[optimizer]])
   expect_equal(returnMessage( mlResult), successMsgs[[optimizer]])
   expect_equal(logLik( summary( mlResult ) ), logLik(mlResult))
   ## individual observations
   expect_silent(mlIndResult <- maxLik( llfInd, start = startVal, method = optimizer))
   expect_stdout(print( summary( mlIndResult )),
                 pattern = paste0(optimizerNames[optimizer], " maximization,.*Estimates:")
                 )
   expect_equal(coef(mlResult), coef(mlIndResult), tolerance = tol)
   expect_equal(stdEr(mlResult), stdEr(mlIndResult), tolerance = 0.01)
   expect_equal(nObs( mlIndResult ), length(x))
   ## with analytic gradients
   expect_silent(mlgResult <- maxLik( llf, gf, start = startVal, method = optimizer))
   expect_equal(coef(mlgResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgResult), stdEr(mlResult), tolerance = 0.01)
   expect_silent(mlgIndResult <- maxLik( llfInd, gfInd, start = startVal,
                                        method = optimizer ))
   expect_equal(coef(mlgIndResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgIndResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients as attribute
   expect_silent(mlGResult <- maxLik( llfGrad, start = startVal,
                                     method = optimizer))
   expect_equal(coef(mlGResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlGResult), stdEr(mlResult), tolerance = 0.01)
   expect_silent(mlGIndResult <- maxLik( llfGradInd, start = startVal, method = optimizer ))
   expect_equal(coef(mlGIndResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlGIndResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients as argument and attribute
   expect_warning(mlgGResult <- maxLik( llfGrad, gf, start = startVal, method = optimizer ))
   expect_equal(coef(mlgGResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgGResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients and Hessians
   expect_silent(mlghResult <- maxLik( llf, gf, hf, start = startVal, method = optimizer ))
   expect_equal(coef(mlghResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlghResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients and Hessian as attribute
   expect_silent(mlGHResult <- maxLik( llfGradHess, start = startVal, method = optimizer ))
   expect_equal(coef(mlGHResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlGHResult), stdEr(mlResult), tolerance = 0.01)
   ## with analytical gradients and Hessian as argument and attribute
   expect_warning(mlgGhHResult <- maxLik( llfGradHess, gf, hf, start = startVal, method = optimizer ))
   expect_equal(coef(mlgGhHResult), coef(mlResult), tolerance = tol)
   expect_equal(stdEr(mlgGhHResult), stdEr(mlResult), tolerance = 0.01)
}


### ---------- with fixed parameters ----------
## start values
startValFix <- c( mu = 1, sigma = 1 )
## fix mu (the mean ) at its start value
isFixed <- c( TRUE, FALSE )
successMsgs <- list(bfgsr = c("successive function values within tolerance limit (tol)"),
                    bfgs = c("successful convergence "),
                           # includes space at end...
                    nm = c("successful convergence "),
                    sann = c("successful convergence "),
                    cg = c("successful convergence ")
                    )
## NR method with fixed parameters
for(optimizer in c("nr", "bfgsr", "bfgs", "sann", "cg")) {
   expect_silent(
      mlFix <- maxLik( llf, start = startValFix, fixed = isFixed, method=optimizer)
   )
   expect_equivalent(coef(mlFix)[1], 1)
   expect_equivalent(stdEr(mlFix)[1], 0)
   expect_silent(
      mlFix3 <- maxLik(llf, start = startValFix, fixed = "mu", method=optimizer)
   )
   expect_equal(coef(mlFix), coef(mlFix3))
   mlFix4 <- maxLik( llf, start = startValFix, fixed = which(isFixed),
                    method=optimizer)
   expect_equal(coef(mlFix), coef(mlFix4), tolerance=tol)
   expect_equivalent(activePar( mlFix ), !isFixed)
   expect_equal(nParam( mlFix ), 2)
   ## with analytical gradients
   mlgFix <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
                    method=optimizer)
   expect_equal(coef(mlgFix), coef(mlFix), tolerance=tol)
   ## with analytical gradients and Hessians
   mlghFix <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
                     method=optimizer)
   expect_equal(coef(mlghFix), coef(mlFix), tolerance=tol)
}
## Repeat the previous for NM as that one does not like 1-D optimization
for(optimizer in c("nm")) {
   expect_warning(
      mlFix <- maxLik( llf, start = startValFix, fixed = isFixed, method=optimizer)
   )
   expect_equivalent(coef(mlFix)[1], 1)
   expect_equivalent(stdEr(mlFix)[1], 0)
   expect_warning(
      mlFix3 <- maxLik(llf, start = startValFix, fixed = "mu", method=optimizer)
   )
   expect_equal(coef(mlFix), coef(mlFix3))
   expect_warning(
      mlFix4 <- maxLik( llf, start = startValFix, fixed = which(isFixed),
                       method=optimizer)
   )
   expect_equal(coef(mlFix), coef(mlFix4), tolerance=tol)
   expect_equivalent(activePar( mlFix ), !isFixed)
   expect_equal(nParam( mlFix ), 2)
   ## with analytical gradients
   expect_warning(
      mlgFix <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
             method=optimizer)
   )
   expect_equal(coef(mlgFix), coef(mlFix), tolerance=tol)
   ## with analytical gradients and Hessians
   expect_warning(
      mlghFix <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
             method=optimizer)
   )
   expect_equal(coef(mlghFix), coef(mlFix), tolerance=tol)
}
## Repeat for BHHH as that one need a different log-likelihood function
for(optimizer in c("bhhh")) {
   expect_silent(
      mlFix <- maxLik( llfInd, start = startValFix, fixed = isFixed, method=optimizer)
   )
   expect_equivalent(coef(mlFix)[1], 1)
   expect_equivalent(stdEr(mlFix)[1], 0)
   expect_silent(
      mlFix3 <- maxLik(llfInd, start = startValFix, fixed = "mu", method=optimizer)
   )
   expect_equal(coef(mlFix), coef(mlFix3))
   expect_silent(
      mlFix4 <- maxLik( llfInd, start = startValFix, fixed = which(isFixed),
                       method=optimizer)
   )
   expect_equal(coef(mlFix), coef(mlFix4), tolerance=tol)
   expect_equivalent(activePar( mlFix ), !isFixed)
   expect_equal(nParam( mlFix ), 2)
   ## with analytical gradients
   expect_silent(
      mlgFix <- maxLik( llf, gfInd, start = startValFix, fixed = isFixed,
             method=optimizer)
   )
   expect_equal(coef(mlgFix), coef(mlFix), tolerance=tol)
   ## with analytical gradients and Hessians
   expect_silent(
      mlghFix <- maxLik( llf, gfInd, hf, start = startValFix, fixed = isFixed,
             method=optimizer)
   )
   expect_equal(coef(mlghFix), coef(mlFix), tolerance=tol)
}

### ---------- inequality constraints ----------
A <- matrix( -1, nrow = 1, ncol = 2 )
inEq <- list( ineqA = A, ineqB = 2.5 )
                           # A theta + B > 0 i.e.
                           # mu + sigma < 2.5
for(optimizer in c("bfgs", "nm", "sann")) {
   expect_silent(
      mlInEq <- maxLik( llf, start = startVal, constraints = inEq,
                       method = optimizer )
   )
   expect_stdout(
      print( summary( mlInEq)),
      pattern = "constrained likelihood estimation. Inference is probably wrong.*outer iterations, barrier value"
   )
   expect_true(sum(coef( mlInEq )) < 2.5)
}

### ---------- equality constraints ----------
eqCon <- list(eqA = A, eqB = 2.5)
                           # A theta + B = 0 i.e.
                           # mu + sigma = 2.5
for(optimizer in c("nr", "bhhh", "bfgs", "nm", "sann")) {
   expect_silent(
      mlEq <- maxLik(llfInd, start = startVal, constraints = eqCon,
                     method = optimizer, SUMTTol = 0)
   )
   expect_stdout(
      print( summary( mlEq)),
      pattern = "constrained likelihood estimation. Inference is probably wrong.*outer iterations, barrier value"
   )
   expect_equal(sum(coef( mlEq )), 2.5, tolerance=1e-4)
}

### ---------- convergence tolerance parameters ----------
a <- maxNR(llf, gf, hf, start=startVal, tol=1e-3, reltol=0, gradtol=0, iterlim=10)
expect_equal(returnCode(a), 2)  # should stop with code 2: tolerance
a <- maxNR(llf, gf, hf, start=startVal, tol=0, reltol=1e-3, gradtol=0, iterlim=10)
expect_equal(returnCode(a), 8)  # 8: relative tolerance
a <- maxNR(llf, gf, hf, start=startVal, tol=0, reltol=0, gradtol=1e-3, iterlim=10)
expect_equal(returnCode(a), 1)  # 1: gradient
a <- maxNR(llf, gf, hf, start=startVal, tol=0, reltol=0, gradtol=0, iterlim=10)
expect_equal(returnCode(a), 4)  # 4: iteration limit
