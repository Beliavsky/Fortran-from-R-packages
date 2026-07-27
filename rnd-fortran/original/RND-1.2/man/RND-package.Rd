\name{RND-package}
\alias{RND-package}
\alias{RND}
\docType{package}
\title{
Risk Neutral Density Extraction Package
}
\description{
This package is a collection of various functions to extract the implied risk neutral density from option.
}
\details{
\tabular{ll}{
Package: \tab RND\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2017-01-10\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Kam Hamidieh <khamidieh@gmail.com>
}
\references{
  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London
}
\keyword{ package }
\keyword{ Risk Neutral Density }
\keyword{ Options }
\examples{

###
### You should see that all methods extract the same density!
###

r     = 0.05
te    = 60/365
s0    = 1000
sigma = 0.25
y     = 0.02

call.strikes.bsm   = seq(from = 500, to = 1500, by = 5)
market.calls.bsm   = price.bsm.option(r =r, te = te, s0 = s0, 
                     k = call.strikes.bsm, sigma = sigma, y = y)$call

put.strikes.bsm    = seq(from = 500, to = 1500, by = 5)
market.puts.bsm    = price.bsm.option(r =r, te = te, s0 = s0, 
                     k = put.strikes.bsm, sigma = sigma, y = y)$put

###
### See where your results will be outputted to...
###

getwd()


###
###  Running this may take a few minutes...
###
###  MOE(market.calls.bsm, call.strikes.bsm, market.puts.bsm, 
###  put.strikes.bsm, s0, r , te, y, "bsm2")
###
}
