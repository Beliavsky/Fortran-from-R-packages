\name{get.point.estimate}
\alias{get.point.estimate}
\title{Point Estimation of the Density}
\description{\code{get.point.estimate} estimates the risk neutral density by center differentiation.}
\usage{
get.point.estimate(market.calls, call.strikes, r, te)
}
\arguments{
  \item{market.calls}{market calls (most expensive to cheapest)}
  \item{call.strikes}{strikes for the calls (smallest to largest)}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
}
\details{This is a non-parametric estimate of the risk neutral density.  Due to center differentiation, 
the density values are not estimated at the highest and lowest strikes.}
\value{
  \item{point.estimates}{values of the estimated density at each strike}
}
\references{

  J. Hull (2011)
  \emph{Options, Futures, and Other Derivatives and DerivaGem Package}
  Prentice Hall, Englewood Cliffs, New Jersey, 8th Edition

}
\author{Kam Hamidieh}
\examples{

###
### Recover the lognormal density based on BSM
###

r     = 0.05
te    = 60/365
s0    = 1000
k     = seq(from = 500, to = 1500, by = 1)
sigma = 0.25
y     = 0.01


bsm.calls = price.bsm.option(r =r, te = te, s0 = s0, k = k, sigma = sigma, y = y)$call
density.est = get.point.estimate(market.calls = bsm.calls, 
              call.strikes = k, r = r , te = te)

len = length(k)-1
### Note, estimates at two data points (smallest and largest strikes) are lost
plot(density.est ~ k[2:len], type = "l")  

}

