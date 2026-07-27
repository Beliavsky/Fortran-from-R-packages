\name{extract.rates}
\alias{extract.rates}
\title{Extract Risk Free Rate and Dividend Yield}
\description{
  \code{extract.rates} extracts the risk free rate and the dividend yield from European options.
}
\usage{
extract.rates(calls, puts, s0, k, te)
}
\arguments{
  \item{calls}{market calls (most expensive to cheapest)}
  \item{puts}{market puts (cheapest to most expensive)}
  \item{s0}{current asset value}
  \item{k}{strikes for the calls (smallest to largest)}
  \item{te}{time to expiration}
}
\details{
  The extraction is based on the put-call parity of the European options.  Shimko (1993) 
  - see below - shows that the slope and intercept of the regression of 
  the calls minus puts onto the strikes contains the risk free and the dividend
  rates.
}
\value{
  \item{risk.free.rate }{extracted risk free rate}
  \item{dividend.yield }{extracted dividend rate}
}
\references{
  D. Shimko (1993)
  Bounds of probability.
  \emph{Risk}, 6, 33-47
}
\author{Kam Hamidieh}
\examples{
#
# Create calls and puts based on BSM
#

r     = 0.05
te    = 60/365
s0    = 1000
k     = seq(from = 900, to = 1100, by = 25)
sigma = 0.25
y     = 0.01

bsm.obj = price.bsm.option(r =r, te = te, s0 = s0, k = k, sigma = sigma, y = y)

calls = bsm.obj$call
puts  = bsm.obj$put

#
# Extract rates should give the values of r and y above:
#

rates = extract.rates(calls = calls, puts = puts, k = k, s0 = s0, te = te)
rates

}

