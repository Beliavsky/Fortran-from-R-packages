\name{compute.implied.volatility}
\alias{compute.implied.volatility}
\title{Compute Impied Volatility}
\description{
  \code{compute.implied.volatility} extracts the implied volatility for a call option.
}
\usage{
compute.implied.volatility(r, te, s0, k, y, call.price, lower, upper)
}
\arguments{
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{s0}{current asset value}
  \item{k}{strike of the call option}
  \item{y}{dividend yield}
  \item{call.price}{call price}
  \item{lower}{lower bound of the implied volatility to look for}
  \item{upper}{upper bound of the implied volatility to look for}
}
\details{The simple R \code{uniroot} function is used to extract the implied volatility.}
\value{
  \item{sigma }{extratced implied volatility}
}
\references{
  J. Hull (2011)
  \emph{Options, Futures, and Other Derivatives and DerivaGem Package}
  Prentice Hall, Englewood Cliffs, New Jersey, 8th Edition

  R. L. McDonald (2013)
  \emph{Derivatives Markets}
  Pearson, Upper Saddle River, New Jersey, 3rd Edition

}
\author{Kam Hamidieh}
\examples{
#
# Create prices from BSM with various sigma's
#

r     =  0.05
y     =  0.02
te    =  60/365
s0    =  400

sigma.range = seq(from = 0.1, to = 0.8, by = 0.05)
k.range     = floor(seq(from = 300, to = 500, length.out = length(sigma.range)))
bsm.calls   = numeric(length(sigma.range))

for (i in 1:length(sigma.range))
{
  bsm.calls[i] = price.bsm.option(r = r, te = te, s0 = s0, k = k.range[i], 
                                  sigma = sigma.range[i], y = y)$call
}
bsm.calls
k.range

#
# Computed implied sigma's should be very close to sigma.range.
#

compute.implied.volatility(r = r, te = te, s0 = s0, k = k.range, y = y, 
                          call.price = bsm.calls, lower = 0.001, upper = 0.999)
sigma.range

}

