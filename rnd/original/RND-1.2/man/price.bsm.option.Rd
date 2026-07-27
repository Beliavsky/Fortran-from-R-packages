\name{price.bsm.option}
\alias{price.bsm.option}
\title{Price BSM Option}
\description{
  \code{bsm.option.price} computes the BSM European option prices.
}
\usage{
price.bsm.option(s0, k, r, te, sigma, y)
}
\arguments{
  \item{s0}{current asset value}
  \item{k}{strike}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{sigma}{volatility}
  \item{y}{dividend yield}
}
\details{
  This function implements the classic Black-Scholes-Merton option pricing model.
}
\value{
  \item{d1 }{value of \code{(log(s0/k) + (r - y + (sigma^2)/2) * te)/(sigma * sqrt(te))}}
  \item{d2 }{value of \code{d1 - sigma * sqrt(te)}}
  \item{call }{call price}
  \item{put }{put price}
}
\references{

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

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
# call should be 4.76, put should be 0.81, from Hull 8th, page 315, 316
#

r     = 0.10
te    = 0.50
s0    = 42
k     = 40
sigma = 0.20
y     = 0

bsm.option = price.bsm.option(r =r, te = te, s0 = s0, k = k, sigma = sigma, y = y)
bsm.option

#
# Make sure put-call parity holds, Hull 8th, page 351
#

(bsm.option$call - bsm.option$put) - (s0 * exp(-y*te) - k * exp(-r*te))

}

