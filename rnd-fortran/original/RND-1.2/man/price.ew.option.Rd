\name{price.ew.option}
\alias{price.ew.option}
\title{Price Options with Edgeworth Approximated Density}
\description{
  \code{price.ew.option} computes the option prices based on Edgeworth approximated densities.
}
\usage{
price.ew.option(r, te, s0, k, sigma, y, skew, kurt)
}
\arguments{
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{s0}{current asset value}
  \item{k}{strike}
  \item{sigma}{volatility}
  \item{y}{dividend rate}
  \item{skew}{normalized skewness}
  \item{kurt}{normalized kurtosis}
}
\details{
Note that this function may produce negative prices if \code{skew} and \code{kurt} are
not well estimated from the data.
}
\value{
  \item{call }{Edgeworth based call}
  \item{put }{Edgeworth based put}
}
\references{
  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

  R. Jarrow and A. Rudd (1982)
  Approximate valuation for arbitrary stochastic processes.
  \emph{Journal of Finanical Economics}, 10, 347-369


  C.J. Corrado and T. Su (1996)
  S&P 500 index option tests of Jarrow and Rudd's approximate option valuation formula.
  \emph{Journal of Futures Markets}, 6, 611-629
}
\author{Kam Hamidieh}
\examples{

#
# Here, the prices must match EXACTLY the BSM prices:
#

r       = 0.05
y       = 0.03
s0      = 1000
sigma   = 0.25
te      = 100/365
k       = seq(from=800, to = 1200, by = 50)
v       = sqrt(exp(sigma^2 * te) - 1)
ln.skew = 3 * v + v^3
ln.kurt = 16 * v^2 + 15 * v^4 + 6 * v^6 + v^8

ew.option.prices  =  price.ew.option(r = r, te = te, s0 = s0, k=k, sigma=sigma, 
                                     y=y, skew = ln.skew, kurt = ln.kurt)
bsm.option.prices =  price.bsm.option(r = r, te = te, s0 = s0, k=k, sigma=sigma, y=y)

ew.option.prices
bsm.option.prices

###
### Now ew prices should be different as we increase the skewness and kurtosis:
###

new.skew = ln.skew * 1.10
new.kurt = ln.kurt * 1.10

new.ew.option.prices  =  price.ew.option(r = r, te = te, s0 = s0, k=k, sigma=sigma, 
                                         y=y, skew = new.skew, kurt = new.kurt)
new.ew.option.prices
bsm.option.prices

}

