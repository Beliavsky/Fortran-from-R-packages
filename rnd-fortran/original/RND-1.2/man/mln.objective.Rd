\name{mln.objective}
\alias{mln.objective}
\title{Objective function for the Mixture of Lognormal}
\description{
   \code{mln.objective} is the objective function to be minimized in \code{extract.mln.density}.
}
\usage{
mln.objective(theta, r, y, te, s0, market.calls, call.strikes, call.weights, 
  market.puts, put.strikes, put.weights, lambda = 1)
}
\arguments{
  \item{theta}{initial values for the optimization}
  \item{r}{risk free rate}
  \item{y}{dividend yield}
  \item{te}{time to expiration}
  \item{s0}{current asset value}
  \item{market.calls}{market calls (most expensive to cheapest)}
  \item{call.strikes}{strikes for the calls (smallest to largest)}
  \item{call.weights}{weights to be used for calls}
  \item{market.puts}{market calls (cheapest to most expensive)}
  \item{put.strikes}{strikes for the puts (smallest to largest)}
  \item{put.weights}{weights to be used for puts}
  \item{lambda}{Penalty parameter to enforce the martingale condition}
}
\details{
  mln is the density f(x) = alpha.1 * g(x) + (1 - alpha.1) * h(x), where g and h are densities of two lognormals with
  parameters (mean.log.1, sdlog.1) and (mean.log.2, sdlog.2) respectively.
}
\value{
  \item{obj }{value of the objective function}
}
\references{
  F. Gianluca and A. Roncoroni (2008)
  \emph{Implementing Models in Quantitative Finance: Methods and Cases}

  B. Bahra (1996):
  Probability distribution of future asset prices implied by option prices.
  \emph{Bank of England Quarterly Bulletin}, August 1996, 299-311

  P. Soderlind and L.E.O. Svensson (1997)
  New techniques to extract market expectations from financial instruments.
  \emph{Journal of Monetary Economics}, 40, 383-429

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{
#
# The mln objective function should be close to zero.
# The weights are automatically set to 1.
#

r  = 0.05
te = 60/365
y  = 0.02
   
meanlog.1 = 6.8
meanlog.2 = 6.95
sdlog.1   = 0.065
sdlog.2   = 0.055
alpha.1   = 0.4

# This is the current price implied by parameter values:
s0 = 981.8815 

call.strikes = seq(from = 800, to = 1200, by = 10)
market.calls = price.mln.option(r=r, y = y, te = te, k = call.strikes, 
               alpha.1 = alpha.1, meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, 
               sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)$call

put.strikes  = seq(from = 805, to = 1200, by = 10)
market.puts  = price.mln.option(r = r, y = y, te = te, k = put.strikes, 
               alpha.1 = alpha.1, meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, 
               sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)$put

mln.objective(theta=c(alpha.1,meanlog.1, meanlog.2 , sdlog.1, sdlog.2), 
               r = r, y = y, te = te, s0 = s0, 
               market.calls = market.calls, call.strikes = call.strikes, 
               market.puts = market.puts, put.strikes = put.strikes, lambda = 1)
}

