\name{extract.mln.density}
\alias{extract.mln.density}
\title{Extract Mixture of Lognormal Densities}
\description{
  \code{mln.extraction} extracts the parameters of the mixture of two lognormals densities.
}
\usage{
extract.mln.density(initial.values = c(NA, NA, NA, NA, NA), r, y, te, s0, 
  market.calls, call.strikes, call.weights = 1, market.puts, put.strikes, 
  put.weights = 1, lambda = 1, hessian.flag = F, cl = list(maxit = 10000))
}
\arguments{
  \item{initial.values}{initial values for the optimization}
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
  \item{hessian.flag}{if F, no hessian is produced}
  \item{cl}{list of parameter values to be passed to the optimization function}
}
\details{
  mln is the density f(x) = alpha.1 * g(x) + (1 - alpha.1) * h(x), where g and h are densities of two lognormals with
  parameters (mean.log.1, sdlog.1) and (mean.log.2, sdlog.2) respectively.
}
\value{
  \item{alpha.1}{extracted proportion of the first lognormal.  Second one is 1 - \code{alpha.1}}
  \item{meanlog.1}{extracted mean of the log of the first lognormal}
  \item{meanlog.2}{extracted mean of the log of the second lognormal}
  \item{sdlog.1}{extracted standard deviation of the log of the first lognormal}
  \item{sdlog.2}{extracted standard deviation of the log of the second lognormal}
  \item{converge.result }{Did the result converge?}
  \item{hessian}{ Hessian matrix}
}
\references{

  F. Gianluca and A. Roncoroni (2008)
  \emph{Implementing Models in Quantitative Finance: Methods and Cases}

  B. Bahra (1996):
  Probability distribution of future asset prices implied by option prices.
  \emph{Bank of England Quarterly Bulletin}, August 1996, 299-311

  P. Soderlind and L.E.O. Svensson (1997)
  New techniques to extract market expectations from financial instruments.
  \emph{Journal of Monetary Economics}, 40, 383-4

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{

#
# Create some calls and puts based on mln and 
# see if we can extract the correct values.
#


r         = 0.05
y         = 0.02
te        = 60/365
meanlog.1 = 6.8
meanlog.2 = 6.95
sdlog.1   = 0.065
sdlog.2   = 0.055
alpha.1   = 0.4


call.strikes = seq(from = 800, to = 1200, by = 10)
market.calls = price.mln.option(r = r, y = y, te = te, k = call.strikes, 
               alpha.1 = alpha.1, meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, 
                                sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)$call

s0 = price.mln.option(r = r, y = y, te = te, k = call.strikes, alpha.1 = alpha.1, 
                      meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, 
                      sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)$s0
s0
put.strikes  = seq(from = 805, to = 1200, by = 10)
market.puts  = price.mln.option(r = r, y = y, te = te, k = put.strikes, 
                                alpha.1 = alpha.1, meanlog.1 = meanlog.1, 
                                meanlog.2 = meanlog.2, sdlog.1 = sdlog.1, 
                                sdlog.2 = sdlog.2)$put

###
### The extracted values should be close to the actual values.
###

extract.mln.density(r = r, y = y, te = te, s0 = s0, market.calls = market.calls, 
               call.strikes = call.strikes, market.puts = market.puts, 
               put.strikes = put.strikes, lambda = 1, hessian.flag = FALSE)


}

