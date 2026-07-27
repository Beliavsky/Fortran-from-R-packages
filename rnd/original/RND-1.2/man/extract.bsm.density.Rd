\name{extract.bsm.density}
\alias{extract.bsm.density}
\title{Extract Lognormal Density}
\description{
  \code{bsm.extraction} extracts the parameters of the lognormal density as implied by the BSM model.
}
\usage{
extract.bsm.density(initial.values = c(NA, NA), r, y, te, s0, market.calls, 
  call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, 
  lambda = 1, hessian.flag = F, cl = list(maxit = 10000))
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
If \code{initial.values} are not specified then the function will attempt to pick them
automatically.  \code{cl} is a list that can be used to pass parameters to the \code{optim}
function.
}
\value{
  Let S_T with the lognormal random variable of the risk neutral density.

  \item{mu}{mean of log(S_T)}
  \item{zeta }{sd of log(S_T)}
  \item{converge.result }{Did the result converge?}
  \item{hessian}{ Hessian matrix}
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
# Create some BSM Based options
#

r     = 0.05
te    = 60/365
s0    = 1000
sigma = 0.25
y     = 0.01

call.strikes   = seq(from = 500, to = 1500, by = 25)
market.calls   = price.bsm.option(r =r, te = te, s0 = s0, 
                 k = call.strikes, sigma = sigma, y = y)$call

put.strikes    = seq(from = 510, to = 1500, by = 25)
market.puts    = price.bsm.option(r =r, te = te, s0 = s0, 
                 k = put.strikes, sigma = sigma, y = y)$put

#
#  Get extract the parameter of the density
#

extract.bsm.density(r = r, y = y, te = te, s0 = s0, market.calls = market.calls, 
               call.strikes = call.strikes,  market.puts = market.puts, 
               put.strikes = put.strikes, lambda = 1, hessian.flag = FALSE)

#
# The extracted parameters should be close to these actual values:
#
actual.mu     = log(s0) + ( r - y - 0.5 * sigma^2) * te
actual.zeta   = sigma * sqrt(te)
actual.mu 
actual.zeta
}

