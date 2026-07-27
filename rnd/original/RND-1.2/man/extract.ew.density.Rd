\name{extract.ew.density}
\alias{extract.ew.density}
\title{Extract Edgeworth Based Density}
\description{
  \code{ew.extraction} extracts the parameters for the density approximated by the Edgeworth expansion method.
}
\usage{
extract.ew.density(initial.values = c(NA, NA, NA), r, y, te, s0, market.calls, 
  call.strikes, call.weights = 1, lambda = 1, hessian.flag = F, 
  cl = list(maxit = 10000))
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
  \item{lambda}{Penalty parameter to enforce the martingale condition}
  \item{hessian.flag}{if F, no hessian is produced}
  \item{cl}{list of parameter values to be passed to the optimization function}
}
\details{If \code{initial.values} are not specified then the function will attempt to pick them
automatically.  \code{cl} in form of a list can be used to pass parameters to the \code{optim}
function.}
\value{
  \item{sigma }{volatility of the underlying lognormal}
  \item{skew }{normalized skewness}
  \item{kurt }{normalized kurtosis}
  \item{converge.result }{Did the result converge?}
  \item{hessian}{ Hessian matrix}
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
#  ln.skew & ln.kurt are the normalized skewness and kurtosis of a true lognormal.
#

r       = 0.05
y       = 0.03
s0      = 1000
sigma   = 0.25
te      = 100/365
strikes = seq(from=600, to = 1400, by = 1)
v       = sqrt(exp(sigma^2 * te) - 1)
ln.skew = 3 * v + v^3
ln.kurt = 16 * v^2 + 15 * v^4 + 6 * v^6 + v^8

#
# Now "perturb" the lognormal
#

new.skew = ln.skew * 1.50
new.kurt = ln.kurt * 1.50

#
# new.skew & new.kurt should not be extracted.
# Note that weights are automatically set to 1.
#

market.calls      =  price.ew.option(r = r, te = te, s0 = s0, k=strikes, sigma=sigma, 
                     y=y, skew = new.skew, kurt = new.kurt)$call
ew.extracted.obj  =  extract.ew.density(r = r, y = y, te = te, s0 = s0, 
                     market.calls = market.calls, call.strikes = strikes, 
                     lambda = 1, hessian.flag = FALSE)
ew.extracted.obj
}

