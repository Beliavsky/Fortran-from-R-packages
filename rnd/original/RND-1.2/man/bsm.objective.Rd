\name{bsm.objective}
\alias{bsm.objective}
\title{BSM Objective Function}
\description{
  \code{bsm.objective} is the objective function to be minimized in \code{extract.bsm.density}.
}
\usage{
bsm.objective(s0, r, te, y, market.calls, call.strikes, call.weights = 1, 
  market.puts, put.strikes, put.weights = 1, lambda = 1, theta)
}

\arguments{
  \item{s0}{current asset value}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{y}{dividend yield}
  \item{market.calls}{market calls (most expensive to cheapest)}
  \item{call.strikes}{strikes for the calls (smallest to largest)}
  \item{call.weights}{weights to be used for calls}
  \item{market.puts}{market calls (cheapest to most expensive)}
  \item{put.strikes}{strikes for the puts (smallest to largest)}
  \item{put.weights}{weights to be used for calls}
  \item{lambda}{Penalty parameter to enforce the martingale condition}
  \item{theta}{initial values for the optimization.  This must be a vector
        of length 2: first component is \eqn{\mu}, the lognormal mean of
        the underlying density, and the second component is \eqn{\sqrt{t}\sigma} which is the
        time scaled volatility parameter of the underlying density.}
}
\details{
This function evaluates the weighted squared differences between the market option values and
values predicted by the Black-Scholes-Merton option pricing formula.
}
\value{
Objective function evalued at a specific set of values.
}
\references{
  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London
}
\author{Kam Hamidieh}
\examples{
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

###
### perfect initial values under BSM framework
###

mu.0     = log(s0) + ( r - y - 0.5 * sigma^2) * te
zeta.0   = sigma * sqrt(te)
mu.0
zeta.0

###
### The objective function should be *very* small
###

bsm.obj.val = bsm.objective(theta=c(mu.0, zeta.0), r = r, y=y, te = te, s0 = s0, 
              market.calls = market.calls, call.strikes = call.strikes, 
              market.puts  = market.puts,  put.strikes = put.strikes, lambda = 1)
bsm.obj.val
}
