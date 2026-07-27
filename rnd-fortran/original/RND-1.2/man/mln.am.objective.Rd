\name{mln.am.objective}
\alias{mln.am.objective}
\title{Objective function for the Mixture of Lognormal of American Options}
\description{
   \code{mln.am.objective} is the objective function to be minimized in \code{extract.am.density}.
}
\usage{
mln.am.objective(theta, s0, r, te, market.calls, call.weights = NA, market.puts, 
  put.weights = NA, strikes, lambda = 1)
}
\arguments{
  \item{theta}{initial values for the optimization}
  \item{s0}{current asset value}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{market.calls}{market calls (most expensive to cheapest)}
  \item{call.weights}{weights to be used for calls}
  \item{market.puts}{market calls (cheapest to most expensive)}
  \item{put.weights}{weights to be used for calls}
  \item{strikes}{strikes for the calls (smallest to largest)}
  \item{lambda}{Penalty parameter to enforce the martingale condition}
}
\details{
  mln is density f(x) = p.1 * f1(x) + p.2 * f2(x) + (1 - p.1 - p.2) * f3(x), where
  f1, f2, and f3 are lognormal densities with log means u.1,u.2, and u.3 and
  standard deviations sigma.1, sigma.2, and sigma.3 respectively.
}
\value{
  \item{obj }{Value of the objective function}
}
\references{

  Melick, W. R. and Thomas, C. P. (1997). Recovering an asset's implied pdf from option prices: An
  application to crude oil during the gulf crisis. \emph{Journal of Financial and Quantitative Analysis},
  32(1), 91-115.

}
\author{Kam Hamidieh}
\examples{

r       = 0.01
te      = 60/365
w.1     = 0.4
w.2     = 0.25
u.1     = 4.2
u.2     = 4.5
u.3     = 4.8
sigma.1 = 0.30
sigma.2 = 0.20
sigma.3 = 0.15
p.1     = 0.25
p.2     = 0.45
theta   = c(w.1,w.2,u.1,u.2,u.3,sigma.1,sigma.2,sigma.3,p.1,p.2)

p.3 = 1 - p.1 - p.2
p.3
expected.f0   =  sum(c(p.1, p.2, p.3) * exp(c(u.1,u.2,u.3) + 
                    (c(sigma.1, sigma.2, sigma.3)^2)/2) )
expected.f0  
 
strikes = 30:170

market.calls = numeric(length(strikes))
market.puts  = numeric(length(strikes))

for (i in 1:length(strikes))
{

  if ( strikes[i] < expected.f0) {
    market.calls[i] = price.am.option(k = strikes[i], r = r, te = te, w = w.1, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$call.value

    market.puts[i]  = price.am.option(k = strikes[i], r = r, te = te, w = w.2, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$put.value 
  }  else {

    market.calls[i] = price.am.option(k = strikes[i], r = r, te = te, w = w.2, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$call.value

    market.puts[i]  = price.am.option(k = strikes[i], r = r, te = te, w = w.1, u.1 = u.1, 
                      u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, 
                      sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2)$put.value 
     }

}

###
### Quickly look at the option values...
###

par(mfrow=c(1,2))
plot(market.calls ~ strikes, type="l")
plot(market.puts  ~ strikes, type="l")
par(mfrow=c(1,1))

###
### ** IMPORTANT **:  The code that follows may take a few seconds.
###                   Copy and paste onto R console the commands
###                   that follow the greater sign >.
###
###
### Next try the objective function.  It should be zero.
### Note: Let weights be the defaults values of 1.
###
#
# > s0      = expected.f0 * exp(-r * te)
# > s0
#
# > mln.am.objective(theta, s0 =s0, r = r, te = te, market.calls = market.calls,  
#                 market.puts = market.puts, strikes = strikes, lambda = 1)
#
###
### Now directly try the optimization with perfect initial values.
###
#
#
# > optim.obj.with.synthetic.data = optim(theta, mln.am.objective, s0 = s0, r=r, te=te, 
#                 market.calls = market.calls, market.puts = market.puts, strikes = strikes, 
#                 lambda = 1, hessian = FALSE , control=list(maxit=10000) )
#
# > optim.obj.with.synthetic.data
#
# > theta
#
###
### It does take a few seconds but the optim converges to exact theta values.
###
}

