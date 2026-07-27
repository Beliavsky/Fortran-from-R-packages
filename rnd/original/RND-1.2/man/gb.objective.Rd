\name{gb.objective}
\alias{gb.objective}
\title{Generalized Beta Objective}
\description{
  \code{gb.objective} is the objective function to be minimized in \code{extract.gb.density}.
}
\usage{
gb.objective(theta, r, te, y, s0, market.calls, call.strikes, call.weights = 1, 
  market.puts, put.strikes, put.weights = 1, lambda = 1)
}
\arguments{
  \item{theta}{initial values for optimization}
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{y}{dividend yield}
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
  This is the function minimized by \code{extract.gb.desnity} function.
}
\value{
  \item{obj }{value of the objective function}
}
\references{
  R.M. Bookstaber and J.B. McDonald (1987)
  A general distribution for describing security price returns.
  \emph{Journal of Business}, 60, 401-424 

  X. Liu and M.B. Shackleton and S.J. Taylor and X. Xu (2007)
  Closed-form transformations from risk-neutral to real-world distributions
  \emph{Journal of Business}, 60, 401-424 

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{

#
# The objective should be very small!
# Note the weights are automatically
# set to 1.
#

r  = 0.03
te = 50/365
k  = seq(from = 800, to = 1200, by = 10)
a  = 10
b  = 1000
v  = 2.85
w  = 2.85
y  = 0.01
s0 = exp((y-r)*te) * b * beta(v + 1/a, w - 1/a)/beta(v,w) 
s0

call.strikes = seq(from = 800, to = 1200, by = 10)
market.calls = price.gb.option(r = r, te = te, s0 = s0, y = y, 
                        k = call.strikes, a = a, b = b, v = v, w = w)$call

put.strikes = seq(from = 805, to = 1200, by = 10)
market.puts = price.gb.option(r = r, te = te, s0 = s0, y = y, 
                        k = put.strikes, a = a, b = b, v = v, w = w)$put

gb.objective(theta=c(a,b,v,w),r = r, te = te, y = y, s0 = s0, 
             market.calls = market.calls, call.strikes = call.strikes,  
             market.puts = market.puts, put.strikes = put.strikes, lambda = 1)


}

