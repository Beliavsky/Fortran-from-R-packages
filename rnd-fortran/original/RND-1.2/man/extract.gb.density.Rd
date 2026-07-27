\name{extract.gb.density}
\alias{extract.gb.density}
\title{Generalized Beta Extraction}
\description{
  \code{extract.gb.density} extracts the generalized beta density from market options.
}
\usage{
extract.gb.density(initial.values = c(NA, NA, NA, NA), r, te, y, s0, market.calls, 
  call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, 
  lambda = 1, hessian.flag = F, cl = list(maxit = 10000))
}
\arguments{
  \item{initial.values}{initial values for the optimization}
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
  \item{hessian.flag}{if F, no hessian is produced}
  \item{cl}{list of parameter values to be passed to the optimization function}
}
\details{
  This function extracts the generalized beta density implied by the options.
}
\value{
  \item{a}{extracted power parameter}
  \item{b}{extracted scale paramter}
  \item{v}{extracted first beta paramter}
  \item{w}{extracted second beta parameter}
  \item{converge.result }{Did the result converge?}
  \item{hessian}{ Hessian matrix}
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
# create some GB based calls and puts
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
market.calls = price.gb.option(r = r, te = te, y = y, s0 = s0, 
                               k = call.strikes, a = a, b = s0, v = v, w = w)$call

put.strikes  = seq(from = 805, to = 1200, by = 10)
market.puts  = price.gb.option(r = r, te = te, y = y, s0 = s0, 
                               k = put.strikes, a = a, b = s0, v = v, w = w)$put


#
# Extraction...should match the a,b,v,w above. You will also get warning messages.
# Weigths are automatically set to 1.
#

extract.gb.density(r=r, te=te, y = y, s0=s0, market.calls = market.calls, 
              call.strikes = call.strikes, market.puts = market.puts, 
              put.strikes = put.strikes, lambda = 1, hessian.flag = FALSE)

}