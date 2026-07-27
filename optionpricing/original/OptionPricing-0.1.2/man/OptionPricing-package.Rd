\name{OptionPricing-package}
\alias{OptionPricing-package}
\alias{OptionPricing}
\docType{package}
\title{
Option Pricing and Greeks Estimation for Asian and European Options
}
\description{
The Price, Delta and Gamma of European and Asian Options under Geometric Brownian Motion
are calculated using the Black-Scholes formula and
Efficient Monte Carlo and Randomized Quasi Monte Carlo Algorithms.}
\details{
The OptionPricing package calculates the Price, Delta and Gamma for European options
using the Black-Scholes formula (see \code{\link{BS_EC}}).
The price, Delta and Gamma for Asian call options under geometric Brownian motion are calculated using
a very efficient Monte Carlo and randomized quasi-Monte Carlo algorithm (see \code{\link{AsianCall}}). 
The function \code{\link{AsianCall_AppLord}} implements
a high-quality approximation for the price of an Asian option.} 
\author{
Kemal Dingec, Wolfgang Hormann
}
\keyword{ package }
\examples{
# standard settings for an efficient simulation using QMC and variance reduction
AsianCall(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method=c("best"),
  sampling=c("QMC"),metpar=list(maxiter=100,tol=1.e-14,cvmethod="splitting"),
  sampar=list(nout=50,n=2039,a=1487,baker=TRUE,genmethod="pca"))

# Calculation of the Price of an Asian option using a good approximation  
  AsianCall_AppLord(T = 1, d = 12, K = 100, r = 0.05, sigma = 0.2, S0 = 100)


# standard settings for an efficient simulation using MC and variance reduction
AsianCall(T=1,d=12,K=170,r=0.05,sigma=0.2,S0=100,method="best",
          sampling="MC",metpar=list(maxiter=100,tol=1.e-14,np=1000),
		  sampar=list(n=10^5))
# Calculation of the approximate price, a bit different to the above result 
  AsianCall_AppLord(T = 1, d = 12, K = 170, r = 0.05, sigma = 0.2, S0 = 100)
  
# Calculation of the Price of an Asian option using a good approximation  
  AsianCall_AppLord(T = 1, d = 12, K = 100, r = 0.05, sigma = 0.2, S0 = 100)

#Price, Delta and Gamma of European options using Black-Scholes  
BS_EC(K=100, r = 0.05, sigma = 0.2, T = 0.25, S0 = 100)
BS_EP(K=100, r = 0.05, sigma = 0.2, T = 0.25, S0 = 100)
}
