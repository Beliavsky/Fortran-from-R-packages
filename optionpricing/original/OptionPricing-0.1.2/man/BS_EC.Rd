\name{BS_EC}
\alias{BS_EC}
\alias{BS_EP}
\title{
Black-Scholes Formula for European Call and Put
}
\description{
Calculates the Price, Delta and Gamma of an European Call or Put option using the Black-Scholes formula.
}
\usage{
BS_EC( T = 0.25, K = 100, r = 0.05, sigma = 0.2, S0 = 100 )
BS_EP( T = 0.25, K = 100, r = 0.05, sigma = 0.2, S0 = 100 )
}
\arguments{
  \item{T}{time to maturity (in years)}
  \item{K}{Strike Price}
  \item{r}{risk-free interest rate}
  \item{sigma}{yearly volatility}
  \item{S0}{Starting Stock Price}
}
\value{
Returns a vector containing the option price, Delta and Gamma
}
\author{Wolfgang Hormann}
\seealso{
\code{\link{OptionPricing-package}}
}
\examples{
BS_EC(K=100, r = 0.05, sigma = 0.2, T = 0.25, S0 = 100)
BS_EP(K=100, r = 0.05, sigma = 0.2, T = 0.25, S0 = 100)
}
