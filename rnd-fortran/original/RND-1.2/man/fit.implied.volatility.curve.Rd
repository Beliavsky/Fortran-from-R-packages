\name{fit.implied.volatility.curve}
\alias{fit.implied.volatility.curve}
\title{Fit Implied Quadratic Volatility Curve}
\description{
  \code{fit.implied.volatility.curve} estimates the coefficients of the quadratic equation for the implied volatilities.
}
\usage{
fit.implied.volatility.curve(x, k)
}
\arguments{
  \item{x}{a set of implied volatilities}
  \item{k}{range of strikes}
}
\details{This function estimates volatility \eqn{\sigma} as a quadratic function of strike \eqn{k} with 
 the coefficents \eqn{a_0, a_1, a_2}: \eqn{\sigma(k) = a_0 + a_1 k + a_2 k^2}  }
\value{
  \item{a0}{constant term in the quadratic ploynomial}
  \item{a1}{coefficient term of k in the quadratic ploynomial}
  \item{a2}{coefficient term of k squared in the quadratic polynomial}
  \item{summary.obj }{statistical summary of the fit}
}
\references{
  D. Shimko (1993)
  Bounds of probability.
  \emph{Risk}, 6, 33-47

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{

#
# Suppose we see the following implied volatilities and strikes:
#

implied.sigma = c(0.11, 0.08, 0.065, 0.06, 0.05)  
strikes       = c(340, 360, 380, 400, 410)
tmp           = fit.implied.volatility.curve(x = implied.sigma, k = strikes)
tmp

strike.range = 340:410
plot(implied.sigma ~ strikes)
lines(strike.range, tmp$a0 + tmp$a1 * strike.range + tmp$a2 * strike.range^2)

}

