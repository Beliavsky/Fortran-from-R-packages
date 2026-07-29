\name{expexp}
\alias{dexpexp}
\alias{pexpexp}
\alias{varexpexp}
\alias{esexpexp}
\title{Exponentiated exponential distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponentiated exponential distribution due to Gupta and Kundu (1999, 2001) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = a \lambda \exp (-\lambda x)
\left[ 1 - \exp (-\lambda x) \right]^{a - 1},
\\
&\displaystyle
F (x) = \left[ 1 - \exp (-\lambda x) \right]^{a},
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{\lambda} \log \left( 1 - p^{1 / a} \right),
\\
&\displaystyle
{\rm ES}_p (X) =  -\frac {1}{p \lambda} \int_0^p \log \left( 1 - v^{1 / a} \right) dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dexpexp(x, lambda=1, a=1, log=FALSE)
pexpexp(x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
varexpexp(p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
esexpexp(p, lambda=1, a=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexpexp(x)
pexpexp(x)
varexpexp(x)
esexpexp(x)}
 
