\name{chen}
\alias{dchen}
\alias{pchen}
\alias{varchen}
\alias{eschen}
\title{Chen distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Chen distribution due to Chen (2000) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \lambda b x^{b - 1} \exp \left( x^b \right) \exp \left[ \lambda - \lambda \exp \left( x^b \right) \right],
\\
&\displaystyle
F (x) = 1 - \exp \left[ \lambda - \lambda \exp \left( x^b \right) \right],
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{ \log \left[ 1 - \frac {\log (1 - p)}{\lambda} \right] \right\}^{1 / b},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p \left\{ \log \left[ 1 - \frac {\log (1 - v)}{\lambda} \right] \right\}^{1 / b} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{b > 0}, the shape parameter, and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dchen(x, b=1, lambda=1, log=FALSE)
pchen(x, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
varchen(p, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
eschen(p, b=1, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dchen(x)
pchen(x)
varchen(x)
eschen(x)}
 
