\name{loglog}
\alias{dloglog}
\alias{ploglog}
\alias{varloglog}
\alias{esloglog}
\title{Loglog distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Loglog distribution due to Pham (2002) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = a \log (\lambda) x^{a - 1}
\lambda^{x^a} \exp \left[ 1 - \lambda^{x^a} \right],
\\
&\displaystyle
F (x) = 1 - \exp \left[ 1 - \lambda^{x^a} \right],
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{ \frac {\log \left[ 1 - \log (1 - p) \right]}{\log \lambda} \right\}^{1/a},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p (\log \lambda)^{1/a}}
\int_0^p \left\{ \log \left[ 1 - \log (1 - v) \right] \right\}^{1/a} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter, and \eqn{\lambda > 1}, the scale parameter.}
\usage{
dloglog(x, a=1, lambda=2, log=FALSE)
ploglog(x, a=1, lambda=2, log.p=FALSE, lower.tail=TRUE)
varloglog(p, a=1, lambda=2, log.p=FALSE, lower.tail=TRUE)
esloglog(p, a=1, lambda=2)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be greater than 1, the default is 2}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dloglog(x)
ploglog(x)
varloglog(x)
esloglog(x)}
 
