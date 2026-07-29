\name{expgeo}
\alias{dexpgeo}
\alias{pexpgeo}
\alias{varexpgeo}
\alias{esexpgeo}
\title{Exponential geometric distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponential geometric distribution due to Adamidis and Loukas (1998) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {\lambda \theta \exp (-\lambda x)}{\left[ 1 - (1 - \theta) \exp (-\lambda x) \right]^2},
\\
&\displaystyle
F (x) = \frac {\theta \exp (-\lambda x)}{1 - (1 - \theta) \exp (-\lambda x)},
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{\lambda} \log \frac {p}{\theta + (1 - \theta) p},
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {\log p}{\lambda} - \frac {\theta \log \theta}{\lambda p (1 - \theta)} +
\frac {\theta + (1 - \theta) p}{\lambda p (1 - \theta)} \log \left[ \theta + (1 - \theta) p \right]
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{0 < \theta < 1}, the first scale parameter, and \eqn{\lambda > 0}, the second scale parameter.}
\usage{
dexpgeo(x, theta=0.5, lambda=1, log=FALSE)
pexpgeo(x, theta=0.5, lambda=1, log.p=FALSE, lower.tail=TRUE)
varexpgeo(p, theta=0.5, lambda=1, log.p=FALSE, lower.tail=TRUE)
esexpgeo(p, theta=0.5, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the first scale parameter, must be in the unit interval, the default is 0.5}
  \item{lambda}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexpgeo(x)
pexpgeo(x)
varexpgeo(x)
esexpgeo(x)}
 
