\name{exppois}
\alias{dexppois}
\alias{pexppois}
\alias{varexppois}
\alias{esexppois}
\title{Exponential Poisson distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponential Poisson distribution due to Kus (2007) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {b \lambda \exp \left[ -b x - \lambda + \lambda \exp (-b x) \right]}{1 - \exp (-\lambda)},
\\
&\displaystyle
F (x) = \frac {1 - \exp \left[ -\lambda + \lambda \exp (-b x) \right]}{1 - \exp (-\lambda)},
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{b} \log \left\{ \frac {1}{\lambda} \log \left[ 1 - p + p \exp (-\lambda) \right] + 1 \right\},
\\
&\displaystyle
{\rm ES}_p (X) =  -\frac {1}{b p} \int_0^p \log \left\{ \frac {1}{\lambda} \log
\left[ 1 - v + v \exp (-\lambda) \right] + 1 \right\} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{b > 0}, the first scale parameter, and \eqn{\lambda > 0}, the second scale parameter.}
\usage{
dexppois(x, b=1, lambda=1, log=FALSE)
pexppois(x, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
varexppois(p, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
esexppois(p, b=1, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{b}{the value of the first scale parameter, must be positive, the default is 1}
  \item{lambda}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexppois(x)
pexppois(x)
varexppois(x)
esexppois(x)}
 

