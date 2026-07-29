\name{exponential}
\alias{dexponential}
\alias{pexponential}
\alias{varexponential}
\alias{esexponential}
\title{Exponential distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponential distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \lambda \exp (-\lambda x),
\\
&\displaystyle
F (x) = 1 - \exp (-\lambda x),
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{\lambda} \log (1 - p),
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {1}{p \lambda} \left\{ \log (1 - p) p - p - \log (1 - p) \right\}
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dexponential(x, lambda=1, log=FALSE)
pexponential(x, lambda=1, log.p=FALSE, lower.tail=TRUE)
varexponential(p, lambda=1, log.p=FALSE, lower.tail=TRUE)
esexponential(p, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexponential(x)
pexponential(x)
varexponential(x)
esexponential(x)}
 
