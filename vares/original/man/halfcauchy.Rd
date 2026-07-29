\name{halfcauchy}
\alias{dhalfcauchy}
\alias{phalfcauchy}
\alias{varhalfcauchy}
\alias{eshalfcauchy}
\title{Half Cauchy distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the half Cauchy distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {2}{\pi} \frac {\sigma}{x^2 + \sigma^2},
\\
&\displaystyle
F (x) = \frac {2}{\pi} \arctan \left( \frac {x}{\sigma} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \sigma \tan \left( \frac {\pi p}{2} \right),
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sigma}{p} \int_0^p \tan \left( \frac {\pi v}{2} \right) dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, and \eqn{\sigma > 0}, the scale parameter.}
\usage{
dhalfcauchy(x, sigma=1, log=FALSE)
phalfcauchy(x, sigma=1, log.p=FALSE, lower.tail=TRUE)
varhalfcauchy(p, sigma=1, log.p=FALSE, lower.tail=TRUE)
eshalfcauchy(p, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dhalfcauchy(x)
phalfcauchy(x)
varhalfcauchy(x)
eshalfcauchy(x)}
 

