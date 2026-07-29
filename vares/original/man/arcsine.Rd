\name{arcsine}
\alias{darcsine}
\alias{parcsine}
\alias{vararcsine}
\alias{esarcsine}
\title{Arcsine distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the arcsine distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\pi \sqrt{(x - a) (b - x)}},
\\
&\displaystyle
F (x) = \frac {2}{\pi} \arcsin \left( \sqrt{\frac {x - a}{b - a}} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = a + (b - a) \sin^2 \left( \frac {\pi p}{2} \right),
\\
&\displaystyle
{\rm ES}_p (X) = a + \frac {b - a}{p} \int_0^p \sin^2 \left( \frac {\pi v}{2} \right) dv
\end{array}}
for \eqn{a \leq x \leq b}, \eqn{0 < p < 1}, \eqn{-\infty < a < \infty}, the  first location parameter, and \eqn{-\infty < a < b < \infty}, the second location parameter.}
\usage{
darcsine(x, a=0, b=1, log=FALSE)
parcsine(x, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
vararcsine(p, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
esarcsine(p, a=0, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first location parameter, can take any real value, the default is zero}
  \item{b}{the value of the second location parameter, can take any real value but must be greater than a, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
darcsine(x)
parcsine(x)
vararcsine(x)
esarcsine(x)}
 
