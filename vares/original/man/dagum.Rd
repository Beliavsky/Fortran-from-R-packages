\name{dagum}
\alias{ddagum}
\alias{pdagum}
\alias{vardagum}
\alias{esdagum}
\title{Dagum distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Dagum distribution due to Dagum (1975, 1977, 1980) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {a c b^a x^{a c - 1}}{\left[ x^a + b^a \right]^{c + 1}},
\\
&\displaystyle
F (x) = \left[ 1 + \left( \frac {b}{x} \right)^a \right]^{-c},
\\
&\displaystyle
{\rm VaR}_p (X) = b \left( 1- p^{-1 / c} \right)^{-1 / a},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {b}{p} \int_0^p \left( 1 - v^{-1 / c} \right)^{-1 / a} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, \eqn{b > 0}, the scale parameter, and \eqn{c > 0}, the second shape parameter.}
\usage{
ddagum(x, a=1, b=1, c=1, log=FALSE)
pdagum(x, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
vardagum(p, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
esdagum(p, a=1, b=1, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{b}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{c}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
ddagum(x)
pdagum(x)
vardagum(x)
esdagum(x)}
 


