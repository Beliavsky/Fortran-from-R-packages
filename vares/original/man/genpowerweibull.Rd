\name{genpowerweibull}
\alias{dgenpowerweibull}
\alias{pgenpowerweibull}
\alias{vargenpowerweibull}
\alias{esgenpowerweibull}
\title{Generalized power Weibull distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the generalized power Weibull distribution due to Nikulin and Haghighi (2006) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = a \theta x^{a - 1} \left[ 1 + x^a \right]^{\theta - 1}
\exp \left\{ 1 - \left[ 1 + x^a \right]^\theta \right\},
\\
&\displaystyle
F (x) = 1 - \exp \left\{ 1 - \left[ 1 + x^a \right]^\theta \right\},
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{ \left[ 1 - \log (1 - p) \right]^{1 / \theta} - 1 \right\}^{1 / a},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p \left\{ \left[ 1 - \log (1 - v) \right]^{1 / \theta} - 1 \right\}^{1 / a} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, and \eqn{\theta > 0}, the second shape parameter.}
\usage{
dgenpowerweibull(x, a=1, theta=1, log=FALSE)
pgenpowerweibull(x, a=1, theta=1, log.p=FALSE, lower.tail=TRUE)
vargenpowerweibull(p, a=1, theta=1, log.p=FALSE, lower.tail=TRUE)
esgenpowerweibull(p, a=1, theta=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{theta}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgenpowerweibull(x)
pgenpowerweibull(x)
vargenpowerweibull(x)
esgenpowerweibull(x)}
 
