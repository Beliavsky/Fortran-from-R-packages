\name{kumexp}
\alias{dkumexp}
\alias{pkumexp}
\alias{varkumexp}
\alias{eskumexp}
\title{Kumaraswamy exponential distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Kumaraswamy exponential distribution due to Cordeiro and de Castro (2011) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = a b \lambda \exp (-\lambda x) \left[ 1 - \exp (-\lambda x) \right]^{a - 1}
\left\{ 1 - \left[ 1 - \exp (-\lambda x) \right]^a \right\}^{b - 1},
\\
&\displaystyle
F (x) = 1 - \left\{ 1 - \left[ 1 - \exp (-\lambda x) \right]^a \right\}^b,
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{\lambda} \log \left\{ 1 - \left[ 1 - (1 - p)^{1 / b} \right]^{1 / a} \right\},
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {1}{p \lambda} \int_0^p \log \left\{ 1 - \left[ 1 - (1 - v)^{1 / b} \right]^{1 / a} \right\} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, \eqn{b > 0}, the second shape parameter, and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dkumexp(x, lambda=1, a=1, b=1, log=FALSE)
pkumexp(x, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varkumexp(p, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
eskumexp(p, lambda=1, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dkumexp(x)
pkumexp(x)
varkumexp(x)
eskumexp(x)}
 
