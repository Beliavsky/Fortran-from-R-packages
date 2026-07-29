\name{genunif}
\alias{dgenunif}
\alias{pgenunif}
\alias{vargenunif}
\alias{esgenunif}
\title{Generalized uniform distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the generalized uniform distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = h k c (x - a)^{c - 1} \left[ 1 - k (x - a)^c \right]^{h - 1},
\\
&\displaystyle
F (x) = 1 - \left[ 1 - k (x - a)^c \right]^h,
\\
&\displaystyle
{\rm VaR}_p (X) = a + k^{-1/c} \left[ 1 - (1 - p)^{1/h} \right]^{1/c},
\\
&\displaystyle
{\rm ES}_p (X) = a + \frac {k^{-1/c}}{p} \int_0^p \left[ 1 -  (1 - v)^{1/h} \right]^{1/c} dv
\end{array}}
for \eqn{a \leq x \leq a + k^{-1/c}}, \eqn{0 < p < 1}, \eqn{-\infty < a < \infty}, the location parameter, \eqn{c > 0}, the first shape parameter, \eqn{k > 0}, the scale parameter, and \eqn{h > 0}, the second shape parameter.}
\usage{
dgenunif(x, a=0, c=1, h=1, k=1, log=FALSE)
pgenunif(x, a=0, c=1, h=1, k=1, log.p=FALSE, lower.tail=TRUE)
vargenunif(p, a=0, c=1, h=1, k=1, log.p=FALSE, lower.tail=TRUE)
esgenunif(p, a=0, c=1, h=1, k=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the location parameter, can take any real value, the default is zero}
  \item{k}{the value of the scale parameter, must be positive, the default is 1}
  \item{c}{the value of the first scale parameter, must be positive, the default is 1}
  \item{h}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgenunif(x)
pgenunif(x)
vargenunif(x)
esgenunif(x)}
 
