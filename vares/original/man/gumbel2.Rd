\name{gumbel2}
\alias{dgumbel2}
\alias{pgumbel2}
\alias{vargumbel2}
\alias{esgumbel2}
\title{Gumbel II distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Gumbel II distribution
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = a b x^{-a - 1} \exp \left( -b x^{-a} \right),
\\
&\displaystyle
F (x) = 1 - \exp \left( -b x^{-a} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = b^{1 / a} \left[ -\log (1 - p) \right]^{-1 / a},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {b^{1 / a}}{p} \int_0^p \left[ -\log (1 - v) \right]^{-1 / a} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter, and \eqn{b > 0}, the scale parameter.}
\usage{
dgumbel2(x, a=1, b=1, log=FALSE)
pgumbel2(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
vargumbel2(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esgumbel2(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the scale parameter, must be positive, the default is 1}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgumbel2(x)
pgumbel2(x)
vargumbel2(x)
esgumbel2(x)
}
 
