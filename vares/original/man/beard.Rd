\name{beard}
\alias{dbeard}
\alias{pbeard}
\alias{varbeard}
\alias{esbeard}
\title{Beard distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Beard distribution due to Beard (1959) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {\displaystyle a \exp (b x) \left[ 1 + a \rho \right]^{\rho^{-1/b}}}
{\displaystyle \left[ 1 + a \rho \exp (b x) \right]^{1 + \rho^{-1/b}}},
\\
&\displaystyle
F (x) = 1 - \frac {\displaystyle \left[ 1 + a \rho \right]^{\rho^{-1/b}}}
{\displaystyle \left[ 1 + a \rho \exp (b x) \right]^{\rho^{-1/b}}},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{b} \log \left[ \frac {1 + a \rho}{a \rho (1 - p)^{\rho^{1 / b}}} - \frac {1}{a \rho} \right],
\\
&\displaystyle
{\rm ES}_p (X) =
\frac {1}{p b} \int_0^p \log \left[ -\frac {1}{a \rho} +
\frac {1 + a \rho}{a \rho (1 - v)^{\rho^{1 / b}}} \right] dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first scale parameter, \eqn{b > 0}, the second scale parameter,
and \eqn{\rho > 0}, the shape parameter.}
\usage{
dbeard(x, a=1, b=1, rho=1, log=FALSE)
pbeard(x, a=1, b=1, rho=1, log.p=FALSE, lower.tail=TRUE)
varbeard(p, a=1, b=1, rho=1, log.p=FALSE, lower.tail=TRUE)
esbeard(p, a=1, b=1, rho=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{b}{the value of the second scale parameter, must be positive, the default is 1}
  \item{rho}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbeard(x)
pbeard(x)
varbeard(x)
esbeard(x)}
 
