\name{stacygamma}
\alias{dstacygamma}
\alias{pstacygamma}
\alias{varstacygamma}
\alias{esstacygamma}
\title{Stacy distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for Stacy distribution due to Stacy (1962) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {c x^{c \gamma - 1} \exp \left[ -(x / \theta)^c \right]}{\theta^{c \gamma} \Gamma (\gamma)},
\\
&\displaystyle
F (x) = 1 - Q \left( \gamma, \left( \frac {x}{\theta} \right)^c \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \theta \left[ Q^{-1} (\gamma, 1 - p) \right]^{1 / c},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\theta}{p} \int_0^p \left[ Q^{-1} (\gamma, 1 - v) \right]^{1 / c} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{\theta > 0}, the scale parameter, \eqn{c > 0}, the first shape parameter, and
\eqn{\gamma > 0}, the second shape parameter.}
\usage{
dstacygamma(x, gamma=1, c=1, theta=1, log=FALSE)
pstacygamma(x, gamma=1, c=1, theta=1, log.p=FALSE, lower.tail=TRUE)
varstacygamma(p, gamma=1, c=1, theta=1, log.p=FALSE, lower.tail=TRUE)
esstacygamma(p, gamma=1, c=1, theta=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the scale parameter, must be positive, the default is 1}
  \item{c}{the value of the first scale parameter, must be positive, the default is 1}
  \item{gamma}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dstacygamma(x)
pstacygamma(x)
varstacygamma(x)
esstacygamma(x)}
 
