\name{xie}
\alias{dxie}
\alias{pxie}
\alias{varxie}
\alias{esxie}
\title{Xie distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Xie distribution due to Xie et al. (2002) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \lambda b \left( \frac {x}{a} \right)^{b - 1}
\exp \left[ (x/a)^b \right] \exp \left( \lambda a  \right)
\exp \left\{ -\lambda a \exp \left[ (x/a)^b \right] \right\},
\\
&\displaystyle
F (x) =
1 - \exp \left( \lambda a \right)
\exp \left\{ -\lambda a \exp \left[ (x/a)^b \right] \right\},
\\
&\displaystyle
{\rm VaR}_p (X) =
a \left\{ \log \left[ 1 - \frac {\log (1 - p)}{\lambda a} \right] \right\}^{1/b},
\\
&\displaystyle
{\rm ES}_p (X) =
\frac {a}{p} \int_0^p \left\{ \log \left[ 1 - \frac {\log (1 - v)}{\lambda a} \right] \right\}^{1/b} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first scale parameter, \eqn{b > 0}, the shape parameter,
and \eqn{\lambda > 0}, the second scale parameter.}
\usage{
dxie(x, a=1, b=1, lambda=1, log=FALSE)
pxie(x, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
varxie(p, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
esxie(p, a=1, b=1, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{lambda}{the value of the second scale parameter, must be positive, the default is 1}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dxie(x)
pxie(x)
varxie(x)
esxie(x)}
 
