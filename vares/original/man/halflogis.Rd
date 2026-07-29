\name{halflogis}
\alias{dhalflogis}
\alias{phalflogis}
\alias{varhalflogis}
\alias{eshalflogis}
\title{Half logistic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the half logistic distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {2 \lambda \exp \left( -\lambda x \right)}
{\left[ 1 + \exp \left( -\lambda x \right) \right]^2},
\\
&\displaystyle
F (x) = \frac {1 - \exp \left( -\lambda x \right)}{1 + \exp \left( -\lambda x \right)},
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{\lambda} \log \frac {1 - p}{1 + p},
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {1}{\lambda} \log \frac {1 - p}{1 + p} + \frac {1}{\lambda p} \log \left( 1 - p^2 \right)
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dhalflogis(x, lambda=1, log=FALSE)
phalflogis(x, lambda=1, log.p=FALSE, lower.tail=TRUE)
varhalflogis(p, lambda=1, log.p=FALSE, lower.tail=TRUE)
eshalflogis(p, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dhalflogis(x)
phalflogis(x)
varhalflogis(x)
eshalflogis(x)}
 
