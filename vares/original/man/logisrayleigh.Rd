\name{logisrayleigh}
\alias{dlogisrayleigh}
\alias{plogisrayleigh}
\alias{varlogisrayleigh}
\alias{eslogisrayleigh}
\title{Logistic Rayleigh distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the  logistic Rayleigh distribution due to Lan and Leemis (2008) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = a \lambda x \exp \left( \lambda x^2 / 2 \right)
\left[ \exp \left( \lambda x^2 / 2 \right) - 1 \right]^{a - 1}
\left\{ 1 + \left[ \exp \left( \lambda x^2 / 2 \right) - 1 \right]^a \right\}^{-2},
\\
&\displaystyle
F(x) = \frac {\left[ \exp \left( \lambda x^2 / 2 \right) - 1 \right]^a}
{1 + \left[ \exp \left( \lambda x^2 / 2 \right) - 1 \right]^a},
\\
&\displaystyle
{\rm VaR}_p (X) = \sqrt{\frac {2}{\lambda}}
\sqrt{\log \left[ 1 + \left( \frac {p}{1 - p} \right)^{1 / a} \right]},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sqrt{2}}{p \sqrt{\lambda}}
\int_0^p \left\{ \log \left[ 1 + \left( \frac {v}{1 - v} \right)^{1 / a} \right] \right\}^{1/2} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter, and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dlogisrayleigh(x, a=1, lambda=1, log=FALSE)
plogisrayleigh(x, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
varlogisrayleigh(p, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
eslogisrayleigh(p, a=1, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dlogisrayleigh(x)
plogisrayleigh(x)
varlogisrayleigh(x)
eslogisrayleigh(x)}
 
