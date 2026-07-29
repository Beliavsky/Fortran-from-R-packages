\name{logistic}
\alias{dlogistic}
\alias{plogistic}
\alias{varlogistic}
\alias{eslogistic}
\title{Logistic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the logistic distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\sigma} \exp \left( -\frac {x - \mu}{\sigma} \right)
\left[ 1 + \exp \left( -\frac {x - \mu}{\sigma} \right) \right]^{-2},
\\
&\displaystyle
F (x) = \frac {1}{1 + \exp \left( -\frac {x - \mu}{\sigma} \right)},
\\
&\displaystyle
{\rm VaR}_p (X) = \mu + \sigma \log \left[ p (1 -  p) \right],
\\
&\displaystyle
{\rm ES}_p (X) =  \mu - 2 \sigma + \sigma \log p - \sigma \frac {1 - p}{p} \log (1 - p)
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, and \eqn{\sigma > 0}, the scale parameter.}
\usage{
dlogistic(x, mu=0, sigma=1, log=FALSE)
plogistic(x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
varlogistic(p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
eslogistic(p, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dlogistic(x)
plogistic(x)
varlogistic(x)
eslogistic(x)}
 
