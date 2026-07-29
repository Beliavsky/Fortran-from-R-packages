\name{RS}
\alias{varRS}
\alias{esRS}
\title{Ramberg-Schmeiser distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Ramber-Schmeiser distribution due to Ramberg and Schmeiser (1974) given by
\deqn{\begin{array}{ll}
&\displaystyle
{\rm VaR}_p (X) = \frac {p^b - (1 - p)^c}{d},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {p^{b}}{d (b + 1)} + \frac {(1 - p)^{c + 1} - 1}{p d (c + 1)}
\end{array}}
for \eqn{0 < p < 1}, \eqn{b > 0}, the first shape parameter,
\eqn{c > 0}, the second shape parameter, and \eqn{d > 0}, the scale parameter.}
\usage{
varRS(p, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
esRS(p, b=1, c=1, d=1)
}
\arguments{
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{d}{the value of the scale parameter, must be positive, the default is 1}
  \item{b}{the value of the first shape parameter, must be positive, the default is 1}
  \item{c}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
varRS(x)
esRS(x)}
 

