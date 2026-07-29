\name{lfr}
\alias{dlfr}
\alias{plfr}
\alias{varlfr}
\alias{eslfr}
\title{Linear failure rate distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the linear failure rate distribution due to Bain (1974) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) =  (a + b x) \exp \left( -a x - b x^2 / 2 \right),
\\
&\displaystyle
F (x) = 1 - \exp \left( -a x - b x^2 / 2 \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {-a + \sqrt{a^2 - 2 b \log (1 - p)}}{b},
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {a}{b} + \frac {1}{b p} \int_0^p \sqrt{a^2 - 2 b \log (1 - v)} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first scale parameter, and \eqn{b > 0}, the second scale parameter.}
\usage{
dlfr(x, a=1, b=1, log=FALSE)
plfr(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varlfr(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
eslfr(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{b}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dlfr(x)
plfr(x)
varlfr(x)
eslfr(x)}
 
