\name{explog}
\alias{dexplog}
\alias{pexplog}
\alias{varexplog}
\alias{esexplog}
\title{Exponential logarithmic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponential logarithmic distribution due to Tahmasbi and Rezaei (2008) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = -\frac {b (1 - a) \exp (-b x)}{\log a \left[ 1 - (1 - a) \exp (-b x) \right]},
\\
&\displaystyle
F (x) = 1 - \frac {\log  \left[ 1 - (1 - a) \exp (-b x) \right]}{\log a},
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{b} \log \left[ \frac {1 - a^{1 - p}}{1 - a} \right],
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {1}{b p} \int_0^p \log \left[ \frac {1 - a^{1 - v}}{1 - a} \right] dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{0 < a < 1}, the first scale parameter, and \eqn{b > 0}, the second scale parameter.}
\usage{
dexplog(x, a=0.5, b=1, log=FALSE)
pexplog(x, a=0.5, b=1, log.p=FALSE, lower.tail=TRUE)
varexplog(p, a=0.5, b=1, log.p=FALSE, lower.tail=TRUE)
esexplog(p, a=0.5, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be in the unit interval, the default is 0.5}
  \item{b}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexplog(x)
pexplog(x)
varexplog(x)
esexplog(x)}
 
