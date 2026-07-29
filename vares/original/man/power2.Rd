\name{power2}
\alias{dpower2}
\alias{ppower2}
\alias{varpower2}
\alias{espower2}
\title{Power function II distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the power function II distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = b (1 - x)^{b - 1},
\\
&\displaystyle
F (x) = 1 - (1 - x)^b,
\\
&\displaystyle
{\rm VaR}_p (X) = 1 - (1 - p)^{1 / b},
\\
&\displaystyle
{\rm ES}_p (X) = 1 + \frac {b \left[ (1 - p)^{1 / b + 1} - 1 \right]}{p (b + 1)}
\end{array}}
for \eqn{0 < x < 1}, \eqn{0 < p < 1},  and \eqn{b > 0}, the shape parameter.}
\usage{
dpower2(x, b=1, log=FALSE)
ppower2(x, b=1, log.p=FALSE, lower.tail=TRUE)
varpower2(p, b=1, log.p=FALSE, lower.tail=TRUE)
espower2(p, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dpower2(x)
ppower2(x)
varpower2(x)
espower2(x)}
 
