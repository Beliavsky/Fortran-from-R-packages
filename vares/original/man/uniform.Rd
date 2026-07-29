\name{uniform}
\alias{duniform}
\alias{puniform}
\alias{varuniform}
\alias{esuniform}
\title{Uniform distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the uniform distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{b - a},
\\
&\displaystyle
F (x) = \frac {x - a}{b - a},
\\
&\displaystyle
{\rm VaR}_p (X) = a + p (b - a),
\\
&\displaystyle
{\rm ES}_p (X) = a + \frac {p}{2} (b - a)
\end{array}}
for \eqn{a < x < b}, \eqn{0 < p < 1}, \eqn{-\infty < a < \infty} , the first location parameter, and \eqn{-\infty < a < b < \infty} , the second location parameter.}
\usage{
duniform(x, a=0, b=1, log=FALSE)
puniform(x, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
varuniform(p, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
esuniform(p, a=0, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first location parameter, can take any real value, the default is zero}
  \item{b}{the value of the second location parameter, can take any real value but must be greater than a, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
duniform(x)
puniform(x)
varuniform(x)
esuniform(x)}
 
