\name{lomax}
\alias{dlomax}
\alias{plomax}
\alias{varlomax}
\alias{eslomax}
\title{Lomax distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Lomax distribution due to Lomax (1954) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {a}{\lambda} \left( 1 + \frac {x}{\lambda} \right)^{-a - 1},
\\
&\displaystyle
F (x) = 1 - \left( 1 + \frac {x}{\lambda} \right)^{-a},
\\
&\displaystyle
{\rm VaR}_p (X) = \lambda \left[ (1 - p)^{-1 / a} - 1 \right],
\\
&\displaystyle
{\rm ES}_p (X) =  -\lambda + \frac {\lambda - \lambda (1 - p)^{1 - 1 / a}}{p - p / a}
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter, and \eqn{\lambda > 0}, the scale parameter.}
\usage{
dlomax(x, a=1, lambda=1, log=FALSE)
plomax(x, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
varlomax(p, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
eslomax(p, a=1, lambda=1)
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
dlomax(x)
plomax(x)
varlomax(x)
eslomax(x)}
 
