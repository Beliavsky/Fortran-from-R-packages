\name{moexp}
\alias{dmoexp}
\alias{pmoexp}
\alias{varmoexp}
\alias{esmoexp}
\title{Marshall-Olkin exponential distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Marshall-Olkin exponential distribution due to Marshall and Olkin (1997) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {\displaystyle \lambda \exp (\lambda x)}
{\displaystyle \left[ \exp (\lambda x) - 1 + a \right]^2},
\\
&\displaystyle
F (x) =  \frac {\displaystyle \exp (\lambda x) - 2 + a}{\displaystyle \exp (\lambda x) - 1 + a},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{\lambda} \log \frac {2 - a - (1 - a) p}{1 - p},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{\lambda} \log
\left[ 2 - a - (1 - a) p \right] - \frac {2 - a}{\lambda (1 - a) p}
\log \frac {2 - a - (1 - a) p}{2 - a} + \frac {1 - p}{\lambda p} \log (1 - p)
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first scale parameter and \eqn{\lambda > 0}, the second scale parameter.}
\usage{
dmoexp(x, lambda=1, a=1, log=FALSE)
pmoexp(x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
varmoexp(p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
esmoexp(p, lambda=1, a=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{lambda}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dmoexp(x)
pmoexp(x)
varmoexp(x)
esmoexp(x)}
 
