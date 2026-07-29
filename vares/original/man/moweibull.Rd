\name{moweibull}
\alias{dmoweibull}
\alias{pmoweibull}
\alias{varmoweibull}
\alias{esmoweibull}
\title{Marshall-Olkin Weibull distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Marshall-Olkin Weibull distribution due to Marshall and Olkin (1997) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = b \lambda^b x^{b -  1} \exp \left[ (\lambda x)^b \right]
\left\{ \exp \left[ (\lambda x)^b \right] - 1 + a \right\}^{-2},
\\
&\displaystyle
F(x) = \frac {\displaystyle \exp \left[ (\lambda x)^b \right] - 2 + a}
{\displaystyle \exp \left[ (\lambda x)^b \right] - 1 + a},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{\lambda} \left[ \log \left( \frac {1}{1 - p} + 1 - a \right) \right]^{1 / b},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{\lambda p} \int_0^p \left[ \log \left( \frac {1}{1 - v} + 1 - a \right) \right]^{1 / b} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first scale parameter, \eqn{b > 0}, the shape parameter,
and \eqn{\lambda > 0}, the second scale parameter.}
\usage{
dmoweibull(x, a=1, b=1, lambda=1, log=FALSE)
pmoweibull(x, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
varmoweibull(p, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
esmoweibull(p, a=1, b=1, lambda=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{lambda}{the value of the second scale parameter, must be positive, the default is 1}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dmoweibull(x)
pmoweibull(x)
varmoweibull(x)
esmoweibull(x)}
 
