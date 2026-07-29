\name{kumnormal}
\alias{dkumnormal}
\alias{pkumnormal}
\alias{varkumnormal}
\alias{eskumnormal}
\title{Kumaraswamy normal distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for Kumaraswamy normal distribution due to Cordeiro and de Castro (2011) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {a b}{\sigma} \phi \left( \frac {x - \mu}{\sigma} \right) \Phi^{a - 1} \left( \frac {x - \mu}{\sigma} \right)
\left[ 1 - \Phi^a \left( \frac {x - \mu}{\sigma} \right) \right]^{b - 1},
\\
&\displaystyle
F (x) = 1 - \left[ 1 - \Phi^a \left( \frac {x - \mu}{\sigma} \right) \right]^b,
\\
&\displaystyle
{\rm VaR}_p (X) = \mu + \sigma \Phi^{-1} \left( \left[ 1 - (1 - p)^{1 / b} \right]^{1 / a} \right),
\\
&\displaystyle
{\rm ES}_p (X) = \mu + \frac {\sigma}{p} \int_0^p \Phi^{-1} \left( \left[ 1 - (1 - v)^{1 / b} \right]^{1 / a} \right) dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, \eqn{\sigma > 0}, the scale parameter,
\eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dkumnormal(x, mu=0, sigma=1, a=1, b=1, log=FALSE)
pkumnormal(x, mu=0, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varkumnormal(p, mu=0, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
eskumnormal(p, mu=0, sigma=1, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dkumnormal(x)
pkumnormal(x)
varkumnormal(x)
eskumnormal(x)}
 
