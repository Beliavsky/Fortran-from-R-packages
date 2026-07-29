\name{expweibull}
\alias{dexpweibull}
\alias{pexpweibull}
\alias{varexpweibull}
\alias{esexpweibull}
\title{Exponentiated Weibull distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponentiated Weibull distribution due to Mudholkar and Srivastava (1993) and Mudholkar et al. (1995) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = a \alpha \sigma^{-\alpha} x^{\alpha - 1}
\exp \left[ -(x / \sigma)^\alpha \right]
\left\{ 1 - \exp \left[ -(x / \sigma)^\alpha \right] \right\}^{a - 1},
\\
&\displaystyle
F (x) = \left\{ 1 - \exp \left[ -(x / \sigma)^\alpha \right] \right\}^a,
\\
&\displaystyle
{\rm VaR}_p (X) = \sigma \left[ -\log \left( 1 - p^{1 / a} \right) \right]^{1 / \alpha},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sigma}{p} \int_0^p \left[ -\log \left( 1 - v^{1 / a} \right) \right]^{1 / \alpha} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter,
\eqn{\alpha > 0}, the second shape parameter, and \eqn{\sigma > 0}, the scale parameter.}
\usage{
dexpweibull(x, a=1, alpha=1, sigma=1, log=FALSE)
pexpweibull(x, a=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
varexpweibull(p, a=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
esexpweibull(p, a=1, alpha=1, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{alpha}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexpweibull(x)
pexpweibull(x)
varexpweibull(x)
esexpweibull(x)}
 
