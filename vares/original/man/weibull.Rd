\name{weibull}
\alias{dWeibull}
\alias{pWeibull}
\alias{varWeibull}
\alias{esWeibull}
\title{Weibull distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Weibull distribution due to Weibull (1951) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {\alpha x^{\alpha - 1}}{\sigma^\alpha}
\exp \left\{ -\left( \frac {x}{\sigma} \right)^{\alpha} \right\},
\\
&\displaystyle
F (x) =  1 - \exp \left\{ -\left( \frac {x}{\sigma} \right)^{\alpha} \right\},
\\
&\displaystyle
{\rm VaR}_p (X) = \sigma \left[ -\log (1 - p) \right]^{1 / \alpha},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sigma}{p} \gamma \left( 1 + 1 / \alpha, - \log (1 - p) \right)
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{\alpha > 0}, the shape parameter, and \eqn{\sigma > 0}, the scale parameter.}
\usage{
dWeibull(x, alpha=1, sigma=1, log=FALSE)
pWeibull(x, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
varWeibull(p, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
esWeibull(p, alpha=1, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{alpha}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dWeibull(x)
pWeibull(x)
varWeibull(x)
esWeibull(x)}
 

