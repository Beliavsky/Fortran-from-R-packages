\name{paretostable}
\alias{dparetostable}
\alias{pparetostable}
\alias{varparetostable}
\alias{esparetostable}
\title{Pareto positive stable distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Pareto positive stable distribution due to Sarabia and Prieto (2009) and Guillen et al. (2011) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {\nu \lambda}{x}
\left[ \log \left( \frac {x}{\sigma} \right) \right]^{\nu - 1}
\exp \left\{ -\lambda \left[ \log \left( \frac {x}{\sigma} \right) \right]^\nu \right\},
\\
&\displaystyle
F (x) = 1 - \exp \left\{ -\lambda \left[ \log \left( \frac {x}{\sigma} \right) \right]^\nu \right\},
\\
&\displaystyle
{\rm VaR}_p (X) = \sigma \exp \left\{ \left[ -\frac {1}{\lambda} \log (1 - p) \right]^{1/\nu} \right\},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sigma}{p} \int_0^p \exp \left\{ \left[ -\frac {1}{\lambda} \log (1 - v) \right]^{1/\nu} \right\} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{\lambda > 0}, the first scale parameter, \eqn{\sigma > 0}, the second scale parameter, and \eqn{\nu > 0}, the shape parameter.}
\usage{
dparetostable(x, lambda=1, nu=1, sigma=1, log=FALSE)
pparetostable(x, lambda=1, nu=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
varparetostable(p, lambda=1, nu=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
esparetostable(p, lambda=1, nu=1, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the first scale parameter, must be positive, the default is 1}
  \item{sigma}{the value of the second scale parameter, must be positive, the default is 1}
  \item{nu}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dparetostable(x)
pparetostable(x)
varparetostable(x)
esparetostable(x)}
 
