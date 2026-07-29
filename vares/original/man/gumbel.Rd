\name{gumbel}
\alias{dgumbel}
\alias{pgumbel}
\alias{vargumbel}
\alias{esgumbel}
\title{Gumbel distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Gumbel distribution given by  due to Gumbel (1954) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\sigma} \exp \left( \frac {\mu - x}{\sigma} \right)
\exp \left[ -\exp \left( \frac {\mu - x}{\sigma} \right) \right],
\\
&\displaystyle
F (x) = \exp \left[ -\exp \left( \frac {\mu - x}{\sigma} \right) \right],
\\
&\displaystyle
{\rm VaR}_p (X) = \mu - \sigma \log (-\log p),
\\
&\displaystyle
{\rm ES}_p (X) = \mu - \frac {\sigma}{p} \int_0^p \log (-\log v) dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, and \eqn{\sigma > 0}, the scale parameter.}
\usage{
dgumbel(x, mu=0, sigma=1, log=FALSE)
pgumbel(x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
vargumbel(p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
esgumbel(p, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgumbel(x)
pgumbel(x)
vargumbel(x)
esgumbel(x)}
 
