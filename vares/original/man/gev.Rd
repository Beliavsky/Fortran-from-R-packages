

\name{gev}
\alias{dgev}
\alias{pgev}
\alias{vargev}
\alias{esgev}
\title{Generalized extreme value distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the generalized extreme value distribution due to Fisher and Tippett (1928) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) =
\frac {1}{\sigma} \left[ 1 + \xi \left( \frac {x - \mu}{\sigma} \right) \right]^{-1/\xi - 1}
\exp \left\{ -\left[ 1 + \xi \left( \frac {x - \mu}{\sigma} \right)
\right]^{-1/\xi} \right\},
\\
&\displaystyle
F(x) = \exp \left\{ -\left[ 1 + \xi \left( \frac {x - \mu}{\sigma} \right) \right]^{-1/\xi} \right\},
\\
&\displaystyle
{\rm VaR}_p (X) = \mu - \frac {\sigma}{\xi} + \frac {\sigma}{\xi} (-\log p)^{-\xi},
\\
&\displaystyle
{\rm ES}_p (X) = \mu - \frac {\sigma}{\xi} + \frac {\sigma}{p \xi} \int_0^p (-\log v)^{-\xi} dv
\end{array}}
for \eqn{x \geq \mu - \sigma / \xi}  if \eqn{\xi > 0},
\eqn{x \leq \mu - \sigma / \xi} if \eqn{\xi < 0},
\eqn{-\infty < x < \infty} if \eqn{\xi = 0},
\eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter,
\eqn{\sigma > 0}, the scale parameter, and \eqn{-\infty < \xi < \infty}, the shape parameter.}
\usage{
dgev(x, mu=0, sigma=1, xi=1, log=FALSE)
pgev(x, mu=0, sigma=1, xi=1, log.p=FALSE, lower.tail=TRUE)
vargev(p, mu=0, sigma=1, xi=1, log.p=FALSE, lower.tail=TRUE)
esgev(p, mu=0, sigma=1, xi=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{xi}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgev(x)
pgev(x)
vargev(x)
esgev(x)}
 


