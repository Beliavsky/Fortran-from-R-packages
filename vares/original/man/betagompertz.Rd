\name{betagompertz}
\alias{dbetagompertz}
\alias{pbetagompertz}
\alias{varbetagompertz}
\alias{esbetagompertz}
\title{Beta Gompertz distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the beta Gompertz distribution due to Cordeiro et al. (2012b) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) =  \frac {b \eta \exp (bx)}{B (c, d)}
\exp \left( d \eta \right)
\exp \left[ -d \eta \exp (bx) \right]
\left\{ 1 - \exp \left[ \eta - \eta \exp (bx) \right] \right\}^{c - 1},
\\
&\displaystyle
F(x) = I_{1 - \exp \left[ \eta - \eta \exp (bx) \right]} (c, d),
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{b} \log \left\{ 1 - \frac {1}{\eta}
\log \left[ 1 - I_p^{-1} (c, d) \right] \right\},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p b} \int_0^p \log \left\{ 1 - \frac {1}{\eta} \log
\left[ 1 - I_v^{-1} (c, d) \right] \right\} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{b > 0}, the first scale parameter, \eqn{\eta > 0}, the second scale parameter,
\eqn{c > 0}, the first shape parameter, and \eqn{d > 0}, the second shape parameter.}
\usage{
dbetagompertz(x, b=1, c=1, d=1, eta=1, log=FALSE)
pbetagompertz(x, b=1, c=1, d=1, eta=1, log.p=FALSE, lower.tail=TRUE)
varbetagompertz(p, b=1, c=1, d=1, eta=1, log.p=FALSE, lower.tail=TRUE)
esbetagompertz(p, b=1, c=1, d=1, eta=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{b}{the value of the first scale parameter, must be positive, the default is 1}
  \item{eta}{the value of the second scale parameter, must be positive, the default is 1}
  \item{c}{the value of the first shape parameter, must be positive, the default is 1}
  \item{d}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbetagompertz(x)
pbetagompertz(x)
varbetagompertz(x)
esbetagompertz(x)}
 
