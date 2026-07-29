\name{betaexp}
\alias{dbetaexp}
\alias{pbetaexp}
\alias{varbetaexp}
\alias{esbetaexp}
\title{Beta exponential distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the beta exponential distribution due to Nadarajah and Kotz (2006) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {\lambda \exp (-b \lambda x)}{B (a, b)}
\left[ 1 - \exp (-\lambda x) \right]^{a - 1},
\\
&\displaystyle
F (x) = I_{1 - \exp (-\lambda x)} (a, b),
\\
&\displaystyle
{\rm VaR}_p (X) = -\frac {1}{\lambda} \log \left[ 1 - I_p^{-1} (a, b) \right],
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {1}{p \lambda} \int_0^p \log \left[ 1 - I_v^{-1} (a, b) \right] dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, \eqn{b > 0}, the second shape parameter, and \eqn{\lambda > 0},
the scale parameter, where \eqn{I_x (a, b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt / B (a, b)} denotes
the incomplete beta function ratio,
\eqn{B (a, b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt} denotes the beta function, and
\eqn{I_x^{-1} (a, b)} denotes the inverse function of \eqn{I_x (a, b)}.}
\usage{
dbetaexp(x, lambda=1, a=1, b=1, log=FALSE)
pbetaexp(x, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varbetaexp(p, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esbetaexp(p, lambda=1, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{lambda}{the value of the scale parameter, must be positive, the default is 1}
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
dbetaexp(x)
pbetaexp(x)
varbetaexp(x)
esbetaexp(x)}
 
