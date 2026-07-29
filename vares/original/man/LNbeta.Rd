\name{LNbeta}
\alias{dLNbeta}
\alias{pLNbeta}
\alias{varLNbeta}
\alias{esLNbeta}
\title{Libby-Novick beta distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Libby-Novick beta distribution due to Libby and Novick (1982) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {\lambda^a x^{a - 1} (1 - x)^{b - 1}}
{B (a, b) \left[ 1 - (1 - \lambda) x \right]^{a + b}},
\\
&\displaystyle
F (x) = I_{\frac {\lambda x}{1 + (\lambda - 1) x}} (a, b),
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {I_p^{-1} (a, b)}{\lambda - (\lambda - 1) I_p^{-1} (a, b)},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p
\frac {I_v^{-1} (a, b)}{\lambda - (\lambda - 1) I_v^{-1} (a, b)} dv
\end{array}}
for \eqn{0 < x < 1}, \eqn{0 < p < 1}, \eqn{\lambda > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dLNbeta(x, lambda=1, a=1, b=1, log=FALSE)
pLNbeta(x, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varLNbeta(p, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esLNbeta(p, lambda=1, a=1, b=1)
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
dLNbeta(x)
pLNbeta(x)
varLNbeta(x)
esLNbeta(x)}
 
