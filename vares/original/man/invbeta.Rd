\name{invbeta}
\alias{dinvbeta}
\alias{pinvbeta}
\alias{varinvbeta}
\alias{esinvbeta}
\title{Inverse beta distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the inverse beta distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {x^{a - 1}}{B (a, b) (1 + x)^{a + b}},
\\
&\displaystyle
F (x) = I_{\frac {x}{1 + x}} (a, b),
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {I_p^{-1} (a, b)}{1 - I_p^{-1} (a, b)},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p \frac {I_v^{-1} (a, b)}{1 - I_v^{-1} (a, b)} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dinvbeta(x, a=1, b=1, log=FALSE)
pinvbeta(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varinvbeta(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esinvbeta(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
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
dinvbeta(x)
pinvbeta(x)
varinvbeta(x)
esinvbeta(x)}
 

