\name{nakagami}
\alias{dnakagami}
\alias{pnakagami}
\alias{varnakagami}
\alias{esnakagami}
\title{Nakagami distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Nakagami distribution due to Nakagami (1960) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {2 m^m}{\Gamma (m) a^m} x^{2 m - 1}
\exp \left( -\frac {m x^2}{a} \right),
\\
&\displaystyle
F (x) = 1 - Q \left( m, \frac {m x^2}{a} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \sqrt{\frac {a}{m}} \sqrt{Q^{-1} (m, 1 - p)},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sqrt{a}}{p \sqrt{m}} \int_0^p \sqrt{Q^{-1} (m, 1 - v)} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the scale parameter, and \eqn{m > 0}, the shape parameter.}
\usage{
dnakagami(x, m=1, a=1, log=FALSE)
pnakagami(x, m=1, a=1, log.p=FALSE, lower.tail=TRUE)
varnakagami(p, m=1, a=1, log.p=FALSE, lower.tail=TRUE)
esnakagami(p, m=1, a=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the scale parameter, must be positive, the default is 1}
  \item{m}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dnakagami(x)
pnakagami(x)
varnakagami(x)
esnakagami(x)}
 
