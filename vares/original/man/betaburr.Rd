\name{betaburr}
\alias{dbetaburr}
\alias{pbetaburr}
\alias{varbetaburr}
\alias{esbetaburr}
\title{Beta Burr distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the beta Burr distribution due to Parana\'{i}ba et al. (2011) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {b a^{bd}}{B (c, d)x^{bd + 1}} \left[ 1 + \left( x / a \right)^{-b} \right]^{-c - d},
\\
&\displaystyle
F (x) = I_{\frac {1}{1 + \left( x / a \right)^{-b}}} (c, d),
\\
&\displaystyle
{\rm VaR}_p (X) = a \left[ I_p^{-1} (c, d) \right]^{1 / b} \left[ 1 - I_p^{-1} (c, d) \right]^{-1 / b},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {a}{p} \int_0^p \left[ I_v^{-1} (c, d) \right]^{1 / b}
\left[ 1 - I_v^{-1} (c, d) \right]^{-1 / b} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the scale parameter, \eqn{b > 0}, the first shape parameter,
\eqn{c > 0}, the second shape parameter, and \eqn{d > 0}, the third shape parameter,
where \eqn{I_x (a, b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt / B (a, b)} denotes the incomplete beta function ratio,
\eqn{B (a, b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt} denotes the beta function,
and \eqn{I_x^{-1} (a, b)} denotes the inverse function of \eqn{I_x (a, b)}.}
\usage{
dbetaburr(x, a=1, b=1, c=1, d=1, log=FALSE)
pbetaburr(x, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
varbetaburr(p, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
esbetaburr(p, a=1, b=1, c=1, d=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the scale parameter, must be positive, the default is 1}
  \item{b}{the value of the first shape parameter, must be positive, the default is 1}
  \item{c}{the value of the second shape parameter, must be positive, the default is 1}
  \item{d}{the value of the third shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbetaburr(x)
pbetaburr(x)
varbetaburr(x)
esbetaburr(x)}
 
