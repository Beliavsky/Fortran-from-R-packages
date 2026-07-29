\name{burr}
\alias{dburr}
\alias{pburr}
\alias{varburr}
\alias{esburr}
\title{Burr distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Burr distribution due to Burr (1942) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {b a^b}{x^{b + 1}} \left[ 1 + \left( x / a \right)^{-b} \right]^{-2},
\\
&\displaystyle
F (x) = \frac {1}{1 + \left( x / a \right)^{-b}},
\\
&\displaystyle
{\rm VaR}_p (X) = a p^{1 / b} (1 - p)^{-1 / b},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {a}{p} B_p \left( 1 / b + 1, 1 - 1 / b \right)
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the scale parameter, and \eqn{b > 0}, the shape parameter,
where  \eqn{B_x (a, b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt} denotes the incomplete beta function.}
\usage{
dburr(x, a=1, b=1, log=FALSE)
pburr(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varburr(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esburr(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the scale parameter, must be positive, the default is 1}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dburr(x)
pburr(x)
varburr(x)
esburr(x)}
 
