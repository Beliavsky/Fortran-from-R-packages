\name{perks}
\alias{dperks}
\alias{pperks}
\alias{varperks}
\alias{esperks}
\title{Perks distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Perks distribution due to Perks (1932) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {\displaystyle a \exp (b x) \left[ 1 + a \right]}
{\displaystyle \left[ 1 + a \exp (b x) \right]^2},
\\
&\displaystyle
F (x) = 1 - \frac {\displaystyle 1 + a}
{\displaystyle 1 + a \exp (b x)},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{b} \log \frac {a + p}{a (1 - p)},
\\
&\displaystyle
{\rm ES}_p (X) = -\left( 1 + \frac {a}{p} \right) \frac {\log a}{b}
+\frac {(a + p) \log (a +p)}{b p}
+ \frac {(1 - p) \log (1 - p)}{b p}
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first scale parameter and \eqn{b > 0}, the second scale parameter.}
\usage{
dperks(x, a=1, b=1, log=FALSE)
pperks(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varperks(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esperks(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{b}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dperks(x)
pperks(x)
varperks(x)
esperks(x)}
 
