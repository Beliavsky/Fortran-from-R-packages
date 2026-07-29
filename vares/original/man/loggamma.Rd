\name{loggamma}
\alias{dloggamma}
\alias{ploggamma}
\alias{varloggamma}
\alias{esloggamma}
\title{Log gamma distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the log gamma distribution due to Consul and Jain (1971) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {a^r x^{a - 1} (-\log x)^{r - 1}}{\Gamma (r)},
\\
&\displaystyle
F (x) = Q (r, -a \log x),
\\
&\displaystyle
{\rm VaR}_p (X) = \exp \left[ -\frac {1}{a} Q^{-1} (r, p) \right],
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p \exp
\left[ -\frac {1}{a} Q^{-1} (r, v) \right] dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, and \eqn{r > 0}, the second shape parameter.}
\usage{
dloggamma(x, a=1, r=1, log=FALSE)
ploggamma(x, a=1, r=1, log.p=FALSE, lower.tail=TRUE)
varloggamma(p, a=1, r=1, log.p=FALSE, lower.tail=TRUE)
esloggamma(p, a=1, r=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be positive, the default is 1}
  \item{r}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dloggamma(x)
ploggamma(x)
varloggamma(x)
esloggamma(x)}
 
