\name{Gamma}
\alias{dGamma}
\alias{pGamma}
\alias{varGamma}
\alias{esGamma}
\title{Gamma distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the gamma distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {b^a x^{a - 1} \exp (-b x)}{\Gamma (a)},
\\
&\displaystyle
F (x) = \frac {\gamma (a, b x)}{\Gamma (a)},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{b} Q^{-1} (a, 1 - p),
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{b p} \int_0^p Q^{-1} (a, 1 - v) dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{b > 0}, the scale parameter, and \eqn{a > 0}, the shape parameter,
where \eqn{\gamma (a, x) = \int_0^x t^{a - 1} \exp \left( -t \right) dt} denotes the incomplete gamma function,
\eqn{Q (a, x) = \int_x^\infty t^{a - 1} \exp \left( -t \right) dt / \Gamma (a)}
denotes the regularized complementary incomplete gamma function,
\eqn{\Gamma (a) = \int_0^\infty t^{a - 1} \exp \left( -t \right) dt} denotes the gamma function,
and \eqn{Q^{-1} (a, x)} denotes the inverse of \eqn{Q (a, x)}.}
\usage{
dGamma(x, a=1, b=1, log=FALSE)
pGamma(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varGamma(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esGamma(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{b}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dGamma(x)
pGamma(x)
varGamma(x)
esGamma(x)}
 

