\name{genpareto}
\alias{dgenpareto}
\alias{pgenpareto}
\alias{vargenpareto}
\alias{esgenpareto}
\title{Generalized Pareto distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the generalized Pareto distribution due to Pickands (1975) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {1}{k} \left( 1 - \frac {c x}{k} \right)^{1 / c - 1},
\\
&\displaystyle
F (x) =  1 - \left( 1 - \frac {c x}{k} \right)^{1 / c},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {k}{c} \left[ 1 - (1 - p)^c \right],
\\
&\displaystyle
{\rm ES}_p (X) = \frac {k}{c} + \frac {k (1 - p)^{c + 1}}{p c (c + 1)} - \frac {k}{p c (c + 1)}
\end{array}}
for \eqn{x < k/c} if \eqn{c > 0}, \eqn{x > k/c} if \eqn{c < 0}, \eqn{x > 0} if \eqn{c = 0}, \eqn{0 < p < 1},
\eqn{k > 0}, the scale parameter and \eqn{-\infty < c < \infty}, the shape parameter.}
\usage{
dgenpareto(x, k=1, c=1, log=FALSE)
pgenpareto(x, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
vargenpareto(p, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
esgenpareto(p, k=1, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{k}{the value of the scale parameter, must be positive, the default is 1}
  \item{c}{the value of the shape parameter, can take any real value, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgenpareto(x)
pgenpareto(x)
vargenpareto(x)
esgenpareto(x)}
 
