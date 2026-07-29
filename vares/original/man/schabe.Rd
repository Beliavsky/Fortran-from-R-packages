\name{schabe}
\alias{dschabe}
\alias{pschabe}
\alias{varschabe}
\alias{esschabe}
\title{Schabe distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Schabe distribution due to Schabe (1994) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {\displaystyle 2 \gamma + (1 - \gamma) x / \theta}{\displaystyle \theta (\gamma + x/\theta)^2},
\\
&\displaystyle
F(x) = \frac {\displaystyle  (1 + \gamma) x}{\displaystyle  x + \gamma \theta},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {p \gamma \theta}{1 + \gamma - p},
\\
&\displaystyle
{\rm ES}_p (X) = -\theta \gamma - \frac {\theta \gamma (1 + \gamma)}{p}
\log \frac {1 + \gamma - p}{1 + \gamma}
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{0  < \gamma < 1}, the first scale parameter, and \eqn{\theta > 0}, the second scale parameter.}
\usage{
dschabe(x, gamma=0.5, theta=1, log=FALSE)
pschabe(x, gamma=0.5, theta=1, log.p=FALSE, lower.tail=TRUE)
varschabe(p, gamma=0.5, theta=1, log.p=FALSE, lower.tail=TRUE)
esschabe(p, gamma=0.5, theta=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{gamma}{the value of the first scale parameter, must be in the unit interval, the default is 0.5}
  \item{theta}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dschabe(x)
pschabe(x)
varschabe(x)
esschabe(x)}
 
