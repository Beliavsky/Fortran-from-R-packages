\name{frechet}
\alias{dfrechet}
\alias{pfrechet}
\alias{varfrechet}
\alias{esfrechet}
\title{Frechet distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Fr\'{e}chet distribution due to Fr\'{e}chet (1927) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {\alpha \sigma^\alpha}{x^{\alpha + 1}} \exp \left\{ -\left( \frac {\sigma}{x} \right)^{\alpha} \right\},
\\
&\displaystyle
F (x) = \exp \left\{ -\left( \frac {\sigma}{x} \right)^{\alpha} \right\},
\\
&\displaystyle
{\rm VaR}_p (X) = \sigma \left[ -\log p \right]^{-1 / \alpha},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sigma}{p} \Gamma \left( 1 - 1 / \alpha, -\log p \right)
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{\alpha > 0}, the shape parameter, and \eqn{\sigma > 0}, the scale parameter,
where \eqn{\Gamma (a, x) = \int_x^\infty t^{a - 1} \exp \left( -t \right) dt}
denotes the complementary incomplete gamma function.}
\usage{
dfrechet(x, alpha=1, sigma=1, log=FALSE)
pfrechet(x, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
varfrechet(p, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
esfrechet(p, alpha=1, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{alpha}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dfrechet(x)
pfrechet(x)
varfrechet(x)
esfrechet(x)}
 
