\name{pareto}
\alias{dpareto}
\alias{ppareto}
\alias{varpareto}
\alias{espareto}
\title{Pareto distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Pareto distribution due to Pareto (1964) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = c K^c x^{-c - 1},
\\
&\displaystyle
F (x) = 1 - \left( \frac {K}{x} \right)^c,
\\
&\displaystyle
{\rm VaR}_p (X) = K (1 - p)^{-1 / c},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {K c}{p (1 - c)} (1 - p)^{1 - 1 / c} - \frac {K c}{p (1 - c)}
\end{array}}
for \eqn{x \geq K}, \eqn{0 < p < 1}, \eqn{K > 0}, the scale parameter, and \eqn{c > 0}, the shape parameter.}
\usage{
dpareto(x, K=1, c=1, log=FALSE)
ppareto(x, K=1, c=1, log.p=FALSE, lower.tail=TRUE)
varpareto(p, K=1, c=1, log.p=FALSE, lower.tail=TRUE)
espareto(p, K=1, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{K}{the value of the scale parameter, must be positive, the default is 1}
  \item{c}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dpareto(x)
ppareto(x)
varpareto(x)
espareto(x)}
 
