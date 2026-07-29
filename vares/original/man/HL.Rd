\name{HL}
\alias{varHL}
\alias{esHL}
\title{Hankin-Lee distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Hankin-Lee distribution due to Hankin and Lee (2006) given by
\deqn{\begin{array}{ll}
&\displaystyle
{\rm VaR}_p (X) = \frac {c p^a}{(1 - p)^b},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {c}{p} B_p (a + 1, 1 - b)
\end{array}}
for \eqn{0 < p < 1}, \eqn{c > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter,
and \eqn{b > 0}, the second shape parameter.}
\usage{
varHL(p, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
esHL(p, a=1, b=1, c=1)
}
\arguments{
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{c}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
varHL(x)
esHL(x)}
 
