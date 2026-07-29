\name{explogis}
\alias{dexplogis}
\alias{pexplogis}
\alias{varexplogis}
\alias{esexplogis}
\title{Exponentiated logistic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponentiated logistic distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =   (a/b) \exp (-x/b)
\left[ 1 + \exp (-x/b) \right]^{-a - 1},
\\
&\displaystyle
F (x) =   \left[ 1 + \exp (-x/b) \right]^{-a},
\\
&\displaystyle
{\rm VaR}_p (X) = -b \log \left[ p^{-1 / a} - 1 \right],
\\
&\displaystyle
{\rm ES}_p (X) = -\frac {b}{p} \int_0^p \log \left[ v^{-1 / a} - 1 \right] dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter, and \eqn{b > 0}, the scale parameter.}
\usage{
dexplogis(x, a=1, b=1, log=FALSE)
pexplogis(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varexplogis(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esexplogis(p, a=1, b=1)
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
dexplogis(x)
pexplogis(x)
varexplogis(x)
esexplogis(x)}
 
