\name{kumhalfnorm}
\alias{dkumhalfnorm}
\alias{pkumhalfnorm}
\alias{varkumhalfnorm}
\alias{eskumhalfnorm}
\title{Kumaraswamy half normal distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Kumaraswamy half normal distribution due to Cordeiro et al. (2012c) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {2a b}{\sigma} \phi \left( \frac {x}{\sigma} \right) \left[ 2 \Phi \left( \frac {x}{\sigma} \right) - 1 \right]^{a - 1}
\left\{ 1 - \left[ 2 \Phi \left( \frac {x}{\sigma} \right) - 1 \right]^a \right\}^{b - 1},
\\
&\displaystyle
F (x) =  1 - \left\{ 1 - \left[ 2 \Phi \left( \frac {x}{\sigma} \right) - 1 \right]^a \right\}^b,
\\
&\displaystyle
{\rm VaR}_p (X) = \sigma \Phi^{-1} \left( \frac {1}{2} + \frac {1}{2}
\left[ 1 - (1 - p)^{1 / b} \right]^{1 / a} \right),
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sigma}{p} \int_0^p \Phi^{-1} \left( \frac {1}{2} + \frac {1}{2}
\left[ 1 - (1 - v)^{1 / b} \right]^{1 / a} \right) dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1},
\eqn{\sigma > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dkumhalfnorm(x, sigma=1, a=1, b=1, log=FALSE)
pkumhalfnorm(x, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varkumhalfnorm(p, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
eskumhalfnorm(p, sigma=1, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
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
dkumhalfnorm(x)
pkumhalfnorm(x)
varkumhalfnorm(x)
eskumhalfnorm(x)}
 
