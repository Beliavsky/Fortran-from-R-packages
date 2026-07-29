\name{kumpareto}
\alias{dkumpareto}
\alias{pkumpareto}
\alias{varkumpareto}
\alias{eskumpareto}
\title{Kumaraswamy Pareto distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Kumaraswamy Pareto distribution due to Pereira et al. (2013) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = a b c K^c x^{-c - 1}
\left[ 1 - \left( \frac {K}{x} \right)^c \right]^{a - 1}
\left\{ 1 - \left[ 1 - \left( \frac {K}{x} \right)^c \right]^a \right\}^{b - 1},
\\
&\displaystyle
F (x) = 1 - \left\{ 1 - \left[ 1 - \left( \frac {K}{x} \right)^c \right]^a \right\}^b,
\\
&\displaystyle
{\rm VaR}_p (X) = K \left\{ 1 - \left[ 1 - (1 - p)^{1 / b} \right]^{1 / a} \right\}^{-1 / c},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {K}{p} \int_0^p \left\{ 1 - \left[ 1 - (1 - v)^{1 / b} \right]^{1 / a} \right\}^{-1 / c} dv
\end{array}}
for \eqn{x \geq K}, \eqn{0 < p < 1}, \eqn{K > 0}, the scale parameter, \eqn{c > 0}, the first shape parameter,
\eqn{a > 0}, the second shape parameter, and \eqn{b > 0}, the third shape parameter.}
\usage{
dkumpareto(x, K=1, a=1, b=1, c=1, log=FALSE)
pkumpareto(x, K=1, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
varkumpareto(p, K=1, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
eskumpareto(p, K=1, a=1, b=1, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{K}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{c}{the value of the third shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dkumpareto(x)
pkumpareto(x)
varkumpareto(x)
eskumpareto(x)}
 
