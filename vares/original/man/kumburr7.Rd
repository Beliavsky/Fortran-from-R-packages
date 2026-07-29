\name{kumburr7}
\alias{dkumburr7}
\alias{pkumburr7}
\alias{varkumburr7}
\alias{eskumburr7}
\title{Kumaraswamy Burr XII distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Kumaraswamy Burr XII distribution due to Parana\'{i}ba et al. (2013) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {a b k c x^{c - 1}}{\left( 1 + x^c \right)^{k + 1}} \left[ 1 - \left( 1 + x^c \right)^{-k} \right]^{a - 1}
\left\{ 1 - \left[ 1 - \left( 1 + x^c \right)^{-k} \right]^a \right\}^{b - 1},
\\
&\displaystyle
F (x) = 1 - \left\{ 1 - \left[ 1 - \left( 1 + x^c \right)^{-k} \right]^a \right\}^b,
\\
&\displaystyle
{\rm VaR}_p (X) = \left[ \left\{ 1 - \left[ 1 - (1 - p)^{1 / b} \right]^{1 / a} \right\}^{-1 / k} - 1 \right]^{1/c},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p \left[ \left\{ 1 - \left[ 1 - (1 - v)^{1 / b} \right]^{1 / a} \right\}^{-1 / k} - 1 \right]^{1/c} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the first shape parameter, \eqn{b > 0}, the second shape parameter,
\eqn{c > 0}, the third shape parameter, and \eqn{k > 0}, the fourth shape parameter.}
\usage{
dkumburr7(x, a=1, b=1, k=1, c=1, log=FALSE)
pkumburr7(x, a=1, b=1, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
varkumburr7(p, a=1, b=1, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
eskumburr7(p, a=1, b=1, k=1, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{c}{the value of the third shape parameter, must be positive, the default is 1}
  \item{k}{the value of the fourth shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dkumburr7(x)
pkumburr7(x)
varkumburr7(x)
eskumburr7(x)}
 
