\name{quad}
\alias{dquad}
\alias{pquad}
\alias{varquad}
\alias{esquad}
\title{Quadratic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the quadratic distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \alpha (x - \beta)^2,
\\
&\displaystyle
F(x) = \frac {\alpha}{3} \left[ (x - \beta)^3 + (\beta - a)^3 \right],
\\
&\displaystyle
{\rm VaR}_p (X) = \beta + \left[ \frac {3 p}{\alpha} - (\beta - a)^3 \right]^{1/3},
\\
&\displaystyle
{\rm ES}_p (X) = \beta + \frac {\alpha}{4 p} \left\{ \left[ \frac {3 p}{\alpha} - (\beta - a)^3 \right]^{4/3} - (\beta - a)^4 \right\}
\end{array}}
for \eqn{a \leq x \leq b}, \eqn{0 < p < 1}, \eqn{-\infty < a < \infty} , the first location parameter, and \eqn{-\infty < a < b < \infty}, the second location parameter, where \eqn{\alpha = \frac {12}{(b - a)^3}} and \eqn{\beta = \frac {a + b}{2}}.}
\usage{
dquad(x, a=0, b=1, log=FALSE)
pquad(x, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
varquad(p, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
esquad(p, a=0, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first location parameter, can take any real value, the default is zero}
  \item{b}{the value of the second location parameter, can take any real value but must be greater than a, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dquad(x)
pquad(x)
varquad(x)
esquad(x)}
 
