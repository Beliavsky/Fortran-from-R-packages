\name{loglogis}
\alias{dloglogis}
\alias{ploglogis}
\alias{varloglogis}
\alias{esloglogis}
\title{Log-logistic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the log-logistic distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {b a^b x^{b - 1}}
{\left( a^b + x^b \right)^2},
\\
&\displaystyle
F (x) =  \frac {x^b}{a^b + x^b},
\\
&\displaystyle
{\rm VaR}_p (X) = a \left( \frac {p}{1 -  p} \right)^{1 / b},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {a}{p} B_p \left( 1 + \frac {1}{b}, 1 - \frac {1}{b} \right)
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{a > 0}, the scale parameter, and \eqn{b > 0}, the shape parameter,
where  \eqn{B_x (a, b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt} denotes the incomplete beta function.}
\usage{
dloglogis(x, a=1, b=1, log=FALSE)
ploglogis(x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
varloglogis(p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
esloglogis(p, a=1, b=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the scale parameter, must be positive, the default is 1}
  \item{b}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dloglogis(x)
ploglogis(x)
varloglogis(x)
esloglogis(x)}
 
