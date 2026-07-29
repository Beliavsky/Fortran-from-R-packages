\name{Hlogis}
\alias{dHlogis}
\alias{pHlogis}
\alias{varHlogis}
\alias{esHlogis}
\title{Hosking logistic distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Hosking logistic distribution due to Hosking (1989, 1990) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {(1 - k x)^{1 / k - 1}}{\left[ 1 + (1 - k x)^{1 / k} \right]^2},
\\
&\displaystyle
F (x) = \frac {1}{1 + (1 - k x)^{1 / k}},
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{k} \left[ 1 - \left( \frac {1 - p}{p} \right)^k \right],
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{k} - \frac {1}{kp} B_p (1 - k, 1 + k)
\end{array}}
for \eqn{x < 1/k} if \eqn{k > 0}, \eqn{x > 1/k} if \eqn{k < 0}, \eqn{-\infty < x < \infty} if \eqn{k = 0},
and \eqn{-\infty < k < \infty}, the shape parameter.}
\usage{
dHlogis(x, k=1, log=FALSE)
pHlogis(x, k=1, log.p=FALSE, lower.tail=TRUE)
varHlogis(p, k=1, log.p=FALSE, lower.tail=TRUE)
esHlogis(p, k=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{k}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dHlogis(x)
pHlogis(x)
varHlogis(x)
esHlogis(x)}
 
