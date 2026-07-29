\name{burr7}
\alias{dburr7}
\alias{pburr7}
\alias{varburr7}
\alias{esburr7}
\title{Burr XII distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Burr XII distribution due to Burr (1942) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {k c x^{c - 1}}{\left( 1 + x^c \right)^{k + 1}},
\\
&\displaystyle
F (x) = 1 - \left( 1 + x^c \right)^{-k},
\\
&\displaystyle
{\rm VaR}_p (X) = \left[ (1 - p)^{-1 / k} - 1 \right]^{1/c},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {1}{p} \int_0^p \left[ (1 - v)^{-1 / k} - 1 \right]^{1/c} dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{c > 0}, the first shape parameter, and \eqn{k > 0}, the second shape parameter.}
\usage{
dburr7(x, k=1, c=1, log=FALSE)
pburr7(x, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
varburr7(p, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
esburr7(p, k=1, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{k}{the value of the first shape parameter, must be positive, the default is 1}
  \item{c}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dburr7(x)
pburr7(x)
varburr7(x)
esburr7(x)}
 
