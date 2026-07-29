\name{MRbeta}
\alias{dMRbeta}
\alias{pMRbeta}
\alias{varMRbeta}
\alias{esMRbeta}
\title{McDonald-Richards beta distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the McDonald-Richards beta distribution due to McDonald and Richards (1987a, 1987b) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {x^{ar - 1} \left( bq^r - x^r \right)^{b - 1}}
{\left( b q^r \right)^{a + b - 1} B (a, b)},
\\
&\displaystyle
F (x) = I_{\frac {x^r}{b q^r}} (a, b),
\\
&\displaystyle
{\rm VaR}_p (X) = b^{1/r} q \left[ I_p^{-1} (a, b) \right]^{1/r},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {b^{1/r} q}{p} \int_0^p \left[ I_v^{-1} (a, b) \right]^{1/r} dv
\end{array}}
for \eqn{0 \leq x \leq b^{1 / r} q}, \eqn{0 < p < 1}, \eqn{q > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter,
\eqn{b > 0}, the second shape parameter, and \eqn{r > 0}, the third shape parameter.}
\usage{
dMRbeta(x, a=1, b=1, r=1, q=1, log=FALSE)
pMRbeta(x, a=1, b=1, r=1, q=1, log.p=FALSE, lower.tail=TRUE)
varMRbeta(p, a=1, b=1, r=1, q=1, log.p=FALSE, lower.tail=TRUE)
esMRbeta(p, a=1, b=1, r=1, q=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{q}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{r}{the value of the third shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dMRbeta(x)
pMRbeta(x)
varMRbeta(x)
esMRbeta(x)}
 
