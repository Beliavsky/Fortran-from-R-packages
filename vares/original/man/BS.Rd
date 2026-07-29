\name{BS}
\alias{dBS}
\alias{pBS}
\alias{varBS}
\alias{esBS}
\title{Birnbaum-Saunders distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Birnbaum-Saunders distribution due to Birnbaum and Saunders (1969a, 1969b) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {x^{1/2} + x^{-1/2}}{2 \gamma x}
\phi \left( \frac {x^{1/2} - x^{-1/2}}{\gamma} \right),
\\
&\displaystyle
F (x) = \Phi \left( \frac {x^{1/2} - x^{-1/2}}{\gamma} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {1}{4} \left\{ \gamma \Phi^{-1} (p) + \sqrt{4 + \gamma^2 \left[ \Phi^{-1} (p) \right]^2} \right\}^2,
\\
&\displaystyle
{\rm ES}_p (X) =
\frac {1}{4 p} \int_0^p \left\{ \gamma \Phi^{-1} (v) + \sqrt{4 + \gamma^2 \left[ \Phi^{-1} (v) \right]^2} \right\}^2 dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, and \eqn{\gamma > 0}, the scale parameter.}
\usage{
dBS(x, gamma=1, log=FALSE)
pBS(x, gamma=1, log.p=FALSE, lower.tail=TRUE)
varBS(p, gamma=1, log.p=FALSE, lower.tail=TRUE)
esBS(p, gamma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{gamma}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dBS(x)
pBS(x)
varBS(x)
esBS(x)}
 
