\name{halfT}
\alias{dhalfT}
\alias{phalfT}
\alias{varhalfT}
\alias{eshalfT}
\title{Half Student's t distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the half Student's \eqn{t} distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {2 \Gamma \left( \frac {n + 1}{2} \right)}{\sqrt{n \pi} \Gamma \left( \frac {n}{2} \right)}
\left( 1 + \frac {x^2}{n} \right)^{-\frac {n + 1}{2}},
\\
&\displaystyle
F (x) =  I_{\frac {x^2}{x^2 + n}} \left( \frac {1}{2}, \frac {n}{2} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \sqrt{\frac {n I_p^{-1} \left( \frac {1}{2}, \frac {n}{2} \right)}
{1 - I_p^{-1} \left( \frac {1}{2}, \frac {n}{2} \right)}},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\sqrt{n}}{p} \int_0^p
\sqrt{\frac {I_v^{-1} \left( \frac {1}{2}, \frac {n}{2} \right)}
{1 - I_v^{-1} \left( \frac {1}{2}, \frac {n}{2} \right)}} dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, and \eqn{n > 0}, the degree of freedom parameter.}
\usage{
dhalfT(x, n=1, log=FALSE)
phalfT(x, n=1, log.p=FALSE, lower.tail=TRUE)
varhalfT(p, n=1, log.p=FALSE, lower.tail=TRUE)
eshalfT(p, n=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{n}{the value of the degree of freedom parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dhalfT(x)
phalfT(x)
varhalfT(x)
eshalfT(x)}
 
