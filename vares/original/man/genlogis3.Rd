
\name{genlogis3}
\alias{dgenlogis3}
\alias{pgenlogis3}
\alias{vargenlogis3}
\alias{esgenlogis3}
\title{Generalized logistic III distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the generalized logistic III distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\sigma B (\alpha, \alpha)} \exp \left( \alpha \frac {x - \mu}{\sigma} \right)
\left\{ 1 + \exp \left( \frac {x - \mu}{\sigma} \right) \right\}^{-2 \alpha},
\\
&\displaystyle
F (x) = I_{\frac {1}{1 + \exp \left( -\frac {x - \mu}{\sigma} \right)}} (\alpha, \alpha),
\\
&\displaystyle
{\rm VaR}_p (X) = \mu - \sigma \log \frac {1 - I_p^{-1} (\alpha, \alpha)}{I_p^{-1} (\alpha, \alpha)},
\\
&\displaystyle
{\rm ES}_p (X) = \mu - \frac {\sigma}{p} \int_0^p \log \frac {1 - I_v^{-1} (\alpha, \alpha)}{I_v^{-1} (\alpha, \alpha)} dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, \eqn{\sigma > 0}, the scale parameter,
and \eqn{\alpha > 0}, the shape parameter.}
\usage{
dgenlogis3(x, alpha=1, mu=0, sigma=1, log=FALSE)
pgenlogis3(x, alpha=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
vargenlogis3(p, alpha=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
esgenlogis3(p, alpha=1, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{alpha}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dgenlogis3(x)
pgenlogis3(x)
vargenlogis3(x)
esgenlogis3(x)}
 
