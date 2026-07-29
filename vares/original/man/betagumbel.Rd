\name{betagumbel}
\alias{dbetagumbel}
\alias{pbetagumbel}
\alias{varbetagumbel}
\alias{esbetagumbel}
\title{Beta Gumbel distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the beta Gumbel distribution due to Nadarajah and Kotz (2004) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\sigma B (a, b)}
\exp \left( \frac {\mu - x}{\sigma} \right)
\exp \left[ -a \exp \frac {\mu - x}{\sigma} \right]
\left\{ 1 - \exp \left[ -\exp \frac {\mu - x}{\sigma} \right] \right\}^{b - 1},
\\
&\displaystyle
F (x) = I_{\exp \left[ -\exp \frac {\mu - x}{\sigma} \right]} (a, b),
\\
&\displaystyle
{\rm VaR}_p (X) = \mu - \sigma \log \left[ -\log I_p^{-1} (a, b) \right],
\\
&\displaystyle
{\rm ES}_p (X) = \mu - \frac {\sigma}{p} \int_0^p \log \left[ -\log I_v^{-1} (a, b) \right] dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, \eqn{\sigma > 0}, the scale parameter,
\eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dbetagumbel(x, a=1, b=1, mu=0, sigma=1, log=FALSE)
pbetagumbel(x, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
varbetagumbel(p, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
esbetagumbel(p, a=1, b=1, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbetagumbel(x)
pbetagumbel(x)
varbetagumbel(x)
esbetagumbel(x)}
 
