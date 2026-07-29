\name{betalognorm}
\alias{dbetalognorm}
\alias{pbetalognorm}
\alias{varbetalognorm}
\alias{esbetalognorm}
\title{Beta lognormal distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the beta lognormal distribution due to Castellares et al. (2013) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\sigma x B (a, b)}
\phi \left( \frac {\log x - \mu}{\sigma} \right)
\Phi^{a - 1} \left( \frac {\log x - \mu}{\sigma} \right)
\Phi^{b - 1} \left( \frac {\mu - \log x}{\sigma} \right),
\\
&\displaystyle
F (x) = I_{\Phi \left( \frac {\log x - \mu}{\sigma} \right)} (a, b),
\\
&\displaystyle
{\rm VaR}_p (X) =  \exp \left[ \mu + \sigma \Phi^{-1} \left( I_p^{-1} (a, b) \right) \right],
\\
&\displaystyle
{\rm ES}_p (X) = \frac {\exp (\mu)}{p} \int_0^p \exp \left[ \sigma \Phi^{-1} \left( I_v^{-1} (a, b) \right) \right] dv
\end{array}}
for \eqn{x > 0}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, \eqn{\sigma > 0}, the scale parameter,
\eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter,
where \eqn{\phi (\cdot)} denotes the pdf of a standard normal random variable,
and \eqn{\Phi (\cdot)} denotes the cdf of a standard normal random variable.}
\usage{
dbetalognorm(x, a=1, b=1, mu=0, sigma=1, log=FALSE)
pbetalognorm(x, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
varbetalognorm(p, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
esbetalognorm(p, a=1, b=1, mu=0, sigma=1)
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
dbetalognorm(x)
pbetalognorm(x)
varbetalognorm(x)
esbetalognorm(x)}
 

