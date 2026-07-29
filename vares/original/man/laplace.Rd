\name{laplace}
\alias{dlaplace}
\alias{plaplace}
\alias{varlaplace}
\alias{eslaplace}
\title{Laplace distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Laplace distribution due to due to Laplace (1774) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{2 \sigma} \exp \left( -\frac {\mid x - \mu \mid}{\sigma} \right),
\\
&\displaystyle
F (x) =  \left\{
\begin{array}{ll}
\displaystyle
\frac {1}{2} \exp \left( \frac {x - \mu}{\sigma} \right), & \mbox{if $x < \mu$,}
\\
\\
\displaystyle
1 - \frac {1}{2} \exp \left( -\frac {x - \mu}{\sigma} \right), & \mbox{if $x \geq \mu$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) =  \left\{
\begin{array}{ll}
\displaystyle
\mu + \sigma \log (2 p), & \mbox{if $p < 1/2$,}
\\
\\
\displaystyle
\mu - \sigma \log \left[ 2 (1 - p) \right], & \mbox{if $p \geq 1/2$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\mu + \sigma \left[ \log (2 p) - 1 \right], & \mbox{if $p < 1/2$,}
\\
\\
\displaystyle
\mu + \sigma - \frac {\sigma}{p} + \sigma \frac {1 - p}{p} \log (1 - p) +\sigma \frac {1 - p}{p} \log 2, & \mbox{if $p \geq 1/2$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, and
\eqn{\sigma > 0}, the scale parameter.}
\usage{
dlaplace(x, mu=0, sigma=1, log=FALSE)
plaplace(x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
varlaplace(p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
eslaplace(p, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dlaplace(x)
plaplace(x)
varlaplace(x)
eslaplace(x)}
 
