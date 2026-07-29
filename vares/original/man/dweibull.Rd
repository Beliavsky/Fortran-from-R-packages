\name{dweibull}
\alias{ddweibull}
\alias{pdweibull}
\alias{vardweibull}
\alias{esdweibull}
\title{Double Weibull distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the double Weibull distribution due to Balakrishnan and Kocherlakota (1985) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \frac {c}{\displaystyle 2 \sigma}
\left | \frac {x - \mu}{\sigma} \right |^{c - 1}
\exp \left\{ -\left | \frac {x - \mu}{\sigma} \right |^c \right\},
\\
&\displaystyle
F (x) = \left\{
\begin{array}{ll}
\displaystyle
\frac {1}{2} \exp \left\{ -\left( \frac {\mu - x}{\sigma} \right)^c \right\}, & \mbox{if $x \leq \mu$,}
\\
\\
\displaystyle
1 - \frac {1}{2} \exp \left\{ -\left( \frac {x - \mu}{\sigma} \right)^c \right\}, & \mbox{if $x > \mu$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) =  \left\{
\begin{array}{ll}
\displaystyle
\mu - \sigma \left[ -\log \left( 2 p \right) \right]^{1 / c}, & \mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\mu + \sigma \left[ -\log \left( 2 (1 - p) \right) \right]^{1 / c}, & \mbox{if $p > 1/2$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) =  \left\{
\begin{array}{ll}
\displaystyle
\mu - \frac {\sigma}{p} \int_0^p \left[ -\log 2 - \log v \right]^{1 / c} dv, & \mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\mu - \frac {\sigma}{p} \int_0^{1/2} \left[ -\log 2 - \log v \right]^{1 / c} dv
\\
\quad
\displaystyle
+ \frac {\sigma}{p}
\int_{1/2}^p \left[ -\log 2 - \log (1 - v) \right]^{1 / c} dv, & \mbox{if $p > 1/2$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter,
\eqn{\sigma > 0}, the scale parameter, and \eqn{c > 0}, the shape parameter.}
\usage{
ddweibull(x, c=1, mu=0, sigma=1, log=FALSE)
pdweibull(x, c=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
vardweibull(p, c=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
esdweibull(p, c=1, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{c}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
ddweibull(x)
pdweibull(x)
vardweibull(x)
esdweibull(x)}
 
