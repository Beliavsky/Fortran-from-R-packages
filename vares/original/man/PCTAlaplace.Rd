\name{PCTAlaplace}
\alias{dPCTAlaplace}
\alias{pPCTAlaplace}
\alias{varPCTAlaplace}
\alias{esPCTAlaplace}
\title{Poiraud-Casanova-Thomas-Agnan Laplace distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Poiraud-Casanova-Thomas-Agnan Laplace distribution due to Poiraud-Casanova and Thomas-Agnan (2000) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \left\{
\begin{array}{ll}
\displaystyle
a \left( 1 - a \right) \exp \left\{ \left( 1 - a \right)
\left( x - \theta \right) \right\}, & \mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
a \left( 1 - a \right) \exp \left\{ a \left( \theta - x \right) \right\}, & \mbox{if $x > \theta$,}
\end{array}
\right.
\\
&\displaystyle
F (x) =
\left\{
\begin{array}{ll}
\displaystyle
a \exp \left\{ \left( 1 - a \right)
\left( x - \theta \right) \right\}, & \mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
1 - \left( 1 - a \right)
\exp \left\{ a \left( \theta - x \right) \right\}, & \mbox{if $x > \theta$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{
\begin{array}{ll}
\displaystyle
\theta + \frac {1}{1 - a} \log \left( \frac {p}{a} \right), & \mbox{if $p \leq a$,}
\\
\\
\displaystyle
\theta - \frac {1}{a} \log \left( \frac {1 - p}{1 - a} \right), & \mbox{if $p > a$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\theta - \frac {\log a}{1 - a} + \frac {\log p - 1}{(1 - a) p}, & \mbox{if $p \leq a$,}
\\
\\
\displaystyle
\theta - \frac {1}{a} + \frac {1}{p} - \frac {a}{(1 - a) p} + \frac {1 - p}{a p}
\log \left( \frac {1 - p}{1 - a} \right), & \mbox{if $p > a$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \theta < \infty}, the location parameter, and
\eqn{a > 0}, the scale parameter.}
\usage{
dPCTAlaplace(x, a=0.5, theta=0, log=FALSE)
pPCTAlaplace(x, a=0.5, theta=0, log.p=FALSE, lower.tail=TRUE)
varPCTAlaplace(p, a=0.5, theta=0, log.p=FALSE, lower.tail=TRUE)
esPCTAlaplace(p, a=0.5, theta=0)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the location parameter, can take any real value, the default is zero}
  \item{a}{the value of the scale parameter, must be in the unit interval, the default is 0.5}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dPCTAlaplace(x)
pPCTAlaplace(x)
varPCTAlaplace(x)
esPCTAlaplace(x)}
 
