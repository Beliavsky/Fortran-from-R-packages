\name{HBlaplace}
\alias{dHBlaplace}
\alias{pHBlaplace}
\alias{varHBlaplace}
\alias{esHBlaplace}
\title{Holla-Bhattacharya Laplace distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the Holla-Bhattacharya Laplace distribution due to Holla and Bhattacharya (1968) given by
\deqn{\begin{array}{ll}
&
f (x) = \left\{ \begin{array}{ll}
\displaystyle
a \phi \exp \left\{ \phi \left( x - \theta \right) \right\}, & \mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
\left( 1 - a \right) \phi \exp \left\{ \phi \left( \theta - x \right) \right\}, & \mbox{if $x > \theta$,}
\end{array}
\right.
\\
&
F (x) = \left\{
\begin{array}{ll}
\displaystyle
a \exp \left( \phi x - \theta \phi \right), & \mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
1 - (1 - a) \exp \left( \theta \phi - \phi x \right), & \mbox{if $x > \theta$,}
\end{array}
\right.
\\
&
{\rm VaR}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\theta + \frac {1}{\phi} \log \left( \frac {p}{a} \right), &
\mbox{if $p \leq a$,}
\\
\\
\displaystyle
\theta - \frac {1}{\phi} \log \left(  \frac {1 - p}{1 - a} \right), &
\mbox{if $p > a$,}
\end{array}
\right.
\\
&
{\rm ES}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\theta - \frac {1}{\phi} + \frac {1}{\phi} \log \frac {p}{a}, & \mbox{if $p \leq a$,}
\\
\\
\displaystyle
\frac {1}{p} \left[ \theta (1 + p - a) + \frac {p - 2a - (1 - a) \log a}{\phi} +
\frac {1 - p}{\phi} \log \frac {1 - p}{1 - a} \right], & \mbox{if $p > a$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \theta < \infty}, the location parameter,
\eqn{0 < a < 1}, the first scale parameter, and \eqn{\phi > 0}, the second scale parameter.}
\usage{
dHBlaplace(x, a=0.5, theta=0, phi=1, log=FALSE)
pHBlaplace(x, a=0.5, theta=0, phi=1, log.p=FALSE, lower.tail=TRUE)
varHBlaplace(p, a=0.5, theta=0, phi=1, log.p=FALSE, lower.tail=TRUE)
esHBlaplace(p, a=0.5, theta=0, phi=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the location parameter, can take any real value, the default is zero}
  \item{a}{the value of the first scale parameter, must be in the unit interval, the default is 0.5}
  \item{phi}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dHBlaplace(x)
pHBlaplace(x)
varHBlaplace(x)
esHBlaplace(x)}
 
