\name{Mlaplace}
\alias{dMlaplace}
\alias{pMlaplace}
\alias{varMlaplace}
\alias{esMlaplace}
\title{McGill Laplace distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the McGill Laplace distribution due to McGill (1962) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \left\{
\begin{array}{ll}
\displaystyle
\frac {1}{2 \psi} \exp \left( \frac {x - \theta}{\psi} \right), & \mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
\frac {1}{2 \phi} \exp \left( \frac {\theta - x}{\phi} \right), & \mbox{if $x > \theta$,}
\end{array}
\right.
\\
&\displaystyle
F (x) =  \left\{
\begin{array}{ll}
\displaystyle
\frac {1}{2} \exp \left( \frac {x - \theta}{\psi} \right), & \mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
1 - \frac {1}{2} \exp \left( \frac {\theta - x}{\phi} \right), & \mbox{if $x > \theta$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{
\begin{array}{ll}
\displaystyle
\theta + \psi \log (2 p), & \mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\theta - \phi \log \left( 2 (1 - p) \right), & \mbox{if $p > 1/2$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) = \left\{
\begin{array}{ll}
\displaystyle
\psi + \theta \log (2 p) - \theta p, & \mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\theta + \phi + \frac {\psi - \phi - 2 \theta}{2 p} + \frac {\phi}{p} \log 2 - \phi \log 2
\\
\displaystyle
\quad
+\frac {\phi}{p} \log (1 - p) - \phi \log (1 - p), & \mbox{if $p > 1/2$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \theta < \infty}, the location parameter,
\eqn{\phi > 0}, the first scale parameter, and \eqn{\psi > 0}, the second scale parameter.}
\usage{
dMlaplace(x, theta=0, phi=1, psi=1, log=FALSE)
pMlaplace(x, theta=0, phi=1, psi=1, log.p=FALSE, lower.tail=TRUE)
varMlaplace(p, theta=0, phi=1, psi=1, log.p=FALSE, lower.tail=TRUE)
esMlaplace(p, theta=0, phi=1, psi=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the location parameter, can take any real value, the default is zero}
  \item{phi}{the value of the first scale parameter, must be positive, the default is 1}
  \item{psi}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dMlaplace(x)
pMlaplace(x)
varMlaplace(x)
esMlaplace(x)}
 
