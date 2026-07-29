\name{clg}
\alias{dclg}
\alias{pclg}
\alias{varclg}
\alias{esclg}
\title{Compound Laplace gamma distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the compound Laplace gamma distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {a b}{2} \left\{ 1 + b \left | x - \theta \right | \right\}^{-\left( a + 1 \right)},
\\
&\displaystyle
F (x) =
\left\{
\begin{array}{ll}
\displaystyle
\frac {1}{2} \left\{ 1 + b \left | x - \theta \right | \right\}^{-a}, &
\mbox{if $x \leq \theta$,}
\\
\\
\displaystyle
1 - \frac {1}{2} \left\{ 1 + b \left | x - \theta \right | \right\}^{-a}, &
\mbox{if $x > \theta$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\theta - \frac {1}{b} - \frac {(2 p)^{-1/a}}{b}, & \mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\theta - \frac {1}{b} + \frac {(2 (1 - p))^{-1/a}}{b}, &
\mbox{if $p > 1/2$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) =  \left\{
\begin{array}{ll}
\displaystyle
\theta - \frac {1}{b} - \frac {(2 p)^{-1/a}}{b (1 - 1/a)}, &
\mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\theta - \frac {1}{b} - \frac {\left[ 2 (1 - p) \right]^{1 - 1/a}}{2 p b (1 - 1/a)}, &
\mbox{if $p > 1/2$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \theta < \infty}, the location parameter,
\eqn{b > 0}, the scale parameter, and \eqn{a > 0}, the shape parameter.}
\usage{
dclg(x, a=1, b=1, theta=0, log=FALSE)
pclg(x, a=1, b=1, theta=0, log.p=FALSE, lower.tail=TRUE)
varclg(p, a=1, b=1, theta=0, log.p=FALSE, lower.tail=TRUE)
esclg(p, a=1, b=1, theta=0)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the location parameter, can take any real value, the default is zero}
  \item{b}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dclg(x)
pclg(x)
varclg(x)
esclg(x)}
 
