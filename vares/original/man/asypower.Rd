\name{asypower}
\alias{dasypower}
\alias{pasypower}
\alias{varasypower}
\alias{esasypower}
\title{Asymmetric power distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the asymmetric power distribution due to Komunjer (2007) given by
\deqn{\begin{array}{ll}
&\displaystyle
f(x) = \left\{ \begin{array}{ll}
\displaystyle
\frac {\displaystyle \delta^{1 / \lambda}}{\displaystyle \Gamma (1 + 1 / \lambda)}
\exp \left[ -\frac {\delta}{a^\lambda} |x|^\lambda \right], & \mbox{if $x \leq 0$},
\\
\\
\displaystyle
\frac {\displaystyle \delta^{1 / \lambda}}{\displaystyle \Gamma (1 + 1 / \lambda)}
\exp \left[ -\frac {\delta}{(1 - a)^\lambda} |x|^\lambda \right], & \mbox{if $x > 0$,}
\end{array}
\right.
\\
&\displaystyle
F (x) = \left\{ \begin{array}{ll}
\displaystyle
a - a {\cal I} \left( \frac {\delta}{a^\lambda} \sqrt{\lambda} |x|^\lambda, 1 / \lambda \right), & \mbox{if $x \leq 0$,}
\\
\\
\displaystyle
a - (1 - a) {\cal I} \left( \frac {\delta}{(1 - a)^\lambda}
\sqrt{\lambda} |x|^\lambda, 1 / \lambda \right), & \mbox{if $x > 0$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
-\left[ \frac {a^\lambda}{\delta \sqrt{\lambda}} \right]^{1 / \lambda}
\left[ {\cal I}^{-1} \left( 1 - \frac {p}{a}, \frac {1}{\lambda} \right) \right]^{1 / \lambda}, & \mbox{if $p \leq a$,}
\\
\\
\displaystyle
-\left[ \frac {(1 - a)^\lambda}{\delta \sqrt{\lambda}} \right]^{1 / \lambda}
\left[ {\cal I}^{-1} \left( 1 - \frac {1 - p}{1 - a}, \frac {1}{\lambda} \right) \right]^{1 / \lambda}, & \mbox{if $p > a$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
-\frac {1}{p} \left[ \frac {a^\lambda}{\delta \sqrt{\lambda}} \right]^{1 / \lambda}
\int_0^p \left[ {\cal I}^{-1} \left( 1 - \frac {v}{a}, \frac {1}{\lambda} \right) \right]^{1 / \lambda} dv, &
\mbox{if $p \leq a$,}
\\
\\
\displaystyle
-\frac {1}{p} \left[ \frac {a^\lambda}{\delta \sqrt{\lambda}} \right]^{1 / \lambda}
\int_0^a \left[ {\cal I}^{-1} \left( 1 - \frac {v}{a}, \frac {1}{\lambda} \right) \right]^{1 / \lambda} dv
\\
\quad
\displaystyle
-\frac {1}{p} \left[ \frac {(1 - a)^\lambda}{\delta \sqrt{\lambda}} \right]^{1 / \lambda}
\int_a^p \left[ {\cal I}^{-1} \left( 1 -
\frac {1 - v}{1 - a}, \frac {1}{\lambda} \right) \right]^{1 / \lambda} dv, & \mbox{if $p > a$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{0 < a < 1}, the first scale parameter, \eqn{\delta > 0}, the second scale parameter,
and \eqn{\lambda > 0}, the shape parameter,
where \eqn{{\cal I} (x, \gamma) = \frac {1}{\Gamma (\gamma)} \int_0^{x \sqrt{\gamma}} t^{\gamma - 1} \exp (-t) dt}.}
\usage{
dasypower(x, a=0.5, lambda=1, delta=1, log=FALSE)
pasypower(x, a=0.5, lambda=1, delta=1, log.p=FALSE, lower.tail=TRUE)
varasypower(p, a=0.5, lambda=1, delta=1, log.p=FALSE, lower.tail=TRUE)
esasypower(p, a=0.5, lambda=1, delta=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first scale parameter, must be in the unit interval, the default is 0.5}
  \item{delta}{the value of the second scale parameter, must be positive, the default is 1}
  \item{lambda}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dasypower(x)
pasypower(x)
varasypower(x)
esasypower(x)}

