\name{tsp}
\alias{dtsp}
\alias{ptsp}
\alias{vartsp}
\alias{estsp}
\title{Two sided power distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the two sided power distribution due to van Dorp and Kotz (2002) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \left\{ \begin{array}{ll}
\displaystyle
a \left( \frac {x}{\theta} \right)^{a - 1}, & \mbox{if $0 < x \leq \theta$,}
\\
\displaystyle
a \left( \frac {1 - x}{1 - \theta} \right)^{a - 1}, & \mbox{if $\theta < x < 1$,}
\end{array}
\right.
\\
&\displaystyle
F (x) = \left\{ \begin{array}{ll}
\displaystyle
\theta \left( \frac {x}{\theta} \right)^a, & \mbox{if $0 < x \leq \theta$,}
\\
\displaystyle
1 - (1 - \theta) \left( \frac {1 - x}{1 - \theta} \right)^a, & \mbox{if $\theta < x < 1$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) =  \left\{ \begin{array}{ll}
\displaystyle
\theta \left( \frac {p}{\theta} \right)^{1 / a}, & \mbox{if $0 < p \leq \theta$,}
\\
\displaystyle
1 - (1 - \theta) \left( \frac {1 - p}{1 - \theta} \right)^{1 / a}, & \mbox{if $\theta < p < 1$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
\frac {a \theta}{a + 1} \left( \frac {p}{\theta} \right)^{1 / a}, & \mbox{if $0 < p \leq \theta$,}
\\
\displaystyle
1 - \frac {\theta}{p} + \frac {a (2 \theta - 1)}{(a + 1) p} + \frac {a (1 - \theta)^2}{(a + 1) p}
\left( \frac {1 - p}{1 - \theta} \right)^{1 + 1 / a}, & \mbox{if $\theta < p < 1$}
\end{array}
\right.
\end{array}}
for \eqn{0 < x < 1}, \eqn{0 < p < 1}, \eqn{a > 0}, the shape parameter, and \eqn{-\infty < \theta < \infty}, the location parameter.}
\usage{
dtsp(x, a=1, theta=0.5, log=FALSE)
ptsp(x, a=1, theta=0.5, log.p=FALSE, lower.tail=TRUE)
vartsp(p, a=1, theta=0.5, log.p=FALSE, lower.tail=TRUE)
estsp(p, a=1, theta=0.5)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the location parameter, must take a value in the unit interval, the default is 0.5}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dtsp(x)
ptsp(x)
vartsp(x)
estsp(x)}
 
