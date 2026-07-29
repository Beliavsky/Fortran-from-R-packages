\name{exppower}
\alias{dexppower}
\alias{pexppower}
\alias{varexppower}
\alias{esexppower}
\title{Exponential power distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the exponential power distribution due to Subbotin (1923) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\displaystyle 2 a^{1/a} \sigma \Gamma \left( 1 + 1/a \right)}
\exp \left\{ -\frac {\mid x - \mu \mid^a}{a \sigma^a} \right\},
\\
&\displaystyle
F (x) =
\left\{
\begin{array}{ll}
\displaystyle
\frac {1}{2} Q \left( \frac {1}{a}, \frac {(\mu - x)^a}{a \sigma^a} \right), & \mbox{if $x \leq \mu$,}
\\
\\
\displaystyle
1 - \frac {1}{2} Q \left( \frac {1}{a}, \frac {(x - \mu)^a}{a \sigma^a} \right), & \mbox{if $x > \mu$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\mu - a^{1/a} \sigma \left[ Q^{-1} \left( \frac {1}{a}, 2 p \right) \right]^{1/a}, & \mbox{if $p \leq 1/2$,}
\\
\\
\mu + a^{1/a} \sigma \left[ Q^{-1} \left( \frac {1}{a}, 2 (1 - p) \right) \right]^{1/a}, & \mbox{if $p > 1/2$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\mu - \frac {a^{1/a} \sigma}{p} \int_0^p \left[ Q^{-1} \left( \frac {1}{a}, 2 v \right) \right]^{1/a} dv, & \mbox{if $p \leq 1/2$,}
\\
\\
\displaystyle
\mu - \frac {a^{1/a} \sigma}{p} \int_0^{1/2} \left[ Q^{-1} \left( \frac {1}{a}, 2 v \right) \right]^{1/a} dv
\\
\displaystyle
\quad
+\frac {a^{1/a} \sigma}{p} \int_{1/2}^p \left[ Q^{-1} \left( \frac {1}{a}, 2 (1 - v) \right) \right]^{1/a} dv, & \mbox{if $p > 1/2$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, \eqn{\sigma > 0}, the scale parameter, and
\eqn{a > 0}, the shape parameter.}
\usage{
dexppower(x, mu=0, sigma=1, a=1, log=FALSE)
pexppower(x, mu=0, sigma=1, a=1, log.p=FALSE, lower.tail=TRUE)
varexppower(p, mu=0, sigma=1, a=1, log.p=FALSE, lower.tail=TRUE)
esexppower(p, mu=0, sigma=1, a=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dexppower(x)
pexppower(x)
varexppower(x)
esexppower(x)}
 
