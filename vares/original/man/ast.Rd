\name{ast}
\alias{dast}
\alias{past}
\alias{varast}
\alias{esast}
\title{Generalized asymmetric Student's t distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the generalized asymmetric Student's \eqn{t} distribution due to Zhu and Galbraith  (2010) given by
\deqn{\begin{array}{ll}
&\displaystyle
\displaystyle
f (x) = \left\{
\begin{array}{ll}
\displaystyle
\frac {\alpha}{\alpha^{*}} K \left( \nu_1 \right) \left[ 1 + \frac {1}{\nu_1}
\left( \frac {x}{2 \alpha^{*}} \right)^2 \right]^{-\frac {\nu_1 + 1}{2}}, & \mbox{if $x \leq 0$,}
\\
\\
\displaystyle
\frac {1 - \alpha}{1 - \alpha^{*}} K \left( \nu_2 \right)
\left[ 1 + \frac {1}{\nu_2} \left( \frac {x}{2 \left( 1 - \alpha^{*} \right)} \right)^2 \right]^{-\frac {\nu_2 + 1}{2}}, &
\mbox{if $x > 0$,}
\end{array}
\right.
\\
&\displaystyle
\displaystyle
F (x) = 2 \alpha F_{\nu_1} \left( \frac {\min (x, 0)}{2 \alpha^{*}} \right) -1 + \alpha + 2 (1 - \alpha)
F_{\nu_2} \left( \frac {\max (x, 0)}{2 - 2 \alpha^{*}} \right),
\\
&\displaystyle
\displaystyle
{\rm VaR}_p (X) = 2 \alpha^{*} F_{\nu_1}^{-1} \left( \frac {\min (p, \alpha)}{2 \alpha} \right) +
2 \left( 1 - \alpha^{*} \right) F_{\nu_2}^{-1} \left( \frac {\max (p, \alpha) + 1 - 2 \alpha}{2 - 2 \alpha} \right),
\\
&\displaystyle
\displaystyle
{\rm ES}_p (X) = \frac {2 \alpha^{*}}{p} \int_0^p F_{\nu_1}^{-1}
\left( \frac {\min (v, \alpha)}{2 \alpha} \right) dv +
\frac {2 \left( 1 - \alpha^{*} \right)}{p} \int_0^p F_{\nu_2}^{-1}
\left( \frac {\max (v, \alpha) + 1 - 2 \alpha}{2 - 2 \alpha} \right) dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{0 < \alpha < 1}, the scale parameter, \eqn{\nu_1 > 0}, the first degree of freedom  parameter,
and \eqn{\nu_2 > 0}, the second degree of freedom parameter, where \eqn{\alpha^{*} = \alpha K \left( \nu_1 \right) / \left\{ \alpha K \left( \nu_1 \right) + (1 - \alpha) K \left( \nu_2 \right) \right\}}, \eqn{K (\nu) = \Gamma \left( (\nu + 1)/2 \right) / \left[ \sqrt{\pi \nu} \Gamma (\nu/2)  \right]}, \eqn{F_\nu(\cdot)} denotes the cdf of a Student's \eqn{t} random variable with \eqn{\nu} degrees  of freedom, and \eqn{F_\nu^{-1} (\cdot)} denotes the inverse of \eqn{F_\nu(\cdot)}.}
\usage{
dast(x, nu1=1, nu2=1, alpha=0.5, log=FALSE)
past(x, nu1=1, nu2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
varast(p, nu1=1, nu2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
esast(p, nu1=1, nu2=1, alpha=0.5)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{alpha}{the value of the scale parameter, must be in the unit interval, the default is 0.5}
  \item{nu1}{the value of the first degree of freedom parameter, must be positive, the default is 1}
  \item{nu2}{the value of the second degree of freedom parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dast(x)
past(x)
varast(x)
esast(x)}
 
