\name{loglaplace}
\alias{dloglaplace}
\alias{ploglaplace}
\alias{varloglaplace}
\alias{esloglaplace}
\title{Log Laplace distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the log Laplace distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \left\{
\begin{array}{ll}
\displaystyle
\frac {a b x^{b - 1}}{\delta^b (a + b)}, & \mbox{if $x \leq \delta$,}
\\
\\
\displaystyle
\frac {a b \delta^a}{x^{a + 1} (a + b)}, & \mbox{if $x > \delta$,}
\end{array}
\right.
\\
&\displaystyle
F (x) = \left\{
\begin{array}{ll}
\displaystyle
\frac {a x^b}{\delta^b (a + b)}, & \mbox{if $x \leq \delta$,}
\\
\\
\displaystyle
1 - \frac {b \delta^a}{x^a (a + b)}, & \mbox{if $x > \delta$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{
\begin{array}{ll}
\displaystyle
\delta \left[ p \frac {a + b}{a} \right]^{1/b}, & \mbox{if $p \leq \frac {a}{a + b}$,}
\\
\\
\displaystyle
\delta \left[ (1 - p) \frac {a + b}{a} \right]^{-1/a}, & \mbox{if $p > \frac {a}{a + b}$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
\frac {\delta b}{b + 1}
\left[ p \frac {a + b}{a} \right]^{1/b}, & \mbox{if $p \leq \frac {a}{a + b}$,}
\\
\\
\displaystyle
\frac {a \delta}{p (1 + 1/b) (a + b)} +
\frac {a^{1/a} b^{1 - 1/a} \delta}{p (a + b) (1 - 1/a)}
\\
\displaystyle
\quad
-\frac {\delta (1 - p)}{p (1 - 1/a)}
\left[ \frac {a}{(a + b) (1 - p)} \right]^{1/a}, &
\mbox{if $p > \frac {a}{a + b}$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{\delta > 0}, the scale parameter,
\eqn{a > 0}, the first shape parameter, and \eqn{b > 0}, the second shape parameter.}
\usage{
dloglaplace(x, a=1, b=1, delta=0, log=FALSE)
ploglaplace(x, a=1, b=1, delta=0, log.p=FALSE, lower.tail=TRUE)
varloglaplace(p, a=1, b=1, delta=0, log.p=FALSE, lower.tail=TRUE)
esloglaplace(p, a=1, b=1, delta=0)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{delta}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{b}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dloglaplace(x)
ploglaplace(x)
varloglaplace(x)
esloglaplace(x)}
 
