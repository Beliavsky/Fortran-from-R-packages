\name{triangular}
\alias{dtriangular}
\alias{ptriangular}
\alias{vartriangular}
\alias{estriangular}
\title{Triangular distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the triangular distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \left\{
\begin{array}{ll}
\displaystyle
0, & \mbox{if $x < a$,}
\\
\\
\displaystyle
\frac {2 (x - a)}{(b - a) (c - a)}, & \mbox{if $a \leq x \leq c$,}
\\
\\
\displaystyle
\frac {2 (b - x)}{(b - a) (b - c)}, & \mbox{if $c < x \leq b$,}
\\
\\
\displaystyle
0, & \mbox{if $b < x$,}
\end{array} \right.
\\
&\displaystyle
F (x) = \left\{
\begin{array}{ll}
\displaystyle
0, & \mbox{if $x < a$,}
\\
\\
\displaystyle
\frac {(x - a)^2}{(b - a) (c - a)}, & \mbox{if $a \leq x \leq c$,}
\\
\\
\displaystyle
1 - \frac {(b - x)^2}{(b - a) (b - c)}, & \mbox{if $c < x \leq b$,}
\\
\\
\displaystyle
1, & \mbox{if $b < x$,}
\end{array} \right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{
\begin{array}{ll}
\displaystyle
a + \sqrt{p (b - a) (c - a)}, & \mbox{if $0 < p < \frac {c - a}{b - a}$,}
\\
\\
\displaystyle
b - \sqrt{(1 - p) (b - a) (b - c)}, & \mbox{if $\frac {c - a}{b - a} \leq p < 1$,}
\end{array} \right.
\\
&\displaystyle
{\rm ES}_p (X) =
\left\{
\begin{array}{ll}
\displaystyle
a + \frac {2}{3} \sqrt{p (b - a) (c - a)}, & \mbox{if $0 < p < \frac {c - a}{b - a}$,}
\\
\\
\displaystyle
b + \frac {a - c}{p} + \frac {2 (2 c - a - b)}{3 p} +2 \sqrt{(b - a) (b - c)} \frac {(1 - p)^{3/2}}{3 p}, &
\mbox{if $\frac {c - a}{b - a} \leq p < 1$}
\end{array} \right.
\end{array}}
for \eqn{a \leq x \leq b}, \eqn{0 < p < 1}, \eqn{-\infty < a < \infty}, the  first location parameter, \eqn{-\infty < a < c < \infty}, the second location parameter, and \eqn{-\infty < c < b < \infty}, the  third location parameter.}
\usage{
dtriangular(x, a=0, b=2, c=1, log=FALSE)
ptriangular(x, a=0, b=2, c=1, log.p=FALSE, lower.tail=TRUE)
vartriangular(p, a=0, b=2, c=1, log.p=FALSE, lower.tail=TRUE)
estriangular(p, a=0, b=2, c=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{a}{the value of the first location parameter, can take any real value, the default is zero}
  \item{c}{the value of the second location parameter, can take any real value but must be greater than a, the default is 1}
  \item{b}{the value of the third location parameter, can take any real value but must be greater than c, the default is 2}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dtriangular(x)
ptriangular(x)
vartriangular(x)
estriangular(x)}
 
