\name{F}
\alias{dF}
\alias{pF}
\alias{varF}
\alias{esF}
\title{F distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the F distribution given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {1}{B \left( \frac {d_1}{2}, \frac {d_2}{2} \right)}
\left( \frac {d_1}{d_2} \right)^{\frac {d_1}{2}}
x^{\frac {d_1}{2} - 1}
\left( 1 + \frac {d_1}{d_2} x \right)^{-\frac {d_1 + d_2}{2}},
\\
&\displaystyle
F (x) = I_{\frac {d_1 x}{d_1 x + d_2}} \left( \frac {d_1}{2}, \frac {d_2}{2} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \frac {d_2}{d_1} \frac {I_p^{-1} \left( \frac {d_1}{2}, \frac {d_2}{2} \right)}
{1 - I_p^{-1} \left( \frac {d_1}{2}, \frac {d_2}{2} \right)},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {d_2}{d_1 p} \int_0^p
\frac {I_v^{-1} \left( \frac {d_1}{2}, \frac {d_2}{2} \right)}
{1 - I_v^{-1} \left( \frac {d_1}{2}, \frac {d_2}{2} \right)} dv
\end{array}}
for \eqn{x \geq K}, \eqn{0 < p < 1}, \eqn{d_1 > 0}, the first degree of freedom parameter,
and \eqn{d_2 > 0}, the second degree of freedom parameter.}
\usage{
dF(x, d1=1, d2=1, log=FALSE)
pF(x, d1=1, d2=1, log.p=FALSE, lower.tail=TRUE)
varF(p, d1=1, d2=1, log.p=FALSE, lower.tail=TRUE)
esF(p, d1=1, d2=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{d1}{the value of the first degree of freedom parameter, must be positive, the default is 1}
  \item{d2}{the value of the second degree of freedom parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dF(x)
pF(x)
varF(x)
esF(x)}
 
