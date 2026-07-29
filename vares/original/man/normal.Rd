\name{normal}
\alias{dnormal}
\alias{pnormal}
\alias{varnormal}
\alias{esnormal}
\title{Normal distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the normal distribution due to de Moivre (1738) and Gauss (1809) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) = \frac {1}{\sigma} \phi \left( \frac {x - \mu}{\sigma} \right),
\\
&\displaystyle
F (x) = \Phi \left( \frac {x - \mu}{\sigma} \right),
\\
&\displaystyle
{\rm VaR}_p (X) = \mu + \sigma \Phi^{-1} (p),
\\
&\displaystyle
{\rm ES}_p (X) = \mu + \frac {\sigma}{p} \int_0^p \Phi^{-1} (v) dv
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \mu < \infty}, the location parameter, and \eqn{\sigma > 0}, the scale parameter,
where \eqn{\phi (\cdot)} denotes the pdf of a standard normal random variable,
and \eqn{\Phi (\cdot)} denotes the cdf of a standard normal random variable.}
\usage{
dnormal(x, mu=0, sigma=1, log=FALSE)
pnormal(x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
varnormal(p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
esnormal(p, mu=0, sigma=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{mu}{the value of the location parameter, can take any real value, the default is zero}
  \item{sigma}{the value of the scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dnormal(x)
pnormal(x)
varnormal(x)
esnormal(x)}
 
