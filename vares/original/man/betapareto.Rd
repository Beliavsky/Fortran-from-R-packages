\name{betapareto}
\alias{dbetapareto}
\alias{pbetapareto}
\alias{varbetapareto}
\alias{esbetapareto}
\title{Beta Pareto distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the beta Pareto distribution due to Akinsete et al. (2008) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \frac {a K^{ad} x^{-ad - 1}}{B (c, d)}
\left[ 1 - \left( \frac {K}{x} \right)^a \right]^{c - 1},
\\
&\displaystyle
F (x) =  I_{1 - \left( \frac {K}{x} \right)^a} (c, d),
\\
&\displaystyle
{\rm VaR}_p (X) =  K \left[ 1 - I_p^{-1} (c, d) \right]^{-1 / a},
\\
&\displaystyle
{\rm ES}_p (X) = \frac {K}{p} \int_0^p \left[ 1 - I_v^{-1} (c, d) \right]^{-1 / a} dv
\end{array}}
for \eqn{x \geq K}, \eqn{0 < p < 1}, \eqn{K > 0}, the scale parameter, \eqn{a > 0}, the first shape parameter, \eqn{c > 0}, the second shape parameter, and \eqn{d > 0}, the third shape parameter.}
\usage{
dbetapareto(x, K=1, a=1, c=1, d=1, log=FALSE)
pbetapareto(x, K=1, a=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
varbetapareto(p, K=1, a=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
esbetapareto(p, K=1, a=1, c=1, d=1)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{K}{the value of the scale parameter, must be positive, the default is 1}
  \item{a}{the value of the first shape parameter, must be positive, the default is 1}
  \item{c}{the value of the second shape parameter, must be positive, the default is 1}
  \item{d}{the value of the third shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dbetapareto(x)
pbetapareto(x)
varbetapareto(x)
esbetapareto(x)}
 
