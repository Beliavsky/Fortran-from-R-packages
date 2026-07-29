\name{asylaplace}
\alias{dasylaplace}
\alias{pasylaplace}
\alias{varasylaplace}
\alias{esasylaplace}
\title{Asymmetric Laplace distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the asymmetric Laplace distribution due to Kotz et al. (2001) given by
\deqn{\begin{array}{ll}
&\displaystyle
\displaystyle
f(x) = \left\{ \begin{array}{ll}
\displaystyle
\frac {\kappa \sqrt{2}}{\tau \left( 1 + \kappa^2 \right)}
\exp \left( -\frac {\kappa \sqrt{2}}{\tau} \left | x - \theta \right | \right), & \mbox{if $x \geq \theta$,}
\\
\\
\displaystyle
\frac {\kappa \sqrt{2}}{\tau \left( 1 + \kappa^2 \right)}
\exp \left( -\frac {\sqrt{2}}{\kappa \tau} \left | x - \theta \right | \right), & \mbox{if $x < \theta$,}
\end{array}
\right.
\\
&\displaystyle
F (x) =  \left\{ \begin{array}{ll}
\displaystyle
1 - \frac {1}{1 + \kappa^2} \exp \left( \frac {\kappa \sqrt{2} (\theta - x)}{\tau} \right), &
\mbox{if $x \geq \theta$,}
\\
\\
\displaystyle
\frac {\kappa^2}{1 + \kappa^2} \exp \left( \frac {\sqrt{2} (x - \theta)}{\kappa \tau} \right), & \mbox{if $x < \theta$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
\theta - \frac {\tau}{\sqrt{2} \kappa}
\log \left[ (1 - p) \left( 1 + \kappa^2 \right) \right], & \mbox{if $p \geq \frac {\kappa^2}{1 + \kappa^2}$,}
\\
\\
\displaystyle
\theta + \frac {\kappa \tau}{\sqrt{2}} \log \left[ p \left( 1 + \kappa^{-2} \right) \right], &
\mbox{if $p < \frac {\kappa^2}{1 + \kappa^2}$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
\frac {\theta}{p} + \theta - \frac {\tau}{\sqrt{2} \kappa} \log \left( 1 + \kappa^2 \right) +
\frac {\sqrt{2} \tau \left( 1 + 2 \kappa^2 \right)}{2 \kappa \left( 1 + \kappa^2 \right) p}
\log \left( 1 + \kappa^2 \right)
\\
\displaystyle
\quad
-\frac {\sqrt{2} \tau \kappa \log \kappa}{\left( 1 + \kappa^2 \right) p} -
\frac {\theta \kappa^2}{\left( 1 + \kappa^2 \right) p} +
\frac {\tau \left( 1 - \kappa^4 \right)}{\sqrt{2} \kappa \left( 1 + \kappa^2 \right) p}
\\
\displaystyle
\quad
-\frac {\tau (1 - p)}{\sqrt{2} \kappa p} + \frac {\tau (1 - p)}{\sqrt{2} \kappa p} \log (1 - p), &
\mbox{if $p \geq \frac {\kappa^2}{1 + \kappa^2}$,}
\\
\\
\displaystyle
\theta + \frac {\kappa \tau}{\sqrt{2}}
\log \left( 1 + \kappa^{-2} \right) +
\frac {\kappa \tau}{\sqrt{2}} (\log p - 1), & \mbox{if $p < \frac {\kappa^2}{1 + \kappa^2}$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{-\infty < \theta < \infty}, the location parameter, \eqn{\kappa > 0}, the first scale parameter, and \eqn{\tau > 0}, the second scale parameter.}
\usage{
dasylaplace(x, tau=1, kappa=1, theta=0, log=FALSE)
pasylaplace(x, tau=1, kappa=1, theta=0, log.p=FALSE, lower.tail=TRUE)
varasylaplace(p, tau=1, kappa=1, theta=0, log.p=FALSE, lower.tail=TRUE)
esasylaplace(p, tau=1, kappa=1, theta=0)
}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{theta}{the value of the location parameter, can take any real value, the default is zero}
  \item{kappa}{the value of the first scale parameter, must be positive, the default is 1}
  \item{tau}{the value of the second scale parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}
}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
dasylaplace(x)
pasylaplace(x)
varasylaplace(x)
esasylaplace(x)}
 
