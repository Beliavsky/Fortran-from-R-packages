\name{dmln}
\alias{dmln}
\title{Density of Mixture Lognormal}
\description{\code{mln} is the probability density function of a mixture of two lognormal densities.}
\usage{
dmln(x, alpha.1, meanlog.1, meanlog.2, sdlog.1, sdlog.2)
}
\arguments{
  \item{x}{value at which the denisty is to be evaluated}
  \item{alpha.1}{proportion of the first lognormal.  Second one is 1 - \code{alpha.1}}
  \item{meanlog.1}{mean of the log of the first lognormal}
  \item{meanlog.2}{mean of the log of the second lognormal}
  \item{sdlog.1}{standard deviation of the log of the first lognormal}
  \item{sdlog.2}{standard deviation of the log of the second lognormal}
}
\details{
  mln is the density f(x) = alpha.1 * g(x) + (1 - alpha.1) * h(x), where g and h are densities of two lognormals with
  parameters (mean.log.1, sdlog.1) and (mean.log.2, sdlog.2) respectively.
}
\value{
  \item{out }{density value at x}
}
\references{
  B. Bahra (1996):
  Probability distribution of future asset prices implied by option prices.
  \emph{Bank of England Quarterly Bulletin}, August 1996, 299-311

  P. Soderlind and L.E.O. Svensson (1997)
  New techniques to extract market expectations from financial instruments.
  \emph{Journal of Monetary Economics}, 40, 383-429

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{

#
# A bimodal risk neutral density!
#

mln.alpha.1   = 0.4
mln.meanlog.1 = 6.3
mln.meanlog.2 = 6.5
mln.sdlog.1   = 0.08
mln.sdlog.2   = 0.06

k  = 300:900
dx = dmln(x = k, alpha.1 = mln.alpha.1, meanlog.1 = mln.meanlog.1, 
         meanlog.2 = mln.meanlog.2, 
         sdlog.1 = mln.sdlog.1, sdlog.2 = mln.sdlog.2)
plot(dx ~ k, type="l")
}

