\name{price.mln.option}
\alias{price.mln.option}
\title{Price Options on Mixture of Lognormals}
\description{
  \code{mln.option.price} gives the price of a call and a put option at a strike when the risk neutral density is a mixture of two lognormals.
}
\usage{
price.mln.option(r, te, y, k, alpha.1, meanlog.1, meanlog.2, sdlog.1, sdlog.2)
}
\arguments{
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{y}{dividend yield}
  \item{k}{strike}
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
  \item{call }{call price}
  \item{put }{put price}
  \item{s0 }{current value of the asset as implied by the mixture distribution}
}
\references{

  F. Gianluca and A. Roncoroni (2008)
  \emph{Implementing Models in Quantitative Finance: Methods and Cases}

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
# Try out a range of options
#

r  = 0.05
te = 60/365
k  = 700:1300
y  = 0.02
meanlog.1 = 6.80
meanlog.2 = 6.95
sdlog.1   = 0.065
sdlog.2   = 0.055
alpha.1   = 0.4


mln.prices = price.mln.option(r = r, y = y, te = te, k = k, alpha.1 = alpha.1, 
  meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)

par(mfrow=c(1,2))
plot(mln.prices$call ~ k)
plot(mln.prices$put  ~ k)
par(mfrow=c(1,1))

}

