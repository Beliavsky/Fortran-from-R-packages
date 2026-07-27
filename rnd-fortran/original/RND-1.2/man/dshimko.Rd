\name{dshimko}
\alias{dshimko}
\title{Density Implied by Shimko Method}
\description{
  \code{dshimko} is the probability density function implied by the Shimko method.
}
\usage{
dshimko(r, te, s0, k, y, a0, a1, a2)
}
\arguments{
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{s0}{current asset value}
  \item{k}{strike at which volatility to be computed}
  \item{y}{dividend yield}
  \item{a0}{constant term in the quadratic polynomial}
  \item{a1}{coefficient term of k in the quadratic polynomial}
  \item{a2}{coefficient term of k squared in the quadratic polynomial}
}
\details{
The implied volatility is modeled as: \eqn{\sigma(k) = a_0 + a_1 k + a_2 k^2}
}
\value{density value at x}
\references{
  D. Shimko (1993)
  Bounds of probability.
  \emph{Risk}, 6, 33-47

  E. Jondeau and S. Poon and M. Rockinger (2007):
  \emph{Financial Modeling Under Non-Gaussian Distributions} 
  Springer-Verlag, London

}
\author{Kam Hamidieh}
\examples{

#
# a0, a1, a2 values come from Shimko's paper.
#

r     =  0.05
y     =  0.02
a0    =  0.892
a1    = -0.00387
a2    =  0.00000445
te    =  60/365
s0    =  400
k     =  seq(from = 250, to = 500, by = 1)
sigma =  0.15

#
# Does it look like a proper density and intergate to one?
#

dx = dshimko(r = r, te = te, s0 = s0, k = k, y = y, a0 = a0, a1 = a1, a2 = a2)
plot(dx ~ k, type="l")

#
# sum(dx) should be about 1 since dx is a density.
#

sum(dx)

}

