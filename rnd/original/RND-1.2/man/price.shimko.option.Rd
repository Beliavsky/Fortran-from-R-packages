\name{price.shimko.option}
\alias{price.shimko.option}
\title{Price Option based on Shimko's Method}
\description{
  \code{price.shimko.option} prices a European option based on the extracted Shimko volatility function.
}
\usage{
price.shimko.option(r, te, s0, k, y, a0, a1, a2)
}
\arguments{
  \item{r}{risk free rate}
  \item{te}{time to expiration}
  \item{s0}{current asset value}
  \item{k}{strike}
  \item{y}{dividend yield}
  \item{a0}{constant term in the quadratic polyynomial}
  \item{a1}{coefficient term of k in the quadratic polynomial}
  \item{a2}{coefficient term of k squared in the quadratic polynomial}
}
\details{This function may produce negative option values when nonsensical
 values are used for a0, a1, and a2.}
\value{
  \item{call }{call price}
  \item{put }{put price}
}
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

r       =  0.05
y       =  0.02
te      =  60/365
s0      =  1000
k       =  950
sigma   =  0.25
a0      =  0.30
a1      =  -0.00387
a2      =  0.00000445

#
#  Note how Shimko price is the same when a0 = sigma, a1=a2=0 but substantially 
#  more when a0, a1, a2 are changed so the implied volatilies are very high!
#

price.bsm.option(r = r, te = te, s0 = s0, k = k, sigma = sigma, y = y)$call
price.shimko.option(r = r, te = te, s0 = s0, k = k, y = y, 
                    a0 = sigma, a1 = 0, a2 = 0)$call
price.shimko.option(r = r, te = te, s0 = s0, k = k, y = y, 
                    a0 = a0, a1 = a1, a2 = a2)$call
}

