\name{rAitchison}
\alias{rAitchison}
\alias{dAitchison}
\alias{AitchisonDistributionIntegrals}
\title{Aitchison Distribution}
\description{
  The Aitchison distribution is a class of distributions the simplex,
  containing the normal and the Dirichlet as subfamilies.  
}
\usage{
dAitchison(x,
           theta=alpha+sigma \%*\% clr(mu),
           beta=-1/2*gsi.svdinverse(sigma),
           alpha=mean(theta),
           mu=clrInv(c(sigma\%*\%(theta-alpha))),
           sigma=-1/2*gsi.svdinverse(beta),
           grid=30,
           realdensity=FALSE,
           expKappa=AitchisonDistributionIntegrals(theta,beta,
                        grid=grid,mode=1)$expKappa)
rAitchison(n,
           theta=alpha+sigma \%*\% clr(mu),
           beta=-1/2*gsi.svdinverse(sigma),
           alpha=mean(theta),
           mu=clrInv(c(sigma\%*\%(theta-alpha))),
           sigma=-1/2*gsi.svdinverse(beta), withfit=FALSE)
AitchisonDistributionIntegrals(
           theta=alpha+sigma \%*\% clr(mu),
           beta=-1/2*gsi.svdinverse(sigma),
           alpha=mean(theta),
           mu=clrInv(c(sigma\%*\%(theta-alpha))),
           sigma=-1/2*gsi.svdinverse(beta),
           grid=30,
           mode=3)
}
\arguments{
  \item{x}{acomp-compositions the density should be computed for.}
  \item{n}{integer: number of datasets to be simulated}
  \item{theta}{numeric vector: Location parameter vector }
  \item{beta}{matrix: Spread parameter matrix (clr or ilr)}
  \item{alpha}{positiv scalar: departure from normality parameter (positive scalar)}
  \item{mu}{acomp-composition, normal reference mean parameter composition}
  \item{sigma}{ matrix: normal reference variance matrix (clr or ilr)}
  \item{grid}{ integer: number of discretisation points along each side of the simplex }
  \item{realdensity}{logical: if true the density is given with respect
    to the Haar measure of the real simplex, if false the density is
    given with respect to the Aitchison measure of the simplex.}
  \item{mode}{integer: desired output:
    -1: Compute nothing, only transform parameters,\cr
    0: Compute only oneIntegral and kappaIntegral,\cr
    1: compute also the clrMean,\cr
    2: compute also the clrSqExpectation,\cr
    3: same as 2, but compute clrVar instead of clrSqExpectation\cr
  }
  \item{expKappa}{The closing divisor of the density}
  \item{withfit}{Should a pre-spliting of the Aitchison density be used for 
    simulation?}
}
\value{
  \item{dAitchison}{Returns the density of the Aitchison distribution
    evaluated at x as a numeric vector.} 
  \item{rAitchison}{Returns a sample of size n of simulated compostions
    as an acomp object.} 
  \item{AitchisondistributionIntegrals}{Returns a list with
    \describe{
      \item{theta}{the theta parameter given or computed}
      \item{beta}{the beta parameter given or computed}
      \item{alpha}{the alpha parameter given or computed}
      \item{mu}{the mu parameter given or computed}
      \item{sigma}{the sigma parameter given or computed}
      \item{expKappa}{the integral over the density without closing
	constant. I.e. the inverse of the closing constant and the exp
	of \eqn{\kappa_{Ait(\theta,\beta)}}{kappa[Ait(theta,beta)]}}
      \item{kappaIntegral}{The expected value of the mean of the logs of
	the components as numerically computed}
      \item{clrMean}{The mean of the clr transformed random variable,
	computed numerically}
      \item{clrSqExpectation}{The expectation of \eqn{clr(X)clr(X)^t}
	computed numerically.}
      \item{clrVar}{The variance covariance matrix of clr(X),
	computed numerically}
      }
  } 
}
\details{
  The Aitchison Distribution is a joint generalisation of the Dirichlet
  Distribution and the additive log-normal distribution (or normal on the
  simplex). It can be parametrized by Ait(theta,beta) or by
  Ait(alpha,mu,Sigma). Actually, beta and Sigma can be easily
  transformed into each other, such that only one of them is
  necessary. Parameter theta is a vector in \eqn{R^D}, alpha is its sum, mu is a
  composition in \eqn{S^D}, and beta and sigma are symmetric matrices, which
  can either be expressed  in ilr or clr space.
  The parameters are transformed as
  \deqn{\beta=-1/2 \Sigma^{-1}}{beta=-1/2 Sigma^{-1}} 
  \deqn{\theta=clr(\mu)\Sigma+\alpha (1,\ldots,1)^t}{theta=clr(mu)Sigma+alpha 1}
  The distribution exists, if either,
  \eqn{\alpha\geq 0}{alpha >= 0} and Sigma is positive definite (or beta
  negative definite) in
  ilr-coordinates, or if each theta is strictly positive and Sigma has
  at least one positive eigenvalue (or beta has at least one negative
  eigenvalue). The simulation procedure currently only works with the
  first case. 
  \cr
  AitchisonDistributionIntegral is a convenience function to compute the
  parameter transformation and several functions of these
  parameters. This is done by numerical integration over a
  multinomial simplex lattice of D parts with \code{grid} many elements
  (see \code{\link[combinat]{xsimplex}}). 
  \cr
  The density of the Aitchison distribution is given by:
  \deqn{f(x,\theta,\beta)=exp((\theta-1)^t \log(x)+ilr(x)^t \beta ilr(x))/exp(\kappa_{Ait(\theta,\beta)})}{f(x,theta,beta)=exp((theta-1)^tlog(x)+ilr(x)^t beta ilr(x))/expKappa(theta,beta)}
  with respect to the classical Haar measure on the simplex, and as
  \deqn{f(x,\theta,\beta)=exp(\theta^t \log(x)+ilr(x)^t \beta
  ilr(x))/exp(\kappa_{Ait(\theta,\beta)})}{f(x,theta,beta)=exp(theta^t
  log(x)+ilr(x)^t beta ilr(x))/expKappa(theta,beta)} with respect to the Aitchison
  measure of the simplex. The closure constant expKappa is computed
  numerically, in \code{AitchisonDistributionIntegrals}.
  
  The random composition generation is done by rejection sampling based on
  an optimally fitted additive logistic normal distribution. Thus, it only
  works if the correponding Sigma in ilr would be positive definite.  
  
}
\note{The simulation procedure currently only works with a
  positive definite Sigma. You need a relatively high grid constant for
  precise values in the numerical integration.}
\author{K.Gerald v.d. Boogaart, R. Tolosana-Delgado \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{runif.acomp}}, \code{\link{rnorm.acomp}},
  \code{\link{rDirichlet.acomp}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

}
\examples{
(erg<-AitchisonDistributionIntegrals(c(-1,3,-2),ilrvar2clr(-diag(c(1,2))),grid=20))

(myvar<-with(erg, -1/2*ilrvar2clr(solve(clrvar2ilr(beta)))))
(mymean<-with(erg,myvar\%*\%theta))

with(erg,myvar-clrVar)
with(erg,mymean-clrMean)

}
\keyword{distribution}
