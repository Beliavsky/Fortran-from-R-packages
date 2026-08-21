\name{CompoundMoments}
\Rdversion{2.1}
\alias{mkl.cgb2}
\alias{moment.cgb2}
\alias{incompl.cgb2}
\title{
Moments of the Compound Distribution based on the GB2
}
\description{
These functions calculate the moment of order \code{k} and incomplete moment of order \code{k} of a GB2 compound random variable \eqn{X} as well as the moment of order \code{k} for each component density.
}
\usage{
mkl.cgb2(k, shape1, scale, shape2, shape3, pl0, decomp="r")
moment.cgb2(k, shape1, scale, shape2 ,shape3, pl0, pl, decomp="r")
incompl.cgb2(x, k, shape1, scale, shape2, shape3, pl0, pl, decomp="r")
}
\arguments{
  \item{x}{numeric; vector of quantiles.}
  \item{k}{numeric; order of the moment.}
  \item{shape1,scale,shape2,shape3}{numeric; positive parameters of the GB2 distribution.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; a vector of mixture probabilities. Sums to one. If \eqn{pl=pl0} we obtain the GB2 distribution.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution.
   By default, \code{decomp} = "r" - right tail decomposition.}
}
\value{ 
\code{mkl.cgb2} returns a vector of the moments of the component densities, \code{moment.cgb2} returns the moment of order \code{k} and \code{incompl.cgb2} - the incomplete moment of order \code{k}.
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\author{
Monique Graf
}
\examples{
#\dontrun{

#\library{cubature}

# GB2 parameters
af <- 5
bf <- 20000
pf <- 0.45 
qf <- 0.75

p0 <- rep(1/3,3)
p1 <- c(0.37,0.43,0.2)

# moments for the component densities
mkl.cgb2(1,af,bf,pf,qf,p0)
mkl.cgb2(-1,af,bf,pf,qf,p0,"l")

#Moment of order k
moment.cgb2(0.5,af,bf,pf,qf,p0,p1)
moment.cgb2(0.5,af,bf,pf,qf,p0,p1,"l")

#Incomplete moment of order k
incompl.cgb2(20000,1,af,bf,pf,qf,p0,p1)
incompl.cgb2(20000,1,af,bf,pf,qf,p0,p1,"l")
#}
}
\keyword{distribution}

