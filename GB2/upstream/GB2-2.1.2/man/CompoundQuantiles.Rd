\name{CompoundQuantiles}
\Rdversion{2.1}
\alias{CompoundQuantiles}
\alias{qcgb2}
\alias{rcgb2}
\title{
Quantiles and random generation of the Compound Distribution based on the GB2
}
\description{
Calculation of the quantiles of a compound GB2 random variable. Random generation of compound GB2 variables.
}
\usage{
qcgb2(prob, shape1, scale, shape2, shape3, pl0, pl, decomp="r", tol=1e-08, ff=1.5, 
debug=FALSE, maxiter=50)
rcgb2(n, shape1, scale, shape2, shape3, pl0, pl, decomp="r", tol=1e-02, maxiter=100,
debug = FALSE)
}
\arguments{
  \item{prob}{numeric; vector of probabilities between 0 and 1.}
  \item{shape1,scale,shape2,shape3}{numeric; positive parameters of the GB2 distribution.}
  \item{n}{numeric; number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; a vector of mixture probabilities. Sums to one. If \eqn{pl=pl0} we obtain the GB2 distribution.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution.
   By default, \code{decomp} = "r" - right tail decomposition.}
  \item{ff}{numeric; a tuning parameter.}
  \item{debug}{logical; By default, \code{debug = FALSE}.}
  \item{maxiter}{numeric; maximum number of iterations to perform.}
  \item{tol}{numeric; tolerance with default 0, meaning to iterate until additional terms do not change the partial sum.}
}
\value{ 
\code{qcgb2} returns a vector of quantiles and \code{rcgb2} return a vector of size \code{n} of GB2 compound random deviates.
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\author{
Monique Graf and Desislava Nedyalkova
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

#Quantiles

qcgb2(0.5,af,bf,pf,qf,p0,p1)
qcgb2(1,af,bf,pf,qf,p0,p1)
qcgb2(c(0.5,0.8),af,bf,pf,qf,p0,p1)

#Random generation
rcgb2(10,af,bf,pf,qf,p0,p1)

#}
}
\keyword{distribution}

