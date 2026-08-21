\name{powerofpsdmatrix}
\alias{powerofpsdmatrix}
\title{power transform of a matrix}
\description{
  Computes the power of a positive semi-definite symmetric matrix.
}
\usage{
          powerofpsdmatrix( M , p,\dots)
          }
\arguments{
  \item{M}{a matrix, preferably symmetric}
  \item{p}{a single number giving the power}
  \item{\dots}{further arguments to the singular value decomposition}
}
\value{
  \code{U\%*\% D^p \%*\% t(P)} where the \code{UDP} is the singular value
  decomposition of M. 
}
\details{
  for a symmetric matrix the computed result can actually be considered as
  a version of the given power of the matrix fullfilling the relation:
  \deqn{M^pM^q=M^{p+q}}
  The symmetry of the matrix is not checked.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
data(SimulatedAmounts)
d <- ilr(sa.lognormals)
var( d \%*\% powerofpsdmatrix(var(d),-1/2)) # Unit matrix
}
\keyword{multivariate}
