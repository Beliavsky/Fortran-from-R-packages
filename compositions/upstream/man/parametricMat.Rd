\name{parametricMat}
\alias{parametricPosdefMat}
\alias{parametricRank1Mat}
\alias{parameterPosdefMat}
\alias{parameterRank1Mat}
\alias{parametricPosdefClrMat}
\alias{parametricRank1ClrMat}
\alias{parameterPosdefClrMat}
\alias{parameterRank1ClrMat}
\title{Unique parametrisations for matrices.}
\description{
  Helper functions to parametrize positive semidefinite matrices in
  multivariate variogram models.
}
\usage{
parametricRank1Mat(p)
parametricPosdefMat(p)
parameterRank1Mat(A)
parameterPosdefMat(A)
parametricRank1ClrMat(p)
parametricPosdefClrMat(p)
parameterRank1ClrMat(A)
parameterPosdefClrMat(A)
          }
\arguments{
  \item{A}{a positiv definit matrix of the given type}
  \item{p}{a vector of parameters describing the matrix, as returned by
    the parameter functions.}
}
\value{
  A or p, depending on what is not given.
}
\details{
  The rank 1 matrix is parametrised by the first eigenvector scaled by
  the square root of the eigenvalue. The positiv semidefinit matrix the
  entries of a upper right triangular matrix R with
  \code{t(R)\%*\%R==A}. The clr matrices are work with the parameters of
  the corresponding ilr matrix. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}},
  \code{\link{CompLinModCoReg}},
  \code{\link{vgmFit}}
}
\examples{
parametricRank1Mat(c(0,0,2))
parametricPosdefMat(c(0,0,1,0,0,0))
parameterRank1Mat(matrix(1,nr=3,nc=3))
parameterPosdefMat(diag(5))
}
\keyword{multivariate}
