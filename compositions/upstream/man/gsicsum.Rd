\name{gsi.csum}
\alias{gsi.csum}
\alias{gsi.rsum}
\title{Internal function: row and column sums of matrices}
\description{
  Abbreviation for collapsing matrices
}
\usage{
gsi.csum(x)
gsi.rsum(x)
}
\arguments{
  \item{x}{the matrix}
}
\value{
  A numeric vector.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\details{
  gsi.csum computes the sums of all the columns of x and gsi.rsum
  computes the sums of the rows of x. For convenience, only finite values
  are added.
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
 \code{\link{sum}}
}
\examples{
#A <- matrix(c(0,1,0,0,0,0),ncol=2)
#b <- diag(3)
#erg <- gsi.svdsolve(A,b)
#erg
#A %*% erg 
#diag(c(0,1,0))  # richtig
}
\keyword{internal}
