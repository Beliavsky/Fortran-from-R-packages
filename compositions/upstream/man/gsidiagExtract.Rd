\name{gsi.diagExtract}
\alias{gsi.diagExtract}
\title{Internal functions: Get the diagonal of a matrix}
\description{
  Get the main diagonal of a matrix, even if the matrix is 1x1.
}
\usage{
gsi.diagExtract(x)
}
\arguments{
  \item{x}{ a matrix}
}
\value{
  a vector containing the main diagonal entries of x.
}
\details{
  The difference to original \code{\link[base]{diag}} is that it always
  gives the diagonal and does nothing flawed in case of a 1x1 matrix or
  a single number considered as such matrix.
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}

\seealso{
  \code{\link{gsi.diagGenerate}}, \code{\link[base]{diag}}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
#data(SimulatedAmounts)
#gsi.diagExtract(var(acomp(sa.lognormals,c(1,2))))
#gsi.diagExtract(var(ilr(acomp(sa.lognormals,c(1,2)))))
#gsi.diagExtract(var(ilt(aplus(sa.lognormals,c(1)))))
}
\keyword{internal}
