\name{gsi.diagGenerate}
\alias{gsi.diagGenerate}
\title{Internal functions: Generate a diagonal matrix }
\description{
  Generate a diagonal matrix from a vector of
  the diagonal entries like.
}
\usage{
gsi.diagGenerate(x)
}
\arguments{
  \item{x}{ a vector}
}
\value{
  a diagonal matrix.
}
\details{
  The difference to original \code{\link[base]{diag}} is that it always
  gives a diagonal matrix and does nothing flawed in case of a length
  one vector.
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi.diagExtract}}, \code{\link[base]{diag}}
}
\examples{
diag(1:3)
#gsi.diagGenerate(1:3)
#gsi.diagGenerate(3)
diag(3)
}
\keyword{internal}
