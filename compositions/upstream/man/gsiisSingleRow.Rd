\name{gsi.isSingleRow}
\alias{gsi.isSingleRow}
\title{Internal function: Can something be considered as a single
  multivariate data item?}
\description{
  Checks wether something can be regarded as a single multivariate item,
  being a matrix or a vector, which is only a row or a column.
}
\usage{
gsi.isSingleRow(X)
}
\arguments{
  \item{X}{the matrix or vector to be checked}
}
\value{
  a logical value 
}
\details{
  It is defined as \code{
    gsi.isSingleRow <- function(X) {
      return( NROW(X) == 1 || NCOL(X) ==1 )
    }
  }
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}
}
\examples{
#gsi.isSingleRow(1:10)
}
\keyword{internal}
