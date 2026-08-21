\name{gsiinternal2}
\alias{gsi.pairrelativeMatrix}
\title{Internal functions of the compositions package}
\description{
  Internal functions
}
\usage{
gsi.pairrelativeMatrix(names)
}
\arguments{
  \item{names}{a character vector provinding names}
}
\value{
  a matrix containing pairwise contrasts
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{gsi}}
}
\examples{
#gsi.pairrelativeMatrix(c("a","b","c"))
}
\keyword{internal}
