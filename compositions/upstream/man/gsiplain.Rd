\name{gsi.plain}
\alias{gsi.plain}
\title{Internal function: Convert to plain vector or matrix}
\description{
  The dataset is converted into a plain vector or matrix: data.frames are
  converted to data matrices and class attributes are removed.
}
\usage{
gsi.plain( x )
}
\arguments{
  \item{x}{The dataset to be converted}
}
\value{
  unclassed object, typically a vector or matrix.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
#gsi.plain(acomp(c(12,3,4)))
}
\keyword{internal}
