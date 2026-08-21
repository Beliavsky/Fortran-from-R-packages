\name{gsi.getD}
\alias{gsi.getD}
\alias{gsi.getN}
\title{Interal function: Get number of samples and number of parts in a compositional object}
\description{
 Get the number of samples N and the number of parts D of the
 compositions in an \code{\link{acomp}}, \code{\link{rcomp}},
   \code{\link{aplus}}, \code{\link{rplus}} object.
}
\usage{
gsi.getD(x)
gsi.getN(x)
          }
\arguments{
  \item{x}{an \code{\link{acomp}}, \code{\link{rcomp}},
   \code{\link{aplus}}, \code{\link{rplus}} object or something that
   could be cast to one of them}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\value{
  an integer giving the number of parts D or the number of samples N.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}

\examples{
gsi.getD(1:5)
#gsi.getN(1:5)
NCOL(1:5)
NROW(1:5)
data(SimulatedAmounts)
gsi.getD(sa.lognormals5)
#gsi.getN(sa.lognormals5)
}
\keyword{internal}
