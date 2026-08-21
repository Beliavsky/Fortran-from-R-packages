\name{gsi.eq}
\alias{gsi.eq}
\alias{gsi.eq}
\title{Internal function: Checking equality of IEEE special numbers}
\description{
  Works like == 
}
\usage{
gsi.eq(x,y)
}
\arguments{
  \item{x}{A numeric vector or matrix}
  \item{y}{A single numerical value, possibly including NaN,Inf,NA,-Inf,-0,0}
}
\value{
  A logical vector of the shape of x with TRUE exactly where the given
  value shows up.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\details{
  unlike \code{==} does not return NA's and all sorts of NA's as an
  actuall value.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{compositions.missing}}
}
\examples{
#gsi.plain(acomp(c(12,3,4)))
}
\keyword{internal}
