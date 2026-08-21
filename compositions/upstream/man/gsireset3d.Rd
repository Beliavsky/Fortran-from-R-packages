\name{gsireset3D}
\alias{gsi.reset3D}
\title{Internal functions of the compositions package}
\description{
  Internal function to reset the rgl device  
}
\usage{
gsi.reset3D(userMatrix=diag(rep(1,4)))
}
\arguments{
  \item{userMatrix}{the new user matrix of the rgl Device}
}
\value{
  none
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}
}
\keyword{internal}
