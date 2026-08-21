\name{gsi.simshape}
\alias{gsi.simshape}
\title{Internal function: Reshape an object to the shape type of another}
\description{
Reshape an object to the shape of another
}
\usage{
          gsi.simshape(x,oldx)
          }
\arguments{
  \item{x}{the data set to be reshaped}
  \item{oldx}{a data set of the intended shape}
}
\value{
  \code{x}, but re-shaped as \code{oldx}
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package. To use a gsi function, attach the \code{gsi}
environment.}
\details{
This function changes the dimension attribute accordingly.    
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
#gsi.simshape(matrix(1:4,nrow=1),1:3)
}
\keyword{internal}
