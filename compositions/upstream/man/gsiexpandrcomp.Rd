\name{gsi.expandrcomp}
\alias{gsi.expandrcomp}
\title{Internal function: Scaling rcomp}
\description{
  This functions tries to compute something similar to a scaling of an
  acomp object in the context of the rcomp geometry.
}
\usage{
gsi.expandrcomp(x,alpha)
          }
\arguments{
  \item{x}{an rcomp object}
  \item{alpha}{a number or a numeric vector between 0 and 1}
}
\value{
  an rcomp-object
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\note{It is better not to use gsi.* functions directly since they are internal
  functions of the package
}  
\examples{
#gsi.expandrcomp(rcomp(1:3),0.5)
}
\keyword{internal}
