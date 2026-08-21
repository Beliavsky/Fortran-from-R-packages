\name{gsiCoorInfo}
\alias{gsi.coorInfo}
\alias{gsi.setCoorInfo}
\alias{gsi.addCoorInfo}
\alias{gsi.getCoorInfo}
\title{Internal functions of the compositions package}
\description{
  Sets or gets the triangular/barycentric coordinate system for
  the current ternary diagrams in the current plot.
}
\usage{
gsi.getCoorInfo()
gsi.setCoorInfo(...,all=list(...))
gsi.addCoorInfo(...,all=list(...))
gsi.plots
}
\arguments{
  \item{\dots}{The info to be set. See details.}
  \item{all}{The complete list of information}
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\details{
  Typically the information set contains a center \code{mean} (an acomp
  object),
  a scaling \code{scale} (a real number),
  the geometry of rescaling \code{geo} (typically "acomp" or "rcomp"),
  the dimension \code{D}, the margin \code{margin} type ("acomp" or
  "rcomp" or a variable number or name) and eventually a numeric vector
  of the two informations \code{d}.

}

\seealso{
  \code{\link{gsi}}, 
}
\keyword{internal}
