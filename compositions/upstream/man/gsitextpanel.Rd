\name{gsi.textpanel}
\alias{gsi.textpanel}
\title{Internal function: A panel displaying a label only}
\description{
  A function useful as a text.panel in pairs.
}
\usage{
gsi.textpanel(x,y,lab,\dots)
}
\arguments{
  \item{x}{discarded}
  \item{y}{discarded}
  \item{lab}{text to be plotted to the middle of the panel}
  \item{\dots}{further graphical parameters passed to \code{\link[graphics]{text}}}
}
\details{
  The function is used against log-scale problems in pairs
  called by function \code{\link{boxplot.acomp}}.
  }
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}
}
\examples{
#data(iris)
#pairs(iris,text.panel=gsi.textpanel)
}
\keyword{internal}
