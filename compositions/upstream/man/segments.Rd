\name{segments}
\alias{segments}
\alias{segments.default}
\alias{segments.rmult}
\alias{segments.acomp}
\alias{segments.rcomp}
\alias{segments.aplus}
\alias{segments.rplus}
\title{Draws straight lines from point to point.}
\description{
  The function draws lines from a points x0 to a point y1 in the given geometry.
}
\usage{
          segments(x0,...)
          \method{segments}{default}(x0,...)
          \method{segments}{acomp}(x0,y1,...,steps=30,aspanel=FALSE)
          \method{segments}{rcomp}(x0,y1,...,steps=30,aspanel=FALSE)
          \method{segments}{aplus}(x0,y1,...,steps=30,aspanel=FALSE)
          \method{segments}{rplus}(x0,y1,...,steps=30,aspanel=FALSE)
          \method{segments}{rmult}(x0,y1,...,steps=30,aspanel=FALSE)
          }
\arguments{
  \item{x0}{dataset of points (of the given type) to draw the line from}
  \item{y1}{dataset of points (of the given type) to draw the line to}
  \item{\dots}{further graphical parameters}
  \item{steps}{the number of discretisation points to draw the segments,
    since the representation might not visually be a straight line}
  \item{aspanel}{Logical, indicates use as slave to do acutal drawing only.}

}
\details{

The default `segments.default(x0,...)` redirects to `segments` in package "graphics".

The other methods add lines to the graphics generated with the corresponding
plot functions of "compositions"-classes.
\cr

Adding to multipaneled plots redraws the plot completely, and is only
possible when the plot has been created with the plotting routines from
this library.

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}},\code{\link{lines.acomp}} 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
segments.acomp(acomp(c(1,2,3)),acomp(c(2,3,1)),col="red")
segments.rcomp(acomp(c(1,2,3)),acomp(c(2,3,1)),col="blue")

plot(aplus(sa.lognormals[,1:2]))
segments.aplus(aplus(c(10,20)),aplus(c(20,10)),col="red")
segments.rplus(rplus(c(10,20)),rplus(c(20,10)),col="blue")

plot(rplus(sa.lognormals[,1:2]))
segments.aplus(aplus(c(10,20)),aplus(c(20,10)),col="red")
segments.rplus(rplus(c(10,20)),rplus(c(20,10)),col="blue")




}
\keyword{multivariate}
