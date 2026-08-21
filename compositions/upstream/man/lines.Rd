\name{lines}
\alias{lines.rmult}
\alias{lines.acomp}
\alias{lines.rcomp}
\alias{lines.aplus}
\alias{lines.rplus}
\title{Draws connected lines from point to point.}
\description{
  Functions taking coordinates given in various ways and
     joining the corresponding points with line segments.
}
\usage{
          \method{lines}{acomp}(x,...,steps=30,aspanel=FALSE)
          \method{lines}{rcomp}(x,...,steps=30,aspanel=FALSE)
          \method{lines}{aplus}(x,...,steps=30,aspanel=FALSE)
          \method{lines}{rplus}(x,...,steps=30,aspanel=FALSE)
          \method{lines}{rmult}(x,...,steps=30,aspanel=FALSE)
          }
\arguments{
  \item{x}{a dataset of the given type}
  \item{\dots}{further graphical parameters}
  \item{steps}{the number of discretisation points to draw the segments,
    which might be not visually straight.}
  \item{aspanel}{Logical, indicates use as slave to do acutal drawing only.}
}
\details{
The functions add lines to the graphics generated with the corresponding
plot functions.
\cr
Adding to multipaneled plots, redraws the plot completely and is only
possible, when the plot has been created with the plotting routines from
this library.
\cr
For the rcomp/rplus geometries the main problem is providing a function
that reasonably works with lines leaving the area. We tried to use a
policy of cuting the line at the actual borders of the (high
dimensional) simplex. That can lead to very strange visual impression
showing lines ending somewhere in the middle of the plot. However these
lines actually hit some border of the simplex that is not shown in the
plot. A hyper dimensional tetrahedron is even more difficult to imagin
than a hyperdimensional cube. 

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{plot.acomp}}, \code{\link{straight}} 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
lines(acomp(sa.lognormals),col="red")
lines(rcomp(sa.lognormals),col="blue")

plot(aplus(sa.lognormals[,1:2]))
lines(aplus(sa.lognormals[,1:2]),col="red")
lines(rplus(sa.lognormals)[,1:2],col="blue")

plot(rplus(sa.lognormals[,1:2]))
tt<-aplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="blue")
tt<-rmult(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="green")

}
\keyword{multivariate}
