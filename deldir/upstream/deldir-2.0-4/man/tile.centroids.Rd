\name{tile.centroids}
\alias{tile.centroids}
\title{
Compute centroids of Dirichlet (Voronoi) tiles
}
\description{
Given a list of Dirichlet tiles, as produced by \code{tile.list()},
produces a data frame consisting of the centroids of those tiles.
}
\usage{
tile.centroids(tl)
}
\arguments{
  \item{tl}{A list of the tiles (produced by \code{tile.list()})
  in a Dirichlet tessellation of a set of planar points.}
}
\value{
  A data frame with two columns named \code{x} and \code{y}.
  Each row of this data frame constitutes the centroid of one
  of the Dirichlet tiles.
}
\references{
URL http://en.wikipedia.org/wiki/Centroid
}
\author{\rolf}

\seealso{ \code{\link{tile.list}()} }
\examples{
set.seed(42)
x <- runif(20)
y <- runif(20)
d <- deldir(x,y)
l <- tile.list(d)
g <- tile.centroids(l)
plot(l,close=TRUE)
points(g,pch=20,col="red")
}
\keyword{ spatial }
