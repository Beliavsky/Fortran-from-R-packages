\name{which.tile}
\alias{which.tile}
\title{
  Determine the tile containing a given point.
}
\description{
  Finds the Dirichlet/Voronoi tile, of a tessellation produced
  by \code{deldir()}, that contains a given point.
}
\usage{
which.tile(x, y, tl)
}
\arguments{
  \item{x}{
  The \code{x} coordinate of the point in question.
}
  \item{y}{
  The \code{y} coordinate of the point in question.
}
  \item{tl}{
  A tile list, as produced by the function \code{\link{tile.list}()}
  from a tessellation produced by \code{\link{deldir}()}.
}
}
\details{
  Just minimises the distance from the point in question to the
  points of the pattern determining the tiles.
}
\value{
  An integer equal to the index of the tile in which the given
  point lies.
}
\author{\rolf}

\seealso{
\code{\link{tile.list}()} \code{\link{deldir}()}.
}
\examples{
set.seed(42)
x   <- runif(20,0,100)
y   <- runif(20,0,100)
dxy <- deldir(x,y)
txy <- tile.list(dxy)
i   <- which.tile(30,50,txy) # The value of i here is 14.
plot(txy,showpoints=FALSE)
text(x,y,labels=1:length(txy),col="red")
points(30,50,pch=20,col="blue")
}
\keyword{ spatial }
