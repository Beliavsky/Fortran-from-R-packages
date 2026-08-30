\name{print.deldir}
\alias{print.deldir}
\title{
   Print some information about a tessellation/triangulation.
}
\description{
   Prints a very brief description of an object of class \code{"deldir"}
   as returned by \code{\link{deldir}()}.
}
\usage{
\method{print}{deldir}(x,digits=NULL,...)
}
\arguments{
  \item{x}{
  A Delaunay triangulation and Dirichlet (Voronoi) tessellation
  of a set of points (object of class \code{"deldir"}).
}
  \item{digits}{
  Integer scalar.  The number of digits to which to round the
  numeric information before printing.  Note this may be give
  negative values.  (See \code{\link{round}()}.)
}
   \item{\dots}{
   Not used.
}
}
\details{
   This is a method for the generic \code{\link{print}()} function.
}
\author{\rolf}
\seealso{
   \code{\link{print}()}
}
\examples{
set.seed(42)
x    <- rnorm(200,0,4)
y    <- rnorm(200,0,4)
dxy1 <- deldir(x,y)
dxy2 <- deldir(x,y,rw=c(-12,12,-11,11))
dxy1
dxy2
print(dxy1,digits=4)
}
\keyword{ utilities }
