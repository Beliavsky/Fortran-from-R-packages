\name{triang.list}
\alias{triang.list}
\title{
  Create a list of Delaunay triangles
}
\description{
  From an object of class \dQuote{deldir} produces a list
  of the Delaunay triangles in the triangulation of a set
  of points in the plane.
}
\usage{
triang.list(object)
}
\arguments{
  \item{object}{
  An object of class \dQuote{deldir} as produced by \code{\link{deldir}()}.
}
}
\value{
  A list each of whose components is a \eqn{3 \times 3}{3 x 3} or
  \eqn{3 \times 4}{3 x 4} data frame corresponding to one of the
  Delaunay triangles specified by \dQuote{object}.  The rows of each
  such data frame correspond to the vertices of the corresponding
  Delaunay triangle. The columns are:
  \itemize{
  \item \code{ptNum} (the index of the point in the original sequence
  of points that is being triangulated.  Note that if a point is
  one of a set of \emph{duplicated} points then \code{ptNum} is
  the \emph{first} of the indices of the points in this set.)
  \item \code{x} (the \eqn{x}-coordinate of the vertex)
  \item \code{y} (the \eqn{y}-coordinate of the vertex)
  \item \code{z} (the \dQuote{auxiliary value} or \dQuote{tag}
  \code{z} associated with the vertex; present only if such values
  were supplied in the call to \code{deldir()})
  }

  The returned value has an attribute \dQuote{rw} consisting of the
  enclosing rectangle of the triangulation.
}
\author{\rolf}
\section{Warning}{
  There may not actually \bold{be} any triangles determined by
  \code{object}, in which case this function returns an empty
  list with an \code{"rw"} attribute.  See \bold{Examples}. 
}
\note{
   The code of this function was taken more-or-less directly
   from code written by Adrian Baddeley for the \dQuote{delaunay()}
   function in the \dQuote{spatstat} package.
}
\seealso{ \code{\link{deldir}()}, \code{\link{plot.triang.list}()},
          \code{\link{tile.list}()}, \code{\link{plot.tile.list}()} }

\examples{
set.seed(42)
x  <- runif(20)
y  <- runif(20)
z  <- sample(1:100,20)
d1 <- deldir(x,y,z=z)
t1 <- triang.list(d1)
# A "triangulation" with no triangles!
d2 <- deldir(x=1:10,y=11:20)
plot(d2)
t2 <- triang.list(d2)
plot(t2,showrect=TRUE,rectcol="blue") # Pretty boring!
}
\keyword{spatial}
