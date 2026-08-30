\name{triMat}
\alias{triMat}
\title{
  Produce matrix of triangle vertex indices.
}
\description{
  Lists the indices of the vertices of each Delaunay triangle in
  the triangulation of a planar point set.  The indices are listed
  (in increasing numeric order) as the rows of an \eqn{n \times 3}{n x 3}
  matrix where \eqn{n} is the number of Delaunay triangles in the
  triangulation.
}

\usage{
triMat(object)
}
\arguments{
  \item{object}{
  An object of class \code{deldir} (as produced by the function
  \code{\link{deldir}()}) specifying the Delaunay triangulation
  and Dirichlet (Voronoi) tessellation of a planar point set.
}
}
\details{
  This function was suggested by Robin Hankin of the School of
  Mathematical and Computing Sciences at Auckland University of
  Technology.
}
\value{
  An \eqn{n \times 3}{n x 3} matrix where \eqn{n} is the number
  of Delaunay triangles in the triangulation specified by \code{object}.
  The \eqn{i^{th}}{i-th} row consists of the indices (in the original
  list of points being triangulated) of vertices of the \eqn{i^{th}}{i-th}
  Delaunay triangle.  The indices are listed in increasing numeric
  order in each row.
}
\section{Note}{
  Earlier versions of this function (prior to release 0.1-14
  of \bold{deldir}) could sometimes give incorrect results.
  This happened if the union of three contiguous Delaunay triangles
  happened to constitute another triangle.  This latter triangle
  would appear in the list of triangles produced by \code{triMat()}
  but is \emph{not} itself a Delaunay triangle.  The updated version
  of \code{triMat()} now checks for this possibility and gives
  (\emph{I think!}) correct results.

  Many thanks to Jay Call, who pointed out this bug to me.
}

\author{\rolf}
\seealso{ \code{\link{deldir}()}
          \code{\link{triang.list}()}
	  \code{\link{plot.triang.list}()}
}
\examples{
# These are the data used by Jay Call to illustrate the bug
# that appeared in a previous incarnation of triMat.
xy <- data.frame(
          x = c(0.048,0.412,0.174,0.472,0.607,0.565,0.005,0.237,0.810,0.023),
          y = c(0.512,0.928,0.955,0.739,0.946,0.134,0.468,0.965,0.631,0.782)
      )
dxy <- deldir(xy)
M   <- triMat(dxy)
plot(dxy,wlines="triang",num=TRUE,axes=FALSE,cmpnt_col=c(1,1,1,1,2,1))
# The triangle with vertices {4,5,8} was listed in the output of
# the previous (buggy) version of triMat().  It is NOT a Delaunay
# triangle and hence should NOT be listed.
}
\keyword{spatial}
