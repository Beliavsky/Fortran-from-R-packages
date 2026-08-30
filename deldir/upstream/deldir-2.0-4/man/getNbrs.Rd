\name{getNbrs}
\alias{getNbrs}
\title{
   List Delaunay neighbours.
}
\description{
   Lists the indices (or identifiers if these are provided) of the
   Delaunay neighbours of each point in the set of points being
   triangulated/tessellated.
}
\usage{
getNbrs(object, interior = NULL)
}
\arguments{
  \item{object}{
  An object of class \code{"deldir"}, as returned by the function
  \code{deldir()}.
}
  \item{interior}{
  Either a rectangle, i.e. a numeric vector of length 4,
  \code{(xmin,xmax,ymin,ymax)}, or a polygon (a list with components
  named \code{"x"} and \code{"y"}, which are numeric vectors
  of equal length, say \code{n}.  The pairs \code{(x[i],y[i])}
  specify the vertices of the polygon and should be presented in
  \emph{anticlockwise} order.  The last vertex should \bold{NOT}
  repeat the first.

  If \code{interior} is given as a rectangle, its entries must
  satisfy \code{xmin < xmax} and \code{ymin < ymax}.  Moreover,
  in this case, \code{interior} must be a subset of the rectangular
  window \code{rw} (see \code{\link{deldir}()}) with respect to
  which the points in question are being triangulated/tessellated.
  That is, its entries must satisfy \code{rw[1] <= xmin}, \code{xmax <= rw[2]},
  \code{rw[3] <= ymin} and \code{ymax <= rw[4]}.
}
}
\details{
  If \code{interior} is specified then Delaunay neighbours are listed
  only for those points which lie in \code{interior}.  Note that
  these \emph{neighbours} need not themselves lie in \code{interior}.

  Note also that it is possible for points \code{i} and \code{j} to
  be neighbours even though the \dQuote{clipped} versions of the tiles
  are discontiguous.
}

\value{
  A (named) list of length equal to the number points in the set
  being triangulated/tessellated, or to the number of such points that
  lie in \code{interior} if that argument was specified.  The names
  of the list are the identifiers of the points as specified by
  \code{id} if \code{id} was specified in the call to \code{deldir()}
  that produced \code{object}.  If \code{id} was not specified, then
  the names are the indices of the points, coerced to character mode.

  The entries of this list are vectors of the identifiers or
  indices of the Delaunay neighbours of the point corresponding to
  that entry.
}
\references{
  See \code{\link{deldir}()} for references.
}
\author{\rolf}
\section{Warning}{
Be careful about addressing the entries of the list returned
by this function.  If \code{id} was not specified in the call
to \code{deldir()} that produced \code{object}, then the names
of this list are the point indices coerced to character mode.
If \code{interior} was specified then the name of \eqn{i}-th entry
of the list will not in general be \code{i}.  E.g. given that
\code{xxx} is the list returned by this function, \code{xxx[[14]]}
will \bold{not} in general give the Delaunay neighbours of point 14.
Instead, specify \code{xxx[["14"]]} or \code{xxx[[id[14]]]} where
\code{id} is the vector of identifiers supplied in the call to
\code{deldir()}.
}

\examples{
set.seed(42)
x <- runif(60)
y <- runif(60)
dxy <- deldir(x,y,rw=c(0,1,0,1))
nbrs <- getNbrs(dxy,interior=c(0.2,0.8,0.2,0.8))
nbrs[["14"]] # Correct.
nbrs[[14]]   # Incorrect.
names(nbrs)[14] # The result is 42.
# Thus nbrs[[14]] actually gives the Delaunay neighbours of point 42.

# Demonstrate that neighbours can have discontiguous clipped tiles.
if(require(polyclip)) {
    x <- c(0.38,0.44,0.04,0.97,0.43,0.96,0.89,0.64,0.97,0.62,0.33,0.35,
           0.40,0.78,0.04,0.75,0.68,0.17,0.26,0.51)
    y <- c(0.68,0.98,0.76,0.57,0.85,0.19,0.27,0.83,0.69,0.24,0.04,0.14,
           0.22,0.48,0.20,0.72,0.01,0.38,0.51,0.00)
    CP <- list(x=c(0.72,0.93,0.76,0.61,-0.03,-0.04,0.41),
               y=c(0.46,0.76,0.94,1.03,1.01,0.37,0.31))
    dxy <- deldir(x,y,rw=c(0,1,0,1))
    TL <- tile.list(dxy)
    plot(TL,labelPts=TRUE)
    plot(TL[16],clipp=CP,fillcol="orange",labelPts=TRUE,add=TRUE)
    polygon(CP,border="red")
    nbrs <- getNbrs(dxy,interior=CP) # Tiles are clipped to CP.
# Note that point 14 is a neighbour of point 16, even though their
# clipped tiles do not meet.
}
}
\keyword{ spatial }
