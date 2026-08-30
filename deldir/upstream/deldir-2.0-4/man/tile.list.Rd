\name{tile.list}
\alias{tile.list}
\title{ Create a list of tiles in a tessellation }
\description{
For each point in the set being tessellated produces a list
entry describing the Dirichlet/Voronoi tile containing that
point.
}
\usage{ tile.list(object,minEdgeLength=NULL,clipp=NULL) }
\arguments{
  \item{object}{ An object of class \code{deldir} as produced
   by the function \code{\link{deldir}()}.}
  \item{minEdgeLength}{
  Positive numeric scalar specifying the minimum length that
  an edge of a tile may have.  It is used to eliminate edges
  that are effectively of zero length, which can cause tiles
  to be \dQuote{invalid}.  This argument defaults to
  \code{sqrt(.Machine$double.eps)} time the diameter (length
  of the diagonal) of the \dQuote{rectangular window}
  associated with the tessellation.  This rectangular window
  is available as the \code{rw} component of \code{object}.
}
  \item{clipp}{An object specifying a polygon to which the
  tessellation, whose tiles are being determined, should be clipped.
  It should consist either of:
  \itemize{
    \item a list containing two components x and y giving the
    coordinates of the vertices of a single polygon. The last vertex
    should not repeat the first vertex. Or:
    \item a list of list(x,y) structures giving the coordinates of
    the vertices of several polygons.
  }
  If this argument is provided then the tiles in the list that
  is produced are \dQuote{clipped} to the polygon specified by
  \code{clipp}.  Empty tiles (those which do not intersect the
  polygon specified by \code{clipp}) are omitted.  The clipping
  process may subdivide some of the tiles into two or more
  discontiguous parts.
}
}
\value{
A list with one entry for each of the points in the set being
tessellated, or for each of the tiles that are retained after clipping
if \code{clipp} is not \code{NULL}.  Each entry is in turn a list
with a number of components.  These components always include:

\item{ptNum}{The index of the point in the original sequence of points
that is being tessellated.  Note that if a point is one of a set
of \emph{duplicated} points then \code{ptNum} is the \emph{first} of the
indices of the points in this set.}
\item{pt}{The coordinates of the point whose tile is being described.}
\item{area}{The area of the tile.}

If the tile in question has \emph{not} been subdivided by the
clipping process then the list components also include:

\item{x}{The \code{x} coordinates of the vertices of the tile, in
anticlockwise order.}
\item{y}{The \code{y} coordinates of the vertices of the tile, in
anticlockwise order.}
\item{bp}{Vector of logicals indicating whether the tile vertex is a
``real'' vertex, or a \emph{boundary point}, i.e. a point where the
tile edge intersects the boundary of the enclosing rectangle.}

If the tile in question \emph{has} been subdivided then the list
does not have the foregoing three components but rather has a
component \code{tileParts} which is in turn a list of length equal
to the number of parts into which the tile was subdivided.  Each
component of \code{tileParts} is yet another list with four
components \code{x}, \code{y}, \code{bp} and \code{area} as described above
and as are appropriate for the part in question. 

\item{z}{The \dQuote{auxiliary value} or \dQuote{tag} associated
with the \code{pt}; present only if such values were supplied in
the call to \code{deldir()}.}
}
\author{\rolf}
\section{Acknowledgement}{
The author expresses sincere thanks to Majid Yazdani who found and
pointed out a serious bug in \code{tile.list} in a previous version
(0.0-5) of the \code{deldir} package.
}

\section{Warning}{
The set of vertices of each tile may be ``incomplete''.  Only
vertices which lie within the enclosing rectangle, and ``boundary
points'' are listed.

Note that the enclosing rectangle may be specified by the user
in the call to \code{\link{deldir}()}.

In contrast to some earlier versions of \code{deldir}, the corners
of the enclosing rectangle are now included as vertices of tiles.
I.e. a tile which in fact extends beyond the rectangular window
and contains a corner of that window will have that corner added
to its list of vertices.  Thus the corresponding polygon is the
intersection of the tile with the enclosing rectangle.
}

\seealso{ \code{\link{deldir}()}, \code{\link{plot.tile.list}()}
          \code{\link{triang.list}()} \code{\link{plot.triang.list}()} }

\examples{
    set.seed(42)
    x <- runif(20)
    y <- runif(20)
    z <- deldir(x,y)
    w <- tile.list(z)

    z <- deldir(x,y,rw=c(0,1,0,1))
    w <- tile.list(z)

    z <- deldir(x,y,rw=c(0,1,0,1),dpl=list(ndx=2,ndy=2))
    w <- tile.list(z)
    if(require(polyclip)) {
        CP <- list(x=c(0.49,0.35,0.15,0.20,0.35,0.42,
                       0.43,0.62,0.46,0.63,0.82,0.79),
                   y=c(0.78,0.86,0.79,0.54,0.58,0.70,
                       0.51,0.46,0.31,0.20,0.37,0.54)) 
        wc <- tile.list(z,clipp=CP) # 10 tiles are retained; the third tile,
                                    # corresponding to point 6, is
                                    # subdivided into two parts.
    # Determine the tiles on the border of a clipping region.
    # Example due to Huimin Wang.
    set.seed(112)
    x      <- runif(100)
    y      <- runif(100)
    dxy    <- deldir(x,y)
    txy    <- tile.list(dxy)
    chind  <- chull(x,y)
    bdry   <- list(x=x[chind],y=y[chind])
    ctxy   <- tile.list(dxy,clipp=bdry)
    inside <- lapply(ctxy,function(tile,bdry){insidePoly(tile$x,tile$y,bdry)},
                                            bdry=bdry)
    border <- sapply(inside,function(x){any(!x) | any(attr(x,"on.boundary"))})
    plot(ctxy[border],main="Border tiles")
    }
}
\keyword{spatial}
