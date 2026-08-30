\name{tilePerim}
\alias{tilePerim}
\title{
   Calculate tile perimeters.
}
\description{
   Calculates the perimeters of all of the Dirichlet (Voronoi)
   tiles in a tessellation of a set of planar points.  Also
   calculates the sum and the mean of these perimeters.
}
\usage{
tilePerim(object,inclbdry=TRUE)
}
\arguments{
  \item{object}{
  An object of class \code{tile.list} (as produced by
  \code{\link{tile.list}()} specifying the Dirichlet (Voronoi)
  tiles in a tessellation of a set of planar points.
}
  \item{inclbdry}{
  Logical scalar.  Should boundary segments (edges of tiles
  at least one of whose endpoints lies on the enclosing
  rectangle \code{rw} (see \code{\link{deldir}()}) be included
  in the perimeter?
}
}
\value{
   A list with components
   \item{perimeters}{
   A vector consisting of the values of the perimeters of
   the Dirichlet tiles in the tessellation.
   }
   \item{totalPerim}{
   The sum of \code{perimeters}.
   }
   \item{meanPerim}{
   The mean of \code{perimeters}.
   }
   \item{perComps}{
   A list whose entries are vectors consisting of the \dQuote{components} of the
   perimeters of each tile.  If/when the tiles are clipped, some
   tiles may be subdivided by the clipping into discontiguous parts.  The
   components referred to above are the perimeters of this parts.  If no
   subdivision has occurred then the vector in question has a single entry
   equal to the perimeter of the corresponding tile.  If subdivision has
   occurred then the perimeter of the tile is the sum of the perimeters
   of the components.
}
}
\author{\rolf}
\note{
   Function added at the request of Haozhe Zhang.
}

\seealso{
  \code{\link{tile.list}()}, \code{\link{plot.tile.list}()}
}
\examples{
    x <- runif(20)
    y <- runif(20)
    z <- deldir(x,y,rw=c(0,1,0,1))
    w <- tile.list(z)
    p1 <- tilePerim(w)
    p0 <- tilePerim(w,inclbdry=FALSE)
    p1$totalPerim - p0$totalPerim # Get 4 = the perimeter of rw.
    ss <- apply(as.matrix(z$dirsgs[,1:4]),1,
                function(x){(x[1]-x[3])^2 + (x[2]-x[4])^2})
    2*sum(sqrt(ss)) - p0$totalPerim # Get 0; in tilePerim() each interior
                                    # edge is counted twice.
    if(require(polyclip)) {
        CP <- list(x=c(0.49,0.35,0.15,0.20,0.35,0.42,
                       0.43,0.62,0.46,0.63,0.82,0.79),
                   y=c(0.78,0.86,0.79,0.54,0.58,0.70,
                       0.51,0.46,0.31,0.20,0.37,0.54))
        wc <- tile.list(z,clipp=CP)
        p2 <- tilePerim(wc) # Doesn't matter here if inclbdry is TRUE or FALSE.
        p2$perComps[["pt.6"]] # The tile for point 6 has got subdivided into
                              # two parts, a tetrahedron and a triangle.
        cul <- rainbow(10)[c(1,7,3:6,2,8:10)] # Rearranging colours to improve
                                              # the contrast between contiguous tiles.
        plot(wc,labelPts=TRUE,fillcol=cul)
    }
}
\keyword{spatial}
