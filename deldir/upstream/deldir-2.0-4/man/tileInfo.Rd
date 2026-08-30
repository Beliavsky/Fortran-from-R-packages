\name{tileInfo}
\alias{tileInfo}
\title{
   Extract information from a tile list.
}
\description{
   Produces a summary of information about the tiles in an
   object of class \code{deldir} as produced by the function
   \code{deldir()}.
}
\usage{
tileInfo(object, bndry = FALSE, clipp=NULL)
}
\arguments{
  \item{object}{ An object of class \code{deldir} as produced
  by the function \code{deldir()}.
}
  \item{bndry}{
  Logical scalar.  If \code{TRUE} then the \dQuote{boundary} tiles
  (those tiles having edges forming part of the \dQuote{rectangular
  window} enclosing the tessellation) are included in the summary.
  Otherwise they are excluded.  Note that if all tiles are
  boundary tiles and \code{bndry} is \code{FALSE}, then the list of
  included tiles is empty.  In this setting the information summary
  would be vacuous and hence an error is thrown.

}
  \item{clipp}{An object specifying a polygon to which the tiles of
  the tessellation should be clipped.  See \code{link{tile.list}()}
  for more information.  Note that tiles which were originally
  boundary tiles may no longer be boundary tiles after being clipped.
}
}
\section{Remark}{
There is a \code{print()} method for class \code{"tileInfo"} which
produces a convenient display of the information returned by this
function.
}
\section{Acknowledgement}{
Thanks to Krisztina Konya of Ruhr-University Bochum, who provided
an example illustrating the need for an error trap in the setting
in which all tiles are boundary tiles.
}
\value{An object of class \code{"tileInfo"} which consists of a
list with components:
  \item{indivTiles}{This is itself a list. If \code{clipp}
  is \code{NULL} then this list has one entry for each tile in
  \dQuote{object}.  If \code{clipp} is not \code{NULL} then tiles
  are retained only if they have non-empty intersection with the
  polygon specified by \code{clipp}.  The list \code{indivTiles} is
  in fact a \emph{named} list, the names being of form \code{pt.n},
  where \code{n} is equal to the value of \code{ptNum} (see below)
  corresponding to the tile.

  The entries of \code{indivTiles} are themselves in turn lists with
  entries
  \itemize{
  \item \code{edges}:  a matrix whose rows consists of the \code{x}
  and \code{y} coordinates of the endpoints of each edge of the tile
  \item \code{edgeLengths}: a vector of the lengths of the edges
  of the tile
  \item \code{area}: a positive number equal to the area of the tile
  \item \code{ptNum} an integer equal to the number of the point
  determining the tile.  Note that this is the number of the point
  in the \emph{original} sequence of points that were tessellated.
  }
  }
  \item{allEdgeCounts}{An integer vector of the number of
  edges for each of the tiles.}
  \item{tabEdgeCounts}{A table of \code{allEdgeCounts}.}
  \item{allEdgeLengths}{A vector of all of the tile edge lengths;
  a catenation of the \code{edgeLengths} components of the entries
  of \code{indivTiles}.  Note that there will be many duplicate
  lengths since each tile edge is, in general, an edge of \emph{two}
  tiles.}
  \item{Areas}{A vector of the areas of the tiles.}
  \item{uniqueEdgeLengths}{A vector of the lengths of the tiles edges
  with the duplicates (which occur in \code{allEdgeLengths}) being
  eliminated. Each tile edge is represented only once.}
  \item{perimeters}{A list, as returned by \code{\link{tilePerim}()}
  containing the perimeters of the tiles, as well as the total and
  the mean of these perimeters.  In addition \code{perimeters} has
  a component \code{perComps} giving the breakdown of the perimeters
  into the perimeters of the parts into which tiles may have been
  subdivided by the clipping process.}
}
\author{\rolf}

\seealso{
 \code{\link{deldir}()} \code{\link{tile.list}()} \code{\link{print.tileInfo}()}
 \code{\link{tilePerim}()}
}
\examples{
set.seed(42)
x    <- runif(20)
y    <- runif(20)
dxy  <- deldir(x,y,rw=c(0,1,0,1))
ixy1 <- tileInfo(dxy)
ixy2 <- tileInfo(dxy,bndry=TRUE)
if(require(polyclip)) {
    CP <- list(x=c(0.49,0.35,0.15,0.20,0.35,0.42,
                   0.43,0.62,0.46,0.63,0.82,0.79),
               y=c(0.78,0.86,0.79,0.54,0.58,0.70,
                   0.51,0.46,0.31,0.20,0.37,0.54))
    ixy3 <- tileInfo(dxy,clipp=CP) # 10 tiles are retained; the third tile,
                                   # corresponding to point 6, is
                                   # subdivided into two parts.
}
}
\keyword{ spatial }
