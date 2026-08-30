\name{print.tileInfo}
\alias{print.tileInfo}
\title{
   Print a summary of tile information.
}
\description{
   Print a reasonably readable summary of an object of class
   \code{tileInfo} as produced by the \code{tileInfo()} function.
}
\usage{
\method{print}{tileInfo}(x, digits = 4, npl = 6, ...)
}
\arguments{
  \item{x}{An object of class \code{tileInfo} as produced by the
  \code{tileInfo()} function.
}
  \item{digits}{
  Integer scalar.  The (maximum) number of decimal digits to which
  the output is to be printed.
}
  \item{npl}{
  Integer scalar.  \dQuote{Number per line}. It specifies the
  (maximum) number of values per line.  Used (only) when printing
  the edge lengths component of \code{x}.  It effects a tidy
  \dQuote{folding} of the printed vector of edge lengths of the
  tile associated with a given point.  If you increase the value of
  \code{digits} you may wish to decrease the value of \code{npl},
  and vice versa.
}
  \item{\dots}{
  Not used.  Present for compatibility with the generic
  \code{print()} function.
}
}
\details{
  The list produced by \code{tileInfo()} is a bit messy and
  hard to comprehend, especially if there is a large number
  of tiles.  This print method produces a screen display which
  is somewhat more perspicuous.
 
  There are four components to the display:

  \itemize{
  \item A list, each entry of which is the vector of edge lengths
  of the tile.  Each edge length is formatted to have a number
  of digits specified by the \code{digits} argument.  Each list
  entry may be displayed over a number of lines. The first of these
  lines is prefixed by an \dQuote{informative} string indicating
  the point that determines the tile whose edge lengths are being
  printed.  The string is formed from the identifier of the point.
  See \code{\link{deldir}()}, \code{\link{plot.deldir}()} and
  \code{\link{getNbrs}()}.  The identifier may consist essentially
  of the index of the point in the sequence of points that is
  being tessellated.

  Succeeding lines, corresponding to the same list entry, are
  prefixed with a number of blanks so as to produce an aesthetically
  pleasing alignment.

  \item A table of the edge counts of the tiles.

  \item A simple print out of the areas of the tiles (rounded
  to a maximum of \code{digits} decimal digits).

  \item A simple print out of the perimeters of the tiles (rounded
  to a maximum of \code{digits} decimal digits).
  }

  This screen display is for \dQuote{looking at} only.  In order
  to do further calculations on the output of \code{tileInfo} it
  is necessary to delve into the bowels of \code{x} and extract
  the relevant bits.
}
\value{
  None.
}
\author{\rolf}

\seealso{
  \code{\link{tileInfo}()}
}
\examples{
set.seed(179)
x    <- runif(100)
y    <- runif(100)
dxy  <- deldir(x,y,rw=c(0,1,0,1))
ixy1 <- tileInfo(dxy)
print(ixy1)
ixy2 <- tileInfo(dxy,bndry=TRUE)
print(ixy2)
if(require(polyclip)) {
    CP <- list(x=c(0.49,0.35,0.15,0.20,0.35,0.42,
                   0.43,0.62,0.46,0.63,0.82,0.79),
               y=c(0.78,0.86,0.79,0.54,0.58,0.70,
                   0.51,0.46,0.31,0.20,0.37,0.54))
    ixy3 <- tileInfo(dxy,clipp=CP)
    options(width=120) # And enlarge the console window.
    print(ixy3) # 33 tiles are retained.
    print(ixy3$perimeters$perComps) # The tiles for points 9 and 94 have
                                    # been split into two components.
}
}
\keyword{ utilities }
