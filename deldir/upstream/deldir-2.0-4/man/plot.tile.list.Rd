\name{plot.tile.list}
\alias{plot.tile.list}
\title{ Plot Dirichlet (Voronoi) tiles }
\description{ A method for \code{plot}.  Plots (sequentially)
the tiles associated with each point in the set being tessellated. }
\usage{
\method{plot}{tile.list}(x, verbose = FALSE, close = FALSE, pch = 1,
                         fillcol = getCol(x,warn=warn), col.pts=NULL,
                         col.lbls=NULL,border=NULL, showpoints = !labelPts,
                         add = FALSE, asp = 1, clipp=NULL, xlab = "x",
                         ylab = "y", main = "", axes=TRUE,warn=TRUE,
                         labelPts=FALSE,adj=NULL,...) 
}
\arguments{
  \item{x}{
  A list of the tiles in a tessellation, as produced
  the function \code{\link{tile.list}()}.
}
  \item{verbose}{
  Logical scalar; if \code{TRUE} the tiles are
  plotted one at a time (with a \dQuote{Go?} prompt after each)
  so that the process can be watched.
}
  \item{close}{
  Logical scalar; if \code{TRUE} the outer edges of of the tiles
  (i.e. the edges which are constitute parts of the edges of the
  enclosing rectangle) are drawn.  Otherwise tiles on the periphery
  of the tessellation are left \dQuote{open}.
}
  \item{pch}{
  The plotting character (or vector of plotting
  characters) with which to plot the points of the pattern which
  was tessellated.  Ignored if \code{showpoints} is \code{FALSE}.
}
  \item{fillcol}{
  Optional vector (possibly of length 1, i.e. a scalar) whose
  entries can be interpreted as colours by \code{\link{col2rgb}()}.
  The \eqn{i}-th entry indicates with which colour to fill the
  \eqn{i}-th tile.  Note that an \code{NA} entry cause the tile
  to be left unfilled.  This argument will be replicated to have
  length equal to the number of tiles.  The default value is
  created (using the tile \dQuote{tags}, i.e. the \code{z}
  components of the tiles, if these exist) by
  the undocumented function \code{getCol()}.  Note that if
  these \code{z} components are not present, then \code{getCol()}
  returns \code{NA}.
}
  \item{col.pts}{
  Optional vector like unto \code{fillcol} whose entries can be
  interpreted as colours by \code{\link{col2rgb}()}.  The \eqn{i}-th
  entry indicates with which colour to plot the \eqn{i}-th point.
  This argument will be replicated to have length equal to the
  number of tiles.  Ignored if \code{showpoints} is \code{FALSE}.
}
  \item{col.lbls}{
  Optional vector like unto \code{col.pts}. Determines the colours
  in which the point labels (see \code{labelPts} below) are plotted.
  This argument will be replicated to have length equal to the
  number of tiles.  Ignored if \code{labelPts} is \code{FALSE}.
}
  \item{border}{
  A scalar that can be interpreted as a colour by \code{col2rgb()},
  indicating the colour with which to plot the tile boundaries.
  Defaults to black unless all of the fill colours specified by
  \code{fillcol} are black, in which case it defaults to white.
  If \code{length(border) > 1} then only its first entry is used.
}
  \item{showpoints}{Logical scalar; if \code{TRUE} the points of
  the pattern which was triangulated/tessellated are plotted.
  If \code{clipp} is specified, then only those points
  associated with non-empty clipped tiles are shown.  Note that
  some such points will be external to the polygon \code{clipp}
  and that these \dQuote{external} points \emph{are} shown.
}
  \item{add}{
  Logical scalar; should the plot of the tiles be added
  to an existing plot?
}
  \item{asp}{The aspect ratio of the plot; integer scalar or
  \code{NA}.  Set this argument equal to \code{NA} to allow the data
  to determine the aspect ratio and hence to make the plot occupy the
  complete plotting region in both \code{x} and \code{y} directions.
  This is inadvisable; see the \bold{Warnings}.
}
  \item{clipp}{An object specifying a polygon to which the tessellation
  being plotted should be clipped.  It should consist either of:
  \itemize{
    \item a list containing two components x and y giving the
    coordinates of the vertices of a single polygon. The last vertex
    should not repeat the first vertex. Or:
    \item a list of list(x,y) structures giving the coordinates of
    the vertices of several polygons.
  }
  If this argument is provided then the plot of the tessellation
  is \dQuote{clipped} to the polygon specified by \code{clipp}.
}
  \item{xlab}{
   Label for the \code{x}-axis (used only if \code{add} is \code{FALSE}).
}
  \item{ylab}{
   Label for the \code{y}-axis (used only if \code{add} is \code{FALSE}).
}
  \item{main}{
   A title for the plot (used only if \code{add} is \code{FALSE}).
}
  \item{axes}{
   Logical scalar.  Should axes (with a marked numeric scale)
   be plotted?  This argument is used only if \code{add} is
   \code{FALSE}.  It allows the plotting of such axes to be
   suppressed.  The argument was added at the request of Damon
   Tutunjian.  To suppress the plotting of axis \emph{labels}
   set \code{xlab=""} and/or \code{ylab=""}.
}
  \item{warn}{Logical scalar passed to the internal function \code{getCol()}.
  Should a warning be issued if the \code{z} components of the entries
  of \code{x} cannot all be interpreted as colours?  (See \bold{Notes}.)
}
  \item{labelPts}{Logical scalar; if \code{TRUE} the labels of the
  points determining the tiles are plotted in the tiles.  Note
  that if \code{labelPts} is \code{TRUE} then \code{showpoints}
  defaults to \code{FALSE}
}
  \item{adj}{The \dQuote{adjustment} argument to \code{text()}.
  If \code{labelPts} and \code{showpoints} are both \code{TRUE} it
  defaults to \code{-1} (so that the labels and point symbols
  are not superimposed).  If \code{labelPts} is \code{TRUE} and
  \code{showpoints} is \code{FALSE} it defaults to \code{0}.
  If \code{labelPts} is \code{FALSE} it is ignored.
}
  \item{...}{Optional arguments; may be passed to \code{points()}
  and \code{text()}.
}
}
  \section{Warnings}{
  \itemize{
  \item As of release 1.0-8 the argument \code{number} of
  \code{\link{plot.deldir}()} and \code{\link{plot.tile.list}()}
  was changed to \code{labelPts}.  As a consequence the argument
  \code{col.num} of this function has had its name changed to
  \code{col.lbls}.

  \item The behaviour of this function with respect to
  \dQuote{clipping} has changed substantially since the previous
  release of \code{deldir}, i.e. 1.1-0.  The argument \code{clipwin}
  has been re-named \code{clipp} (\dQuote{p} for \dQuote{polygon}).
  Clipping is now effected via the new package \code{polyclip}.
  The \code{spatstat} package is no longer used.  The argument
  \code{use.gpclib} has been eliminated, since \code{gpclib} (which
  used to be called upon by \code{spatstat} has been superseded by
  \code{polyclip} which has an unrestrictive license.

  \item As of release 0.1-1 of the \code{deldir} package, the
  argument \code{fillcol} to this function \emph{replaces} the old
  argument \code{polycol}, but behaves somewhat differently.

  \item The argument \code{showrect} which was present in versions
  of this function prior to release 0.1-1 has been eliminated.
  It was redundant.

  \item As of release 0.1-1 the \code{col.pts} argument \emph{might}
  behave somewhat differently from how it behaved in the past.

  \item The arguments \code{border}, \code{clipp}, and \code{warn}
  are new as of release 0.1-1.

  \item Users, unless they \emph{really} understand what they are
  doing and why they are doing it, are \emph{strongly advised}
  not to set the value of \code{asp} but rather to leave \code{asp}
  equal to its default value of \code{1}.  Any other value distorts
  the tesselation and destroys the perpendicular appearance of
  lines which are indeed perpendicular.  (And conversely can cause
  lines which are not perpendicular to appear as if they are.)
  The argument \code{asp} was added at the request of Zubin Dowlaty
  (who presumably knows what he's doing!).
  }
}
\section{Notes}{
  \itemize{
  \item If \code{clipp} is not \code{NULL} and \code{showpoints}
  is \code{TRUE} then it is possible that some of the points
  \dQuote{shown} will not fall inside any of the plotted tiles.
  (This will happen if the parts of the tiles in which they fall
  have been \dQuote{clipped} out.)  If a tile is clipped out
  \emph{completely} then the point which determines that tile is
  \emph{not} plotted irrespective of the value of \code{showpoints}.

  \item If the \code{z} components of the entries of \code{x}
  exist but cannot all be interpreted as colours then the internal
  function \code{getCol()} returns \code{NA}.  If \code{warn} is
  \code{TRUE} then a warning is issued.  The function \code{getCol()}
  will also return \code{NA} (no warning is issued in this case)
  if there \emph{aren't} any \code{z} components. This circumstance
  will arise if no \code{z} argument was supplied in the call to
  \code{deldir()}.  An \code{NA} value of \code{fillcol} results
  (as is indicated by the argument list entry for \code{fillcol})
  in (all of) the tiles being left unfilled.

  \item The change from argument \code{polycol} to argument
  \code{fillcol}, and the resulting change in the way in which
  plotted tiles are filled with colours, was made as a result of
  a request from Chris Triggs.  Likewise the argument \code{clipp}
  was added due to a request from Chris Triggs.
  }
}

\value{The list of tiles being plotted.  This will be the input
list of  tiles specified by argument \code{x}, or this list
clipped to the polygon \code{clipp} if the latter was specified.
}

\author{\rolf}
\seealso{ \code{\link{deldir}()}, \code{\link{tile.list}()},
          \code{\link{triang.list}()}, \code{\link{plot.triang.list}()} }
\examples{
  set.seed(42)
  x <- runif(20)
  y <- runif(20)
  z <- deldir(x,y,rw=c(0,1,0,1))
  w <- tile.list(z)
  plot(w)
  ccc <- heat.colors(20) # Or topo.colors(20), or terrain.colors(20)
                         # or cm.colors(20), or rainbow(20).
  plot(w,fillcol=ccc,close=TRUE)
  if(require(polyclip)) {
     CP  <- list(x=c(0.49,0.35,0.15,0.20,0.35,0.42,
                     0.43,0.62,0.46,0.63,0.82,0.79),
                 y=c(0.78,0.86,0.79,0.54,0.58,0.70,
                     0.51,0.46,0.31,0.20,0.37,0.54))
     cul <- rainbow(10)[c(1,7,3:6,2,8:10)] # Rearranging colours to improve
                                           # the contrast between contiguous tiles.
     plot(w,clipp=CP,showpoints=FALSE,fillcol=cul)
  }
  plot(w,labelPts=TRUE,col.lbls="red")
  plot(w,labelPts=TRUE,col.lbls="red",cex=0.5)
  plot(w,showpoints=TRUE,labelPts=TRUE,col.pts="green",col.lbls="red")
  plot(w,axes=FALSE,xlab="",ylab="",close=TRUE)
}
\keyword{ hplot }
