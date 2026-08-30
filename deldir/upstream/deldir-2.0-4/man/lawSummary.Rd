\name{lawSummary}
\alias{lawSummary}
\title{
  Produce a Lewis-Aboav-Weaire summary of a tessellation.
}
\description{
  Produce a summary of a Dirichlet (Voronoi) tessellation in terms
  of parameters relevant to Lewis's law and Aboav-Weaire's law.
  Note that \dQuote{law} in the function name corresponds to
  \dQuote{Lewis-Aboav-Weaire}.
}

\usage{
lawSummary(object)
}
\arguments{
  \item{object}{
  An object of class \code{"deldir"} as returned by the function
  \code{deldir()}.
}
}
\details{
  Tiles are stripped away from the tessellation in \dQuote{layers}.
  Layer 1 consists of \dQuote{boundary} tiles, i.e. tiles having
  at least one vertex on the enclosing rectangle (determined by the
  \code{rw} argument of \code{\link{deldir}()}).  Layer 2 consists
  of tiles which are neighbours of tiles in layer 1 (i.e. tiles
  determined by points that are Delaunay neighbours of points
  determining the tiles in layer 1).  Layer 3 consists of tiles
  which are neighbours of tiles in layer 2.

  The parameters of interest in respect of the Lewis-Aboav-Weaire summary are
  then calculated in terms of the tiles that remain after the three layers
  have been stripped away, which will be referred to as \dQuote{interior}
  tiles.  These parameters are:
  \itemize{
    \item the areas of each of the interior tiles
    \item the number of edges of each of the interior tiles
    \item the number of edges of all neighbouring tiles of each of
    the interior tiles.
  }
  Note that the neighbouring tiles of the interior tiles may include
  tiles which are \emph{not themselves} interior tiles (i.e. tiles
  which are in layer 3).

  This function was created at the request of Kai Xu
  (Fisheries College, Jimei University, Xiamen, Fujian, China 361021).
}
\value{
  If no tiles remain after the three layers have been stripped
  away, then the returned value is \code{NULL}.  Otherwise the
  returned value is a list with components calculated in terms of
  the remaining (\dQuote{interior}) tiles.  These components are:
  \itemize{
     \item \code{tile.vertices} A list whose entries are data frames
     giving the coordinates of the vertices of the interior tiles.
     \item \code{tile.areas} A vector of the areas of the interior tiles
     in the tessellation in question.
     \item{tile.tags} A vector or factor whose values are the \dQuote{tags}
     of the interior tiles.  The \dQuote{original} of this object (the
     \dQuote{tags} associated with all of the tiles) is provided
     as the \code{z} argument to \code{deldir()}.  The \code{tile.tags}
     component of the value returned by \code{lawSummary()} is present
     only if \code{deldir()} was called with a (non-\code{NULL}) value
     of the \code{z} argument.
     \item \code{num.edges} A vector of the number of edges of each
     such tile.
     \item \code{num.nbr.edges} A list with a component for each
     point, in the set being tessellated, whose corresponding tile
     is an interior tile.  Each component of this
     list is the vector of the number of edges of the interior tiles
     determined by points which are Delaunay neighbours of
     the point corresponding to the list component in question.
     \item \code{totnum.nbr.edges} A vector whose entries consist
     of the sums of the vectors in the foregoing list.
  }
  The returned list also has attributes as follows:
  \itemize{
     \item \code{i1} An integer vector whose entries are in the indices
     of the tiles in layer 1.
     \item \code{i2} An integer vector whose entries are in the indices
     of the tiles in layer 2.
     \item \code{i3} An integer vector whose entries are in the indices
     of the tiles in layer 3.
     \item \code{i.kept} An integer vector whose entries are in the indices
     of the tiles that are kept, i.e. those that remain after the three layers
     have been stripped away.
  }
}

\author{\rolf}

\seealso{ \code{\link{tile.list}()} \code{\link{tile.centroids}()} }

\examples{
# A random pattern:
set.seed(42)
xy1   <- data.frame(x=runif(400,0,20),y=runif(400,0,20))
dxy1  <- deldir(xy1)
ldxy1 <- lawSummary(dxy1)
tl1   <- tile.list(dxy1)
plot(0,0,type="n",xlim=c(-2,35),ylim=c(0,20),asp=1,xlab="x",ylab="y",bty="l")
plot(tl1,showpoints=FALSE,add=TRUE)
points(xy1[attr(ldxy1,"i1"),],pch=20,col="yellow")
points(xy1[attr(ldxy1,"i2"),],pch=20,col="blue")
points(xy1[attr(ldxy1,"i3"),],pch=20,col="green")
points(xy1[attr(ldxy1,"i.kept"),],pch=20,col="red")
legend("right",pch=20,col=c("yellow","blue","green","red"),
       legend=c("layer 1","layer 2","layer 3","interior"))

# A highly structured pattern (example due to Kai Xu):
set.seed(115)
x  <- c(rep(1:20,10),rep((1:20)+0.5,10))
y  <- c(rep(1:10,each=20),rep((1:10)+0.5,each=20))*sqrt(3)
a  <- runif(400,0,2*pi)
b  <- runif(400,-1,1)
x  <- x+0.1*cos(a)*b
y  <- y+0.1*sin(a)*b
xy2 <- data.frame(x,y)
dxy2  <- deldir(xy2)
ldxy2 <- lawSummary(dxy2)
tl2   <- tile.list(dxy2)
plot(0,0,type="n",xlim=c(-2,35),ylim=c(0,20),asp=1,xlab="x",ylab="y",bty="l")
plot(tl2,showpoints=FALSE,add=TRUE)
points(xy2[attr(ldxy2,"i1"),],pch=20,col="yellow")
points(xy2[attr(ldxy2,"i2"),],pch=20,col="blue")
points(xy2[attr(ldxy2,"i3"),],pch=20,col="green")
points(xy2[attr(ldxy2,"i.kept"),],pch=20,col="red")
legend("right",pch=20,col=c("yellow","blue","green","red"),
       legend=c("layer 1","layer 2","layer 3","interior"))

}
\keyword{spatial}
