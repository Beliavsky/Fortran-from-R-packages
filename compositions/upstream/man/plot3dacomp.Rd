\name{plot3Dacomp}
\alias{plot3D.acomp}
\alias{plot3D.rcomp}
\title{3D-plot of compositional data}
\description{
  3D-plot of compositional data. The plot is
  mainly an exploratory tool, not intended for exact display of data.
}
\usage{
 \method{plot3D}{acomp}(x, parts=1:min(ncol(X),4),..., 
             lwd=2, axis.col="gray", add=FALSE, cex=2, 
             vlabs=colnames(x), vlabs.col=axis.col, center=FALSE, 
             scale=FALSE, log=FALSE, bbox=FALSE, axes=TRUE, size=cex,col=1)
\method{plot3D}{rcomp}(x,parts=1:min(ncol(X),4),...,
             lwd=2,axis.col="gray",add=FALSE,cex=2,
             vlabs=colnames(x),vlabs.col=axis.col,center=FALSE,
             scale=FALSE,log=FALSE,bbox=FALSE,axes=TRUE,size=cex,col=1)
}
\arguments{
  \item{x}{an aplus object to be plotted}
  \item{parts}{a numeric xor character vector of length 3 coding the
    columns to be plotted}
  \item{\dots}{additional plotting parameters as described in
    %\code{\link[rgl]{rgl.material}}
    \code{rgl::material3d}
    }
  \item{add}{logical, adding  or new plot}
  \item{cex}{size of the plotting symbols}
  \item{lwd}{line width}
  \item{axis.col}{color of the axis}
  \item{vlabs}{the column names to be plotted, if missing defaults to
    the column names of the selected columns of X}
  \item{vlabs.col}{color of the labels}
  \item{center}{logical, should the data be centered}
  \item{scale}{logical, should the data be scaled}
  \item{log}{logical, indicating wether to plot in log scale}
  \item{bbox}{logical, whether to add a bounding box}
  \item{axes}{logical, whether plot a coordinate cross}
  \item{size}{size of the plotting symbols}
  \item{col}{the color used for dots, defaults to black.}

}
\value{
  Called for its side effect of a 3D plot of an acomp object in an rgl plot. 
  It invisibly returns the 3D plotting coordinates of the objects displayed

}
\details{
  The routine behaves different when 3 or four components should be
  plotted. In case of four components:
  \cr
  If log is TRUE the data is plotted in \code{\link{ilr}}
  coordinates. This is the isometric view of the data.
  \cr
  If log is FALSE the data is plotted in \code{\link{ipt}} coordinates
  and a tetrahedron is plotted around it if \code{coors == TRUE}. This
  can be used to do a tetrahedron plot.
  \cr
  In case of three components:
  \cr
  If log is TRUE the data is plotted in \code{\link{clr}}
  coordinates. This can be used to visualize the clr plane.
  \cr
  If log is FALSE the data is plotted as is, showing the embedding of
  the
  three-part simplex in the three-dimensional space. 
  \cr
  In all cases:
  If \code{coors} is true, coordinate arrows are plotted
  of length 1 in the origin of the space, except in the tetrahedron case.
}
\note{
The function \code{\link{kingTetrahedron}} provides an alternate way of
tetrahedron plots, based on a more advanced viewer, which must
be downloaded separately.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{kingTetrahedron}}
    %\code{\link[rgl]{points3d}},
   %\code{\link[graphics]{plot}}, 
  \code{rgl::points3d},
  \code{graphics::plot},
  \code{\link{plot3D}},
  \code{\link{plot3D.rmult}},
 
  \code{\link{plot3D.rcomp}},
  \code{\link{plot3D.aplus}},\code{\link{plot3D.rplus}}
}
\examples{
data(SimulatedAmounts)
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(acomp(sa.lognormals5),1:3,col="green")
  plot3D(acomp(sa.lognormals5),1:3,log=TRUE,col="green")
  plot3D(acomp(sa.lognormals5),1:4,col="green")
  plot3D(acomp(sa.lognormals5),1:4,log=TRUE,col="green")
} ## this function requires package 'rgl'
}
\keyword{multivariate}
