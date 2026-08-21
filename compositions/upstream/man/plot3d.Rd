\name{plot3D}
\alias{plot3D}
\alias{plot3D.default}
\title{plot in 3D based on rgl}
\description{
  3-dimensional plots, which can be rotated and zoomed in/out
}
\usage{
plot3D(x,...)
\method{plot3D}{default}(x,\dots,add=FALSE,bbox=TRUE,axes=FALSE,
           cex=1,size=cex,col=1)
}
\arguments{
  \item{x}{an object to be plotted, e.g. a data frame or a data matrix}
  \item{\dots}{additional plotting parameters as described in
  %\code{\link[rgl]{rgl.material}}
  \code{rgl::material3d}  
  }
  \item{add}{logical, adding  or new plot}
  \item{bbox}{logical, whether to add a bounding box}
  \item{axes}{logical, whether to plot an axes of coordinates}
  \item{cex}{size of the plotting symbol}
  \item{size}{size of the plotting symbol, only size or cex should be
    used}
  \item{col}{the color used for dots, defaults to black.}
}
\value{
  the 3D plotting coordinates of the objects displayed, returned invisibly
}
\details{
  The function provides a generic interface for 3-dimensional plotting
  in analogy to the 2d-plotting interface of plot, using rgl package.
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  %\code{\link[rgl]{points3d}},
  %\code{\link[graphics]{plot}}, 
  \code{rgl::points3d},
  \code{graphics::plot},
  \code{\link{plot3D.rmult}},
  
  \code{\link{plot3D.acomp}},\code{\link{plot3D.rcomp}},
  \code{\link{plot3D.aplus}},\code{\link{plot3D.rplus}}
}
\examples{
x <- cbind(rnorm(10),rnorm(10),rnorm(10))
data(SimulatedAmounts)
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(x)
  plot3D(sa.lognormals,cex=4,col=1:nrow(sa.lognormals))
} ## this function requires package 'rgl'
}
\keyword{multivariate}
