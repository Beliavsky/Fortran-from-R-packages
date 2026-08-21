\name{axis3D}
\alias{axis3D}
\title{Drawing a 3D coordiante system to a plot, based on package rgl}
\description{
  Adds a coordinate system to a 3D rgl graphic. In future releases,
  functionality to add tickmarks will be (hopefully) provided. 
  Now, it is just a system of arrows giving the
  directions of the three axes.
}
\usage{
axis3D(axis.origin=c(0,0,0),axis.scale=1,axis.col="gray",vlabs=c("x","y","z"),
       vlabs.col=axis.col,bbox=FALSE,axis.lwd=2,axis.len=mean(axis.scale)/10,
       axis.angle=30,orth=c(1,0.0001,0.000001),axes=TRUE,...)
}
\arguments{
  \item{axis.origin}{The location where to put the origin of the
    coordinate arrows typicall either 0, the minimum or the mean of the
    dataset}
  \item{axis.scale}{either a number or a 3D vector giving the length of
    the arrows for the axis in the coordiantes of the plot}
  \item{axis.col}{Color to plot the coordinate system}
  \item{vlabs}{The names of the axes, plotted at the end}
  \item{vlabs.col}{color for the axes labels}
  \item{bbox}{boolean, whether to plot a bounding box}
  \item{axis.lwd}{line width of the axes}
  \item{axis.angle}{angle of the arrow heads}
  \item{axis.len}{length of the arrow heads}
  \item{orth}{the orth argument of \code{\link{arrows3D}}  }
  \item{axes}{a boolean, wether to plot the axes}
  \item{\dots}{these arguments are passed to arrows3D as
    %\code{\link[rgl]{rgl.material}} 
    \code{rgl::material3d} arguments}
}
\value{
  Nothing
}
\details{
  The function is called to plot a coordiante system consisting of
  arrows into an rgl plot.  
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  %\code{\link[rgl]{points3d}},
  %\code{\link[graphics]{plot}}, 
  \code{rgl::points3d}, \code{graphics::plot},
  \code{\link{plot3D}},\code{\link{arrows3D}}
}
\examples{
x <- cbind(rnorm(10),rnorm(10),rnorm(10))
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(x)
  x0 <- x*0
  axis3D()
} ## this function requires package 'rgl'
}
\keyword{multivariate}
