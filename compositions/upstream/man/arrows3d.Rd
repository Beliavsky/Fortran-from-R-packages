\name{arrows3D}
\alias{arrows3D}
\alias{arrows3D.default}
\title{arrows in 3D, based on package rgl}
\description{
  adds 3-dimensional arrows to an rgl plot.
}
\usage{
arrows3D(\dots)
\method{arrows3D}{default}(x0,x1,\dots,length=0.25,
                     angle=30,code=2,col="black",
                     lty=NULL,lwd=2,orth=c(1,0.0001,0.0000001),
                     labs=NULL,size=lwd)
}
\arguments{
  \item{x0}{a matrix or vector giving the starting points of the arrows}
  \item{x1}{a matrix or vector giving the end points of the arrows}
  \item{\dots}{additional plotting parameters as described in \code{rgl::material3d}
       %\code{\link[rgl]{rgl.material}}
   }
  \item{length}{a number giving the length of the arrowhead}
  \item{angle}{numeric giving the angle of the arrowhead}
  \item{code}{0=no arrowhead,1=arrowhead at x0,2=arrowhead at
    x1,3=double headed}
  \item{col}{the color of the arrow}
  \item{lty}{Not implemented, here for compatibility reasons with
    \code{arrows}}
  \item{lwd}{line width in pixels}
  \item{orth}{the flat side of the arrow is not unique by x0 and
    x1. This ambiguity is solved in a way that the arrow seams as wide
    as possible from the viewing direction orth.}
  \item{labs}{labels to be plotted to the endpoints of the arrows}
  \item{size}{size of the plotting symbol}
}
\value{
    the 3D plotting coordinates of the tips of the arrows displayed, 
    returned invisibly
}
\details{
  The function is called to plot arrows into an rgl plot. The size of
  the arrow head is given in a absolute way. Therefore it is
  important to give the right scale for the length, to see the arrow head and
  that it does not fill the whole window. 
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  %\code{\link[rgl]{points3d}},
  %\code{\link[graphics]{plot}}, 
  \code{\link{plot3D}},
  \code{rgl::points3d}, \code{graphics::plot}
}
\examples{
x <- cbind(rnorm(10),rnorm(10),rnorm(10))
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(x)
  x0 <- x*0
  arrows3D(x0,x)
} ## this function requires package 'rgl'
}
\keyword{multivariate}
