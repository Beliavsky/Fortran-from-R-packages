\name{plot3Drplus}
\alias{plot3D.rplus}
\title{plot in 3D based on rgl}
\description{
  3-dimensional plots, which can be rotated and zoomed in/out
}
\usage{
\method{plot3D}{rplus}(x,parts=1:3,...,vlabs=NULL,add=FALSE,bbox=FALSE,
                           cex=1,size=cex,axes=TRUE,col=1)
}
\arguments{
  \item{x}{an rplus object to be plotted}
  \item{parts}{the variables in the rplus object to be plotted}
  \item{\dots}{additional plotting parameters as described in
        %\code{\link[rgl]{rgl.material}}
    \code{rgl::material3d}
    }
  \item{vlabs}{the labels used for the variable axes}
  \item{add}{logical, adding  or new plot}
  \item{bbox}{logical, whether to add a bounding box}
  \item{cex}{size of the plotting symbol (as character expansion factor)}
  \item{size}{size of the plotting symbol, only size or cex should be used}
  \item{axes}{logical, whether to plot a coordinate cross}
  \item{col}{the color used for dots, defaults to black.}

}
\value{
  the 3D plotting coordinates of the objects displayed, returned invisibly
}
\details{
  The function plots rplus objects in a 3D coordinate system, in an rgl plot. 
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

  
  \code{\link{plot3D.acomp}},\code{\link{plot3D.rcomp}},
  \code{\link{plot3D.aplus}},\code{\link{plot3D}}
}
\examples{
x <- cbind(rnorm(10),rnorm(10),rnorm(10))
data(SimulatedAmounts)
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(rplus(exp(x)))
  plot3D(rplus(sa.lognormals),cex=4,col=1:nrow(sa.lognormals))
} ## this function requires package 'rgl'
}
\keyword{multivariate}
