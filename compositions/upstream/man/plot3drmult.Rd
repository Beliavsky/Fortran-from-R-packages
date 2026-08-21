\name{plot3Drmult}
\alias{plot3D.rmult}
\title{plot in 3D based on rgl}
\description{
  3-dimensional plots, which can be rotated and zoomed in/out
}
\usage{
\method{plot3D}{rmult}(x,parts=1:3,...,
                 center=FALSE,scale=FALSE,add=FALSE,axes=!add,
                 cex=2,vlabs=colnames(x),size=cex,bbox=FALSE,col=1)
}
\arguments{
  \item{x}{an object to be plotted, e.g. a data frame or a data matrix}
  \item{parts}{the variables in the rmult object to be plotted}
  \item{\dots}{additional plotting parameters as described in
        %\code{\link[rgl]{rgl.material}}
    \code{rgl::material3d}
    }
  \item{center}{logical, center the data? This might be necessary to
    stay within the openGL-arithmetic used in rgl.}
  \item{scale}{logical, scale the data? This might be necessary to
    stay within the openGL-arithmetic used in rgl.}
  \item{add}{logical, adding  or new plot}
  \item{bbox}{logical, whether to add a bounding box}
  \item{axes}{logical, whether to plot a coordinate cross}
  \item{cex}{size of the plotting symbol (as expanding factor)}
  \item{vlabs}{labels for the variables}
  \item{size}{size of the plotting symbol, only size or cex should be used}
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
  \code{\link{kingTetrahedron}}
    %\code{\link[rgl]{points3d}},
   %\code{\link[graphics]{plot}}, 
  \code{rgl::points3d},
  \code{graphics::plot},
  \code{\link{plot3D}},
  
  \code{\link{plot3D.acomp}},\code{\link{plot3D.rcomp}},
  \code{\link{plot3D.aplus}},\code{\link{plot3D.rplus}}
}
\examples{
x <- cbind(rnorm(10),rnorm(10),rnorm(10))
data(SimulatedAmounts)
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(x)
  plot3D(rmult(sa.lognormals),cex=4,col=1:nrow(sa.lognormals))
} ## this function requires package 'rgl'
}
\keyword{multivariate}
