\name{plot3Daplus}
\alias{plot3D.aplus}
\title{3D-plot of positive data}
\description{
  3D-plot of positive data typically in log-log-log scale. The plot is
  mainly an exploratory tool, and not intended for exact display of data.
}
\usage{
\method{plot3D}{aplus}(x,parts=1:3,...,
                vlabs=NULL,add=FALSE,log=TRUE,bbox=FALSE,axes=TRUE,col=1)
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
  \item{vlabs}{the column names to be plotted, if missing defaults to
    the column names of the selected columns of X}
  \item{log}{logical, indicating wether to plot in log scale}
  \item{bbox}{logical, whether to add a bounding box}
  \item{axes}{logical, plot a coordinate system}
  \item{col}{the color used for dots, defaults to black.}
}
\value{
  Called for its side effect of a 3D plot of an aplus object in an rgl plot. 
  It invisibly returns the 3D plotting coordinates of the objects displayed
}
\details{
  If log is TRUE the data is plotted in \code{\link{ilt}}
  coordinates. If \code{coors} is true, coordinate arrows are plotted
  of length 1 and in the (aplus-)mean of the dataset.
  \cr
  If log is FALSE the data is plotted with plot.rplus
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
  \code{\link{plot3D.rplus}}
}
\examples{
data(SimulatedAmounts)
if(requireNamespace("rgl", quietly = TRUE)) {
  plot3D(aplus(sa.lognormals),size=2)
} ## this function requires package 'rgl'

}
\keyword{multivariate}
