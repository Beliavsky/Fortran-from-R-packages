\name{plot.divchain}
\alias{plot.divchain}
\title{
   Plot a dividing chain.
}
\description{
   Plot the dividing chain of a Dirichlet tessellation.
   The tessellation must have been created from a set of points
   having associated categorical \dQuote{tags}.  The dividing
   chain consists of those edges of Dirichlet tiles which separate
   points having different values of the given tags.
}
\usage{
   \method{plot}{divchain}(x, add = FALSE, ...)
}
\arguments{
  \item{x}{
     An object of class \dQuote{divchain}.  See \code{\link{divchain.deldir}()}
     for details.
}
  \item{add}{
     Logical scalar.  It \code{add=TRUE} the plot of the dividing chain
     is added to an existing plot.
}
  \item{\dots}{
     Graphical parameters such as \code{main}, \code{xlab}, \code{col.main},
     \code{col.lab}.  In particular if \code{bty} is supplied
     (as a value other than \code{n}) a \dQuote{box} will be drawn
     around the plot that is formed when \code{add=FALSE}.
     Also a non-standard graphical parameter \code{boxcol} may be supplied
     which will be taken to be the colour with which the box is drawn.
     If a \code{col} argument is supplied, this determines the colour
     for plotting the segments constituting the dividing chain.
}
}
\value{
   None.
}
\author{\rolf}
\note{
  This function was created in response to a question asked
  on \code{stackoverflow.com} by a user named \dQuote{Dan}.
}

\seealso{
  \code{\link{divchain}()} \code{\link{divchain.default}()}
  \code{\link{divchain.deldir}()} \code{\link{deldir}()}
}
\examples{
   set.seed(42)
   x   <- runif(50)
   y   <- runif(50)
   z   <- factor(kmeans(cbind(x,y),centers=4)$cluster)
   dc  <- divchain(x,y,z,rw=c(0,1,0,1))
   plot(dc,lwd=2,col="blue",bty="o")
}
\keyword{ spatial }
