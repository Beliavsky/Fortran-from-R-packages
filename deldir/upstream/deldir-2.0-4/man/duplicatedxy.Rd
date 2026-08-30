\name{duplicatedxy}
\alias{duplicatedxy}
\title{
   Determine duplicated points.
}
\description{
    Find which points among a given set are duplicates of others.
}
\usage{
duplicatedxy(x, y)
}
\arguments{
  \item{x}{
  Either a vector of \code{x} coordinates of a set of (two
  dimensional) points, or a list (or data frame) with columns
  \code{x} and \code{y} giving the coordinates of a set of such
  points.
}
  \item{y}{
  A vector of \code{y} coordinates of a set of (two
  dimensional) points.  Ignored if \code{x} is a list or data
  frame.
}
}
\details{
   Often it is of interest to associate each Dirichlet tile in a
   tessellation of a planar point set with the point determining
   the tile.  This becomes problematic if there are \emph{duplicate}
   points in the set being tessellated/triangulated.  Duplicated
   points are automatically eliminated \dQuote{internally} by
   \code{deldir()}.  The association between tiles and the indices
   of the original set of points is now preserved by the component
   \code{ind.orig} of the object returned by \code{deldir()}.
   However confusion could still arise.

   If it is of interest to associate Dirichlet tiles with the
   points determining them, then it is better to proceed by
   eliminating duplicate points to start with.  This function
   (\code{duplicatedxy()}) provides a convenient way of doing so.
}
\value{
   A logical vector of length equal to the (original) number
   of points being considered, with entries \code{TRUE} if the
   corresponding point is a duplicate of a point with a smaller
   index, and \code{FALSE} otherwise.
}
\author{\rolf}
\note{
  The real work is done by the base \bold{R} function \code{\link{duplicated}()}.
}

\section{Warning}{
   Which indices will be considered to be indices of duplicated
   points (i.e. get \code{TRUE} values) will of course depend on
   the order in which the points are presented.
}

\seealso{
\code{\link{duplicated}()}, \code{\link{deldir}()}
}
\examples{
   set.seed(42)
   xy  <- data.frame(x=runif(20),y=runif(20))
   # Lots of duplicated points.
   xy  <- rbind(xy,xy[sample(1:20,20,TRUE),])
   # Scramble.
   ii  <- sample(1:40,40)
   x   <- xy$x[ii]
   y   <- xy$y[ii]
   # Unduplicate!
   iii <- !duplicatedxy(x,y)
   xu  <- x[iii]
   yu  <- y[iii]
   # The i-th tile is determined by (xu[i],yu[i]):
   dxy <- deldir(xu,yu)
}
\keyword{ utilities }
