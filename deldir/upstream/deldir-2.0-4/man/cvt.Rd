\encoding{UTF-8}
\name{cvt}
\alias{cvt}
\title{
    Centroidal Voronoi (Dirichlet) tessellation.
}
\description{
    Calculates the centroidal Voronoi (Dirichlet) tessellation using
    Lloyd's algorithm.
}
\usage{
    cvt(object, stopcrit = c("change", "maxit"), tol = NULL,
       maxit = 100, verbose = FALSE)
}
\arguments{
  \item{object}{
  An object of class either \code{"deldir"} (as returned by \code{\link{deldir}()})
  or \code{"tile.list"} (as returned by \code{\link{tile.list}()}).
}
  \item{stopcrit}{
  Text string specifying the stopping criterion for the algorithm. If
  this is \code{"change"} then the algorithm halts when the maximum
  change in in the distances between corresponding centroids,
  between iterations, is less than \code{tol} (see below).
  It \code{stopcrit} is \code{"maxit"} then the algorithm halts
  after a specified number of iterations (\code{maxit}; see below)
  have been completed.  This argument may be abbreviated, e.g. to
  \code{"c"} or \code{"m"}.
}
  \item{tol}{
  The tolerance used when the stopping criterion is \code{"change"}.
  Defaults to \code{.Machine$double.eps}.
}
  \item{maxit}{
  The maximum number of iterations to perform when the stopping criterion
  is \code{"maxit"}.
}
  \item{verbose}{
  Logical scalar.  If \code{verbose} is \code{TRUE} then rudimentary
  \dQuote{progress reports} are printed out, every 10 iterations if
  the stopping criterion is \code{"change"} or every iteration if the
  stopping criterion is \code{"maxit"}.
}
}
\section{Note}{
This function was added to the \code{deldir} package at the
suggestion of Dr. Michaël Aupetit, Senior Scientist at the
Qatar Computing Research Institute, Hamad Bin Khalifa University.
}
\details{The algorithm iteratively tessellates a set of points and
  then replaces the points with the centroids of the tiles associated
  with those points.  \dQuote{Eventually} (at convergence) points
  and the centroids of their associated tiles coincide.
}
\value{
  A list with components:
  \item{centroids}{A data frame with columns \code{"x"} and
  \code{"y"} specifying the coordinates of the limiting locations
  of the tile centroids.}
  \item{tiles}{An object of class \code{"tile.list"} specifying
  the Dirichlet (Voronoi) tiles in the tessellation of the points
  whose coordinates are given in \code{centroids}.  Note the tile
  associated with the \eqn{i}th point has centroid \emph{equal}
  to that point.}
}
\references{
\url{https://en.wikipedia.org/wiki/Lloyd's_algorithm}

Lloyd, Stuart P. (1982). Least squares quantization in PCM.
\emph{IEEE Transactions on Information Theory} \bold{28} (2),
pp. 129--137, doi:10.1109/TIT.1982.1056489.
}

\author{\rolf}
\seealso{
\code{\link{deldir}()} \code{\link{tile.list}()}
}
\examples{
\dontrun{ # Takes too long.
    set.seed(42)
    x <- runif(20)
    y <- runif(20)
    dxy <- deldir(x,y,rw=c(0,1,0,1))
    cxy1 <- cvt(dxy,verb=TRUE)
    plot(cxy1$tiles)
    with(cxy1$centroids,points(x,y,pch=20,col="red"))
    cxy2 <- cvt(dxy,stopcrit="m",verb=TRUE)
    plot(cxy2$tiles)
    with(cxy2$centroids,points(x,y,pch=20,col="red"))
# Visually indistinguishable from the cxy1 plot.
# But ...
    all.equal(cxy1$centroids,cxy2$centroids) # Not quite.
    cxy3 <- cvt(dxy,stopcrit="m",verb=TRUE,maxit=250)
    all.equal(cxy1$centroids,cxy3$centroids) # Close, but no cigar.
    cxy4 <- cvt(dxy,verb=TRUE,tol=1e-14)
    cxy5 <- cvt(dxy,stopcrit="m",verb=TRUE,maxit=600)
    all.equal(cxy4$centroids,cxy5$centroids) # TRUE
# Takes a lot of iterations or a really small tolerance
# to get "good" convergence.  But this is almost surely
# of no practical importance.
    txy <- tile.list(dxy)
    cxy6 <- cvt(txy)
    all.equal(cxy6$centroids,cxy1$centroids) # TRUE
}
}
\keyword{spatial}
