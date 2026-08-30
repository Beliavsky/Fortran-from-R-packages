\name{seaweed}
\alias{seaweed}
\docType{data}
\title{
    seaweed
}
\description{
   A data frame whose columns are the coordinates of the centroids
   of the cells in a seaweed frond. The points are estimates-by-eye
   of where the centroids of the cells occur.
}

\usage{data("seaweed")}

\format{
  A data frame with 266 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{The \eqn{x}-coordinates of the cell centroids.}
    \item{\code{y}}{The \eqn{y}-coordinates of the cell centroids.}
  }
}

\source{
  These data were kindly supplied by Dr. John Bothwell of the
  Department of Biosciences, Durham University.  The data were
  collected by Kevin Yun and Georgia Campbell, members of Dr.
  Bothwell's research group.
}
\examples{
# data(seaweed)
# It is unnecessary to use \code{data} since \code{seaweed} is
# a "first class object".  It is "lazily loaded".

dsw <- deldir(seaweed)
isw <- tileInfo(dsw)
# Expand the width of the terminal window.
options(width=120)
isw
tsw <- tile.list(dsw)
plot(tsw,labelPts=TRUE,col.lbls="red",cex=0.5,adj=0.5)
}

\keyword{datasets}
