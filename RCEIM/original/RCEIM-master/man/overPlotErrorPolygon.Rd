\name{overPlotErrorPolygon}
\alias{overPlotErrorPolygon}

\title{
Overplot an error polygon
}
\description{
A simple function to overplot an error polygon around a curve. Note that the error is considered as symmetric, and exclusively on y. The polygon will be created from the coordinate tuples \eqn{(x, (y - err_y))} and \eqn{(x, (y + err_y))}.
}
\usage{
overPlotErrorPolygon(x, y, err_y, col = "grey", logPlot = FALSE, ...)
}
\arguments{
  \item{x}{
  A vector containing the x coordinate of the data.
}
  \item{y}{
  A vector containing the y coordinate of the data.
}
  \item{err_y}{
  A vector containing the error in y.
}
  \item{col}{
  The color that will be used for filling the polygon.
}
  \item{logPlot}{
  A boolean indicating if the plot is in logscale.
}
  \item{\dots}{
  Further arguments to be passed to polygon().
}
}

\value{
A polygon is overplotted in the active graphics device.
}
\author{
Alberto Krone-Martins
}

\examples{
# Shows a simple random curve and overplots a randomly created error bar.
dev.new()
xx <- 1:10
yy <- (1:10)/5 + 4 + (runif(10)-0.5)/0.5*2
plot(xx, yy, type="l", xlab="x", ylab="y", ylim=c(0,10))
err_yy <- 1.5 + (runif(10)-0.5)/0.5
overPlotErrorPolygon(xx,yy,err_yy, col=rgb(0,0,1,0.3), border=NA)
}

\keyword{ aplot }