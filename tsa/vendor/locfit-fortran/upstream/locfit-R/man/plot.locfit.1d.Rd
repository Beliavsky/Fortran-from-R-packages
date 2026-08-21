\name{plot.locfit.1d}
\alias{plot.locfit.1d}
\title{
Plot a one dimensional preplot.locfit object.
}
\usage{
\method{plot}{locfit.1d}(x, add=FALSE, main="", xlab="default", ylab=x$yname,
  type="l", ylim, lty=1, col=1, \dots)
}
\arguments{
  \item{x}{One dimensional \code{preplot.locfit} object.}
  \item{add}{If \code{TRUE}, the plot will be added to the existing
    plot.}
  \item{main, xlab, ylab, type, ylim, lty, col}{Graphical parameters
    passed on to \code{\link{plot}} (only if \code{add=FALSE}).}
  \item{...}{Additional graphical parameters to the \code{plot} function
    (only if \code{add=FALSE}).}
}
\description{
This function is not usually called directly. It will be called automatically
when plotting a one-dimensional \code{locfit} or \code{preplot.locfit}
object.
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}, \code{\link{preplot.locfit}}
}
\keyword{methods}
